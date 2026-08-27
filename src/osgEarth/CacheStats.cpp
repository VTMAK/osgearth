/* osgEarth
 * Copyright 2025 Pelican Mapping
 * MIT License
 */
#include <osgEarth/CacheStats>
#include <osgEarth/Cache>
#include <osgEarth/DateTime>
#include <osgEarth/FileUtils>
#include <osgEarth/Notify>
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mutex>
#include <sstream>

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

using namespace osgEarth;

namespace
{
    constexpr std::size_t HISTOGRAM_BUCKETS = 64u;

    const char* metricName(CacheStatistics::Metric metric)
    {
        static const char* names[] = {
            "read", "write", "backend_read", "backend_write", "serialize",
            "deserialize", "reader_wait", "transaction", "remove", "touch", "clear"
        };
        return names[static_cast<std::size_t>(metric)];
    }

    std::string jsonEscape(const std::string& input)
    {
        std::ostringstream out;
        for (unsigned char c : input)
        {
            switch (c)
            {
            case '\\': out << "\\\\"; break;
            case '"':  out << "\\\""; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default:
                if (c < 0x20u)
                    out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<unsigned>(c) << std::dec;
                else
                    out << static_cast<char>(c);
            }
        }
        return out.str();
    }

    std::uint64_t processID()
    {
#ifdef _WIN32
        return static_cast<std::uint64_t>(::GetCurrentProcessId());
#else
        return static_cast<std::uint64_t>(::getpid());
#endif
    }

    bool parseBool(const std::string& value, bool& output)
    {
        std::string lower = value;
        std::transform(lower.begin(), lower.end(), lower.begin(),
            [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (lower == "1" || lower == "true" || lower == "yes" || lower == "on" || lower == "stdout")
        {
            output = true;
            return true;
        }
        if (lower == "0" || lower == "false" || lower == "no" || lower == "off")
        {
            output = false;
            return true;
        }
        return false;
    }

    void updateMaximum(std::atomic<std::uint64_t>& target, std::uint64_t value)
    {
        std::uint64_t current = target.load(std::memory_order_relaxed);
        while (current < value &&
            !target.compare_exchange_weak(current, value, std::memory_order_relaxed)) { }
    }

    std::string replaceAll(std::string value, const std::string& token, const std::string& replacement)
    {
        for (std::size_t pos = 0u; (pos = value.find(token, pos)) != std::string::npos; )
        {
            value.replace(pos, token.size(), replacement);
            pos += replacement.size();
        }
        return value;
    }
}

struct CacheStatistics::Impl
{
    struct MetricData
    {
        std::atomic<std::uint64_t> count{ 0u };
        std::atomic<std::uint64_t> success{ 0u };
        std::atomic<std::uint64_t> miss{ 0u };
        std::atomic<std::uint64_t> error{ 0u };
        std::atomic<std::uint64_t> bytes{ 0u };
        std::atomic<std::uint64_t> totalNanos{ 0u };
        std::atomic<std::uint64_t> maxNanos{ 0u };
        std::array<std::atomic<std::uint64_t>, HISTOGRAM_BUCKETS> histogram;

        MetricData()
        {
            for (auto& bucket : histogram)
                bucket.store(0u, std::memory_order_relaxed);
        }
    };

    std::string driver;
    std::string cachePath;
    std::string outputPath;
    unsigned intervalSeconds = 0u;
    std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();
    std::atomic<std::int64_t> lastReportSeconds{ 0 };
    std::array<MetricData, static_cast<std::size_t>(Metric::Count)> metrics;
    std::atomic<std::uint64_t> events{ 0u };
    std::atomic<std::uint64_t> busyRetries{ 0u };
    std::atomic<std::uint64_t> transactionMutations{ 0u };
    std::atomic<std::uint64_t> queueHighWaterBytes{ 0u };
    std::atomic<std::uint64_t> queueHighWaterItems{ 0u };
    std::atomic_bool finalReported{ false };
    mutable std::mutex reportMutex;
};

std::shared_ptr<CacheStatistics> CacheStatistics::create(
    const CacheOptions& options,
    const std::string& driver,
    const std::string& cachePath)
{
    bool enabled = options.collectStats().get();
    std::string outputPath = options.statsPath().get();
    if (const char* env = std::getenv(OSGEARTH_ENV_CACHE_STATS))
    {
        const std::string value(env);
        bool parsed = false;
        if (parseBool(value, parsed))
        {
            enabled = parsed;
        }
        else if (!value.empty())
        {
            enabled = true;
            outputPath = value;
        }
    }

    if (!enabled)
        return nullptr;

    auto impl = std::make_unique<Impl>();
    impl->driver = driver;
    impl->cachePath = cachePath;
    impl->outputPath = outputPath;
    impl->intervalSeconds = options.statsInterval().get();
    impl->lastReportSeconds.store(0, std::memory_order_relaxed);
    return std::shared_ptr<CacheStatistics>(new CacheStatistics(std::move(impl)));
}

CacheStatistics::CacheStatistics(std::unique_ptr<Impl> impl) : _impl(std::move(impl)) { }

CacheStatistics::~CacheStatistics()
{
    report(true);
}

void CacheStatistics::record(Metric metric, Outcome outcome, std::uint64_t nanos, std::uint64_t bytes)
{
    auto& data = _impl->metrics[static_cast<std::size_t>(metric)];
    data.count.fetch_add(1u, std::memory_order_relaxed);
    data.bytes.fetch_add(bytes, std::memory_order_relaxed);
    data.totalNanos.fetch_add(nanos, std::memory_order_relaxed);
    updateMaximum(data.maxNanos, nanos);
    if (outcome == Outcome::Success) data.success.fetch_add(1u, std::memory_order_relaxed);
    else if (outcome == Outcome::Miss) data.miss.fetch_add(1u, std::memory_order_relaxed);
    else data.error.fetch_add(1u, std::memory_order_relaxed);

    std::size_t bucket = 0u;
    for (std::uint64_t value = nanos; value > 1u && bucket + 1u < HISTOGRAM_BUCKETS; value >>= 1u)
        ++bucket;
    data.histogram[bucket].fetch_add(1u, std::memory_order_relaxed);

    const std::uint64_t event = _impl->events.fetch_add(1u, std::memory_order_relaxed) + 1u;
    if (_impl->intervalSeconds > 0u && (event & 1023u) == 0u)
    {
        const auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::steady_clock::now() - _impl->started).count();
        auto previous = _impl->lastReportSeconds.load(std::memory_order_relaxed);
        if (elapsed - previous >= static_cast<std::int64_t>(_impl->intervalSeconds) &&
            _impl->lastReportSeconds.compare_exchange_strong(
                previous, elapsed, std::memory_order_relaxed))
        {
            report(false);
        }
    }
}

void CacheStatistics::addBusyRetry(std::uint64_t count)
{
    _impl->busyRetries.fetch_add(count, std::memory_order_relaxed);
}

void CacheStatistics::addTransactionMutations(std::uint64_t count)
{
    _impl->transactionMutations.fetch_add(count, std::memory_order_relaxed);
}

void CacheStatistics::updateQueueHighWater(std::uint64_t bytes, std::uint64_t items)
{
    updateMaximum(_impl->queueHighWaterBytes, bytes);
    updateMaximum(_impl->queueHighWaterItems, items);
}

std::string CacheStatistics::toJSON(bool finalReport) const
{
    const auto elapsedNanos = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now() - _impl->started).count();

    std::ostringstream out;
    out << std::fixed << std::setprecision(3)
        << "{\"schema\":1"
        << ",\"timestamp\":\"" << jsonEscape(DateTime().asISO8601()) << "\""
        << ",\"final\":" << (finalReport ? "true" : "false")
        << ",\"driver\":\"" << jsonEscape(_impl->driver) << "\""
        << ",\"cache_path\":\"" << jsonEscape(_impl->cachePath) << "\""
        << ",\"pid\":" << processID()
        << ",\"elapsed_seconds\":" << (static_cast<double>(elapsedNanos) / 1.0e9)
        << ",\"busy_retries\":" << _impl->busyRetries.load(std::memory_order_relaxed)
        << ",\"transaction_mutations\":" << _impl->transactionMutations.load(std::memory_order_relaxed)
        << ",\"queue_high_water_bytes\":" << _impl->queueHighWaterBytes.load(std::memory_order_relaxed)
        << ",\"queue_high_water_items\":" << _impl->queueHighWaterItems.load(std::memory_order_relaxed)
        << ",\"metrics\":{";

    for (std::size_t m = 0u; m < static_cast<std::size_t>(Metric::Count); ++m)
    {
        if (m > 0u) out << ',';
        const auto& data = _impl->metrics[m];
        const std::uint64_t count = data.count.load(std::memory_order_relaxed);
        const std::uint64_t total = data.totalNanos.load(std::memory_order_relaxed);

        auto percentileUS = [&](double fraction)
        {
            if (count == 0u) return 0.0;
            const std::uint64_t target = static_cast<std::uint64_t>(std::ceil(count * fraction));
            std::uint64_t cumulative = 0u;
            for (std::size_t bucket = 0u; bucket < HISTOGRAM_BUCKETS; ++bucket)
            {
                cumulative += data.histogram[bucket].load(std::memory_order_relaxed);
                if (cumulative >= target)
                {
                    const double upperNanos = std::ldexp(1.0, static_cast<int>(bucket + 1u));
                    return upperNanos / 1000.0;
                }
            }
            return static_cast<double>(data.maxNanos.load(std::memory_order_relaxed)) / 1000.0;
        };

        out << '\"' << metricName(static_cast<Metric>(m)) << "\":{" 
            << "\"count\":" << count
            << ",\"success\":" << data.success.load(std::memory_order_relaxed)
            << ",\"miss\":" << data.miss.load(std::memory_order_relaxed)
            << ",\"error\":" << data.error.load(std::memory_order_relaxed)
            << ",\"bytes\":" << data.bytes.load(std::memory_order_relaxed)
            << ",\"mean_us\":" << (count ? static_cast<double>(total) / count / 1000.0 : 0.0)
            << ",\"max_us\":" << static_cast<double>(data.maxNanos.load(std::memory_order_relaxed)) / 1000.0
            << ",\"p50_us\":" << percentileUS(0.50)
            << ",\"p95_us\":" << percentileUS(0.95)
            << ",\"p99_us\":" << percentileUS(0.99)
            << '}';
    }
    out << "}}";
    return out.str();
}

void CacheStatistics::report(bool finalReport)
{
    if (finalReport && _impl->finalReported.exchange(true, std::memory_order_relaxed))
        return;

    std::lock_guard<std::mutex> reportLock(_impl->reportMutex);
    const std::string json = toJSON(finalReport);
    if (_impl->outputPath.empty())
    {
        OE_NOTICE << "[CacheStats] " << json << std::endl;
        return;
    }

    std::string path = replaceAll(_impl->outputPath, "{driver}", _impl->driver);
    path = replaceAll(path, "{pid}", std::to_string(processID()));
    osgEarth::makeDirectoryForFile(path);

    static std::mutex s_fileMutex;
    std::lock_guard<std::mutex> fileLock(s_fileMutex);
    std::ofstream output(path, std::ios::out | std::ios::app);
    if (output.is_open())
        output << json << '\n';
    else
        OE_WARN << "[CacheStats] Failed to append report to \"" << path << "\"" << std::endl;
}

CacheStatsScope::CacheStatsScope(
    const std::shared_ptr<CacheStatistics>& stats,
    CacheStatistics::Metric metric) :
    _stats(stats),
    _metric(metric)
{
    if (_stats)
        _start = std::chrono::steady_clock::now();
}

CacheStatsScope::~CacheStatsScope()
{
    finish();
}

void CacheStatsScope::finish()
{
    if (_stats)
    {
        const auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - _start).count();
        _stats->record(_metric, _outcome, static_cast<std::uint64_t>(nanos), _bytes);
        _stats.reset();
    }
}
