/* osgEarth
 * Copyright 2025 Pelican Mapping
 * MIT License
 */
#include <benchmark/benchmark.h>
#include <osgEarth/Cache>
#include <osgEarth/FileUtils>
#include <osgEarth/IOTypes>
#include <algorithm>
#include <atomic>
#include <string>
#include <thread>
#include <vector>

using namespace osgEarth;

namespace
{
    namespace fs
    {
        bool remove_all(const std::string& path)
        {
            return Util::removeDirectory(path);
        }
    }

    enum class Backend { FileSystem, SQLite3 };

    const char* driverName(Backend backend)
    {
        return backend == Backend::FileSystem ? "filesystem" : "sqlite3";
    }

    std::string uniquePath(Backend backend, const char* workload)
    {
        static std::atomic<std::uint64_t> s_sequence{ 0u };
        return std::string("cache_benchmark_") + driverName(backend) + '_' + workload + '_' +
            std::to_string(s_sequence.fetch_add(1u, std::memory_order_relaxed));
    }

    CacheOptions makeOptions(
        Backend backend,
        const std::string& path,
        unsigned backgroundThreads,
        bool separateBins = true)
    {
        Config config;
        config.set("path", path);
        config.set("threads", backgroundThreads);
        config.set("reader_connections", 8u);
        config.set("separate_bins", separateBins);
        config.set("stats", false);
        CacheOptions options(config);
        options.setDriver(driverName(backend));
        return options;
    }

    struct CacheHandle
    {
        osg::ref_ptr<Cache> cache;
        osg::ref_ptr<CacheBin> bin;

        bool open(const CacheOptions& options)
        {
            cache = CacheFactory::create(options);
            if (!cache.valid() || cache->getStatus().isError())
                return false;
            bin = cache->getOrCreateDefaultBin();
            return bin.valid();
        }

        void close()
        {
            bin = nullptr;
            cache = nullptr;
        }
    };

    struct MultiBinHandle
    {
        osg::ref_ptr<Cache> cache;
        std::vector<osg::ref_ptr<CacheBin>> bins;

        bool open(const CacheOptions& options, std::size_t binCount)
        {
            cache = CacheFactory::create(options);
            if (!cache.valid() || cache->getStatus().isError())
                return false;

            bins.reserve(binCount);
            for (std::size_t i = 0u; i < binCount; ++i)
            {
                osg::ref_ptr<CacheBin> bin = cache->addBin("layer_" + std::to_string(i));
                if (!bin.valid())
                {
                    close();
                    return false;
                }
                bins.push_back(std::move(bin));
            }
            return true;
        }

        void close()
        {
            bins.clear();
            cache = nullptr;
        }
    };

    std::string makePayload(std::size_t size)
    {
        static constexpr char alphabet[] =
            "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz-_";
        std::uint32_t random = 0x9e3779b9u;
        std::string result(size, 'x');
        for (std::size_t i = 0u; i < result.size(); ++i)
        {
            random = random * 1664525u + 1013904223u;
            result[i] = alphabet[(random >> 24u) & 63u];
        }
        return result;
    }

    std::vector<std::string> makeKeys(std::size_t count, const char* prefix = "key_")
    {
        std::vector<std::string> result;
        result.reserve(count);
        for (std::size_t i = 0u; i < count; ++i)
            result.emplace_back(std::string(prefix) + std::to_string(i));
        return result;
    }

    template<typename Function>
    void parallelFor(std::size_t count, unsigned threadCount, Function&& function)
    {
        threadCount = std::max(1u, threadCount);
        if (threadCount == 1u)
        {
            for (std::size_t i = 0u; i < count; ++i)
                function(i);
            return;
        }

        std::atomic<std::size_t> next{ 0u };
        std::atomic_bool start{ false };
        std::vector<std::thread> threads;
        threads.reserve(threadCount);
        for (unsigned t = 0u; t < threadCount; ++t)
        {
            threads.emplace_back([&]()
            {
                while (!start.load(std::memory_order_acquire))
                    std::this_thread::yield();
                for (;;)
                {
                    const std::size_t i = next.fetch_add(1u, std::memory_order_relaxed);
                    if (i >= count)
                        break;
                    function(i);
                }
            });
        }
        start.store(true, std::memory_order_release);
        for (auto& thread : threads)
            thread.join();
    }

    bool verifyStrings(
        Backend backend,
        const std::string& path,
        const std::vector<std::string>& keys,
        const std::string& expected)
    {
        CacheHandle handle;
        if (!handle.open(makeOptions(backend, path, 0u)))
            return false;

        bool valid = true;
        for (const auto& key : keys)
        {
            ReadResult result = handle.bin->readString(key, nullptr);
            const StringObject* value = result.get<StringObject>();
            if (!result.succeeded() || !value || value->getString() != expected)
            {
                valid = false;
                break;
            }
        }
        handle.close();
        return valid;
    }

    bool verifyMultiBinStrings(
        const std::string& path,
        bool separateBins,
        std::size_t binCount,
        const std::vector<std::string>& keys,
        const std::string& expected)
    {
        MultiBinHandle handle;
        if (!handle.open(makeOptions(Backend::SQLite3, path, 0u, separateBins), binCount))
            return false;

        bool valid = true;
        for (const auto& bin : handle.bins)
        {
            for (const auto& key : keys)
            {
                ReadResult result = bin->readString(key, nullptr);
                const StringObject* value = result.get<StringObject>();
                if (!result.succeeded() || !value || value->getString() != expected)
                {
                    valid = false;
                    break;
                }
            }
            if (!valid)
                break;
        }
        handle.close();
        return valid;
    }

    bool populateMultiBin(
        MultiBinHandle& handle,
        const std::vector<std::string>& keys,
        const osg::Object* object)
    {
        for (const auto& bin : handle.bins)
        {
            for (const auto& key : keys)
            {
                if (!bin->write(key, object, nullptr))
                    return false;
            }
        }
        return true;
    }

    template<Backend backend, bool durable>
    void cacheWrite(benchmark::State& state)
    {
        const std::size_t payloadBytes = static_cast<std::size_t>(state.range(0));
        const std::size_t operations = static_cast<std::size_t>(state.range(1));
        const unsigned callers = static_cast<unsigned>(state.range(2));
        const unsigned backgroundThreads = static_cast<unsigned>(state.range(3));

        std::uint64_t completed = 0u;
        for (auto _ : state)
        {
            state.PauseTiming();
            const std::string path = uniquePath(backend, durable ? "durable_write" : "submit_write");
            fs::remove_all(path);
            const std::string payload = makePayload(payloadBytes);
            osg::ref_ptr<StringObject> object = new StringObject(payload);
            const auto keys = makeKeys(operations);
            const CacheOptions options = makeOptions(backend, path, backgroundThreads);
            CacheHandle handle;

            if (!durable && !handle.open(options))
            {
                state.SkipWithError("Failed to open cache");
                fs::remove_all(path);
                return;
            }

            std::atomic<std::uint64_t> failures{ 0u };
            state.ResumeTiming();

            if (durable && !handle.open(options))
            {
                state.PauseTiming();
                state.SkipWithError("Failed to open cache");
                fs::remove_all(path);
                return;
            }

            parallelFor(operations, callers, [&](std::size_t i)
            {
                if (!handle.bin->write(keys[i], object.get(), nullptr))
                    failures.fetch_add(1u, std::memory_order_relaxed);
            });

            if (durable)
                handle.close(); // include queue drain and durable persistence

            state.PauseTiming();
            if (!durable)
                handle.close(); // submission latency intentionally excludes the drain

            if (failures.load(std::memory_order_relaxed) != 0u ||
                !verifyStrings(backend, path, keys, payload))
            {
                state.SkipWithError("Cache write validation failed");
                fs::remove_all(path);
                return;
            }

            fs::remove_all(path);
            completed += operations;
            state.ResumeTiming();
        }

        state.SetItemsProcessed(static_cast<std::int64_t>(completed));
        state.SetBytesProcessed(static_cast<std::int64_t>(completed * payloadBytes));
        state.counters["payload_bytes"] = static_cast<double>(payloadBytes);
        state.counters["caller_threads"] = static_cast<double>(callers);
        state.counters["background_threads"] = static_cast<double>(backgroundThreads);
    }

    template<Backend backend>
    void cacheWarmRead(benchmark::State& state)
    {
        const std::size_t payloadBytes = static_cast<std::size_t>(state.range(0));
        const std::size_t records = static_cast<std::size_t>(state.range(1));
        const std::size_t operations = static_cast<std::size_t>(state.range(2));
        const unsigned callers = static_cast<unsigned>(state.range(3));
        std::uint64_t completed = 0u;

        for (auto _ : state)
        {
            state.PauseTiming();
            const std::string path = uniquePath(backend, "warm_read");
            fs::remove_all(path);
            const std::string payload = makePayload(payloadBytes);
            osg::ref_ptr<StringObject> object = new StringObject(payload);
            const auto keys = makeKeys(records);
            const CacheOptions options = makeOptions(backend, path, 0u);

            CacheHandle populate;
            if (!populate.open(options))
            {
                state.SkipWithError("Failed to open cache for population");
                fs::remove_all(path);
                return;
            }
            bool populated = true;
            for (const auto& key : keys)
                populated = populate.bin->write(key, object.get(), nullptr) && populated;
            populate.close();

            CacheHandle handle;
            if (!populated || !handle.open(options))
            {
                state.SkipWithError("Failed to populate/reopen cache");
                fs::remove_all(path);
                return;
            }

            std::atomic<std::uint64_t> failures{ 0u };
            state.ResumeTiming();
            parallelFor(operations, callers, [&](std::size_t i)
            {
                ReadResult result = handle.bin->readString(keys[i % records], nullptr);
                const StringObject* value = result.get<StringObject>();
                if (!result.succeeded() || !value || value->getString().size() != payloadBytes)
                    failures.fetch_add(1u, std::memory_order_relaxed);
                benchmark::DoNotOptimize(value);
            });
            state.PauseTiming();

            handle.close();
            fs::remove_all(path);
            if (failures.load(std::memory_order_relaxed) != 0u)
            {
                state.SkipWithError("Warm cache read validation failed");
                return;
            }
            completed += operations;
            state.ResumeTiming();
        }

        state.SetItemsProcessed(static_cast<std::int64_t>(completed));
        state.SetBytesProcessed(static_cast<std::int64_t>(completed * payloadBytes));
        state.counters["payload_bytes"] = static_cast<double>(payloadBytes);
        state.counters["caller_threads"] = static_cast<double>(callers);
    }

    template<Backend backend>
    void cacheMiss(benchmark::State& state)
    {
        const std::size_t operations = static_cast<std::size_t>(state.range(0));
        const unsigned callers = static_cast<unsigned>(state.range(1));
        std::uint64_t completed = 0u;

        for (auto _ : state)
        {
            state.PauseTiming();
            const std::string path = uniquePath(backend, "miss");
            fs::remove_all(path);
            const auto keys = makeKeys(operations, "missing_");
            CacheHandle handle;
            if (!handle.open(makeOptions(backend, path, 0u)))
            {
                state.SkipWithError("Failed to open cache");
                fs::remove_all(path);
                return;
            }
            std::atomic<std::uint64_t> unexpectedHits{ 0u };
            state.ResumeTiming();
            parallelFor(operations, callers, [&](std::size_t i)
            {
                ReadResult result = handle.bin->readString(keys[i], nullptr);
                if (result.succeeded())
                    unexpectedHits.fetch_add(1u, std::memory_order_relaxed);
            });
            state.PauseTiming();
            handle.close();
            fs::remove_all(path);
            if (unexpectedHits.load(std::memory_order_relaxed) != 0u)
            {
                state.SkipWithError("Cache miss workload returned a hit");
                return;
            }
            completed += operations;
            state.ResumeTiming();
        }

        state.SetItemsProcessed(static_cast<std::int64_t>(completed));
        state.counters["caller_threads"] = static_cast<double>(callers);
    }

    template<Backend backend>
    void cacheMixed(benchmark::State& state)
    {
        const std::size_t payloadBytes = static_cast<std::size_t>(state.range(0));
        const std::size_t records = static_cast<std::size_t>(state.range(1));
        const std::size_t operations = static_cast<std::size_t>(state.range(2));
        const unsigned callers = static_cast<unsigned>(state.range(3));
        const unsigned backgroundThreads = static_cast<unsigned>(state.range(4));
        std::uint64_t completed = 0u;

        for (auto _ : state)
        {
            state.PauseTiming();
            const std::string path = uniquePath(backend, "mixed");
            fs::remove_all(path);
            const std::string payload = makePayload(payloadBytes);
            osg::ref_ptr<StringObject> object = new StringObject(payload);
            const auto keys = makeKeys(records);
            CacheHandle populate;
            if (!populate.open(makeOptions(backend, path, 0u)))
            {
                state.SkipWithError("Failed to open cache for population");
                fs::remove_all(path);
                return;
            }
            for (const auto& key : keys)
            {
                if (!populate.bin->write(key, object.get(), nullptr))
                {
                    state.SkipWithError("Failed to populate cache");
                    populate.close();
                    fs::remove_all(path);
                    return;
                }
            }
            populate.close();

            CacheHandle handle;
            if (!handle.open(makeOptions(backend, path, backgroundThreads)))
            {
                state.SkipWithError("Failed to reopen cache");
                fs::remove_all(path);
                return;
            }
            std::atomic<std::uint64_t> failures{ 0u };
            state.ResumeTiming();
            parallelFor(operations, callers, [&](std::size_t i)
            {
                const std::string& key = keys[i % records];
                if ((i % 20u) == 0u)
                {
                    if (!handle.bin->write(key, object.get(), nullptr))
                        failures.fetch_add(1u, std::memory_order_relaxed);
                }
                else
                {
                    ReadResult result = handle.bin->readString(key, nullptr);
                    if (!result.succeeded())
                        failures.fetch_add(1u, std::memory_order_relaxed);
                }
            });
            handle.close(); // total throughput includes durable drain of the 5% writes
            state.PauseTiming();

            if (failures.load(std::memory_order_relaxed) != 0u ||
                !verifyStrings(backend, path, keys, payload))
            {
                state.SkipWithError("Mixed workload validation failed");
                fs::remove_all(path);
                return;
            }
            fs::remove_all(path);
            completed += operations;
            state.ResumeTiming();
        }

        state.SetItemsProcessed(static_cast<std::int64_t>(completed));
        state.SetBytesProcessed(static_cast<std::int64_t>(completed * payloadBytes));
        state.counters["payload_bytes"] = static_cast<double>(payloadBytes);
        state.counters["read_percent"] = 95.0;
        state.counters["caller_threads"] = static_cast<double>(callers);
        state.counters["background_threads"] = static_cast<double>(backgroundThreads);
    }

    template<bool separateBins>
    void sqliteMultiBinDurableWrite(benchmark::State& state)
    {
        const std::size_t payloadBytes = static_cast<std::size_t>(state.range(0));
        const std::size_t binCount = static_cast<std::size_t>(state.range(1));
        const std::size_t recordsPerBin = static_cast<std::size_t>(state.range(2));
        const unsigned callers = static_cast<unsigned>(state.range(3));
        const unsigned backgroundThreads = static_cast<unsigned>(state.range(4));
        const std::size_t operations = binCount * recordsPerBin;
        std::uint64_t completed = 0u;

        for (auto _ : state)
        {
            state.PauseTiming();
            const std::string path = uniquePath(
                Backend::SQLite3, separateBins ? "multibin_write_separate" : "multibin_write_shared");
            fs::remove_all(path);
            const std::string payload = makePayload(payloadBytes);
            osg::ref_ptr<StringObject> object = new StringObject(payload);
            const auto keys = makeKeys(recordsPerBin);
            const CacheOptions options = makeOptions(
                Backend::SQLite3, path, backgroundThreads, separateBins);
            MultiBinHandle handle;
            std::atomic<std::uint64_t> failures{ 0u };

            state.ResumeTiming();
            if (!handle.open(options, binCount))
            {
                state.PauseTiming();
                state.SkipWithError("Failed to open multi-bin cache");
                fs::remove_all(path);
                return;
            }
            parallelFor(operations, callers, [&](std::size_t i)
            {
                const std::size_t binIndex = i % binCount;
                const std::size_t keyIndex = i / binCount;
                if (!handle.bins[binIndex]->write(keys[keyIndex], object.get(), nullptr))
                    failures.fetch_add(1u, std::memory_order_relaxed);
            });
            handle.close(); // includes all database drains and closes
            state.PauseTiming();

            if (failures.load(std::memory_order_relaxed) != 0u ||
                !verifyMultiBinStrings(path, separateBins, binCount, keys, payload))
            {
                state.SkipWithError("Multi-bin durable write validation failed");
                fs::remove_all(path);
                return;
            }
            fs::remove_all(path);
            completed += operations;
            state.ResumeTiming();
        }

        state.SetItemsProcessed(static_cast<std::int64_t>(completed));
        state.SetBytesProcessed(static_cast<std::int64_t>(completed * payloadBytes));
        state.counters["bins"] = static_cast<double>(binCount);
        state.counters["database_files"] = static_cast<double>(separateBins ? binCount : 1u);
        state.counters["caller_threads"] = static_cast<double>(callers);
        state.counters["background_threads"] = static_cast<double>(backgroundThreads);
    }

    template<bool separateBins>
    void sqliteMultiBinWarmRead(benchmark::State& state)
    {
        const std::size_t payloadBytes = static_cast<std::size_t>(state.range(0));
        const std::size_t binCount = static_cast<std::size_t>(state.range(1));
        const std::size_t recordsPerBin = static_cast<std::size_t>(state.range(2));
        const std::size_t operations = static_cast<std::size_t>(state.range(3));
        const unsigned callers = static_cast<unsigned>(state.range(4));
        std::uint64_t completed = 0u;

        for (auto _ : state)
        {
            state.PauseTiming();
            const std::string path = uniquePath(
                Backend::SQLite3, separateBins ? "multibin_read_separate" : "multibin_read_shared");
            fs::remove_all(path);
            const std::string payload = makePayload(payloadBytes);
            osg::ref_ptr<StringObject> object = new StringObject(payload);
            const auto keys = makeKeys(recordsPerBin);
            const CacheOptions options = makeOptions(Backend::SQLite3, path, 0u, separateBins);

            MultiBinHandle populate;
            if (!populate.open(options, binCount) ||
                !populateMultiBin(populate, keys, object.get()))
            {
                state.SkipWithError("Failed to populate multi-bin cache");
                populate.close();
                fs::remove_all(path);
                return;
            }
            populate.close();

            MultiBinHandle handle;
            if (!handle.open(options, binCount))
            {
                state.SkipWithError("Failed to reopen multi-bin cache");
                fs::remove_all(path);
                return;
            }
            std::atomic<std::uint64_t> failures{ 0u };
            state.ResumeTiming();
            parallelFor(operations, callers, [&](std::size_t i)
            {
                const std::size_t binIndex = i % binCount;
                const std::size_t keyIndex = (i / binCount) % recordsPerBin;
                ReadResult result = handle.bins[binIndex]->readString(keys[keyIndex], nullptr);
                const StringObject* value = result.get<StringObject>();
                if (!result.succeeded() || !value || value->getString().size() != payloadBytes)
                    failures.fetch_add(1u, std::memory_order_relaxed);
                benchmark::DoNotOptimize(value);
            });
            state.PauseTiming();

            handle.close();
            fs::remove_all(path);
            if (failures.load(std::memory_order_relaxed) != 0u)
            {
                state.SkipWithError("Multi-bin warm read validation failed");
                return;
            }
            completed += operations;
            state.ResumeTiming();
        }

        state.SetItemsProcessed(static_cast<std::int64_t>(completed));
        state.SetBytesProcessed(static_cast<std::int64_t>(completed * payloadBytes));
        state.counters["bins"] = static_cast<double>(binCount);
        state.counters["database_files"] = static_cast<double>(separateBins ? binCount : 1u);
        state.counters["caller_threads"] = static_cast<double>(callers);
    }

    template<bool separateBins>
    void sqliteMultiBinMixed(benchmark::State& state)
    {
        const std::size_t payloadBytes = static_cast<std::size_t>(state.range(0));
        const std::size_t binCount = static_cast<std::size_t>(state.range(1));
        const std::size_t recordsPerBin = static_cast<std::size_t>(state.range(2));
        const std::size_t operations = static_cast<std::size_t>(state.range(3));
        const unsigned callers = static_cast<unsigned>(state.range(4));
        const unsigned backgroundThreads = static_cast<unsigned>(state.range(5));
        std::uint64_t completed = 0u;

        for (auto _ : state)
        {
            state.PauseTiming();
            const std::string path = uniquePath(
                Backend::SQLite3, separateBins ? "multibin_mixed_separate" : "multibin_mixed_shared");
            fs::remove_all(path);
            const std::string payload = makePayload(payloadBytes);
            osg::ref_ptr<StringObject> object = new StringObject(payload);
            const auto keys = makeKeys(recordsPerBin);

            MultiBinHandle populate;
            if (!populate.open(
                    makeOptions(Backend::SQLite3, path, 0u, separateBins), binCount) ||
                !populateMultiBin(populate, keys, object.get()))
            {
                state.SkipWithError("Failed to populate multi-bin cache");
                populate.close();
                fs::remove_all(path);
                return;
            }
            populate.close();

            MultiBinHandle handle;
            if (!handle.open(
                    makeOptions(Backend::SQLite3, path, backgroundThreads, separateBins), binCount))
            {
                state.SkipWithError("Failed to reopen multi-bin cache");
                fs::remove_all(path);
                return;
            }
            std::atomic<std::uint64_t> failures{ 0u };
            state.ResumeTiming();
            parallelFor(operations, callers, [&](std::size_t i)
            {
                // Keep every twentieth operation a write while rotating those writes
                // across every bin, including power-of-two bin counts.
                const std::size_t binIndex = (i * 7u + i / 20u) % binCount;
                const std::size_t keyIndex = (i / binCount) % recordsPerBin;
                CacheBin* bin = handle.bins[binIndex].get();
                if ((i % 20u) == 0u)
                {
                    if (!bin->write(keys[keyIndex], object.get(), nullptr))
                        failures.fetch_add(1u, std::memory_order_relaxed);
                }
                else
                {
                    ReadResult result = bin->readString(keys[keyIndex], nullptr);
                    if (!result.succeeded())
                        failures.fetch_add(1u, std::memory_order_relaxed);
                }
            });
            handle.close(); // includes durable drain of the 5% writes
            state.PauseTiming();

            if (failures.load(std::memory_order_relaxed) != 0u ||
                !verifyMultiBinStrings(path, separateBins, binCount, keys, payload))
            {
                state.SkipWithError("Multi-bin mixed workload validation failed");
                fs::remove_all(path);
                return;
            }
            fs::remove_all(path);
            completed += operations;
            state.ResumeTiming();
        }

        state.SetItemsProcessed(static_cast<std::int64_t>(completed));
        state.SetBytesProcessed(static_cast<std::int64_t>(completed * payloadBytes));
        state.counters["payload_bytes"] = static_cast<double>(payloadBytes);
        state.counters["bins"] = static_cast<double>(binCount);
        state.counters["database_files"] = static_cast<double>(separateBins ? binCount : 1u);
        state.counters["read_percent"] = 95.0;
        state.counters["caller_threads"] = static_cast<double>(callers);
        state.counters["background_threads"] = static_cast<double>(backgroundThreads);
    }

    void writeArguments(benchmark::internal::Benchmark* benchmark)
    {
        benchmark
            ->Args({ 65536, 256, 4, 4 })
            ->Args({ 131072, 128, 4, 4 })
            ->Args({ 262144, 64, 4, 4 })
            ->ArgNames({ "bytes", "operations", "callers", "writers" })
            ->Iterations(1)->Repetitions(3)->UseRealTime();
    }

    void submitArguments(benchmark::internal::Benchmark* benchmark)
    {
        benchmark
            ->Args({ 65536, 256, 4, 4 })
            ->Args({ 131072, 128, 4, 4 })
            ->Args({ 262144, 64, 4, 4 })
            ->ArgNames({ "bytes", "operations", "callers", "writers" })
            ->Iterations(1)->Repetitions(3)->UseRealTime();
    }

    void readArguments(benchmark::internal::Benchmark* benchmark)
    {
        benchmark
            ->Args({ 65536, 256, 2048, 8 })
            ->Args({ 131072, 128, 1024, 8 })
            ->Args({ 262144, 64, 512, 8 })
            ->ArgNames({ "bytes", "records", "operations", "callers" })
            ->Iterations(1)->Repetitions(3)->UseRealTime();
    }

    void missArguments(benchmark::internal::Benchmark* benchmark)
    {
        benchmark
            ->Args({ 5000, 1 })
            ->Args({ 5000, 8 })
            ->ArgNames({ "operations", "callers" })
            ->Iterations(1)->Repetitions(3)->UseRealTime();
    }

    void mixedArguments(benchmark::internal::Benchmark* benchmark)
    {
        benchmark
            ->Args({ 65536, 256, 2048, 8, 4 })
            ->Args({ 131072, 128, 1024, 8, 4 })
            ->Args({ 262144, 64, 512, 8, 4 })
            ->ArgNames({ "bytes", "records", "operations", "callers", "writers" })
            ->Iterations(1)->Repetitions(3)->UseRealTime();
    }

    void multiBinWriteArguments(benchmark::internal::Benchmark* benchmark)
    {
        benchmark
            ->Args({ 65536, 2, 128, 8, 4 })
            ->Args({ 65536, 8, 32, 8, 4 })
            ->Args({ 65536, 16, 16, 8, 4 })
            ->Args({ 131072, 2, 64, 8, 4 })
            ->Args({ 131072, 8, 16, 8, 4 })
            ->Args({ 131072, 16, 8, 8, 4 })
            ->Args({ 262144, 2, 32, 8, 4 })
            ->Args({ 262144, 8, 8, 8, 4 })
            ->Args({ 262144, 16, 4, 8, 4 })
            ->ArgNames({ "bytes", "bins", "records_per_bin", "callers", "writers" })
            ->Iterations(1)->Repetitions(3)->UseRealTime();
    }

    void multiBinReadArguments(benchmark::internal::Benchmark* benchmark)
    {
        benchmark
            ->Args({ 65536, 2, 128, 2048, 8 })
            ->Args({ 65536, 8, 32, 2048, 8 })
            ->Args({ 65536, 16, 16, 2048, 8 })
            ->Args({ 131072, 2, 64, 1024, 8 })
            ->Args({ 131072, 8, 16, 1024, 8 })
            ->Args({ 131072, 16, 8, 1024, 8 })
            ->Args({ 262144, 2, 32, 512, 8 })
            ->Args({ 262144, 8, 8, 512, 8 })
            ->Args({ 262144, 16, 4, 512, 8 })
            ->ArgNames({ "bytes", "bins", "records_per_bin", "operations", "callers" })
            ->Iterations(1)->Repetitions(3)->UseRealTime();
    }

    void multiBinMixedArguments(benchmark::internal::Benchmark* benchmark)
    {
        benchmark
            ->Args({ 65536, 2, 128, 2048, 8, 4 })
            ->Args({ 65536, 8, 32, 2048, 8, 4 })
            ->Args({ 65536, 16, 16, 2048, 8, 4 })
            ->Args({ 131072, 2, 64, 1024, 8, 4 })
            ->Args({ 131072, 8, 16, 1024, 8, 4 })
            ->Args({ 131072, 16, 8, 1024, 8, 4 })
            ->Args({ 262144, 2, 32, 512, 8, 4 })
            ->Args({ 262144, 8, 8, 512, 8, 4 })
            ->Args({ 262144, 16, 4, 512, 8, 4 })
            ->ArgNames({ "bytes", "bins", "records_per_bin", "operations", "callers", "writers" })
            ->Iterations(1)->Repetitions(3)->UseRealTime();
    }
}

BENCHMARK_TEMPLATE(cacheWrite, Backend::FileSystem, true)
    ->Name("Cache/DurableWrite/filesystem")->Apply(writeArguments);
BENCHMARK_TEMPLATE(cacheWrite, Backend::SQLite3, true)
    ->Name("Cache/DurableWrite/sqlite3")->Apply(writeArguments);

BENCHMARK_TEMPLATE(cacheWrite, Backend::FileSystem, false)
    ->Name("Cache/SubmitWrite/filesystem")->Apply(submitArguments);
BENCHMARK_TEMPLATE(cacheWrite, Backend::SQLite3, false)
    ->Name("Cache/SubmitWrite/sqlite3")->Apply(submitArguments);

BENCHMARK_TEMPLATE(cacheWarmRead, Backend::FileSystem)
    ->Name("Cache/WarmRead/filesystem")->Apply(readArguments);
BENCHMARK_TEMPLATE(cacheWarmRead, Backend::SQLite3)
    ->Name("Cache/WarmRead/sqlite3")->Apply(readArguments);

BENCHMARK_TEMPLATE(cacheMiss, Backend::FileSystem)
    ->Name("Cache/Miss/filesystem")->Apply(missArguments);
BENCHMARK_TEMPLATE(cacheMiss, Backend::SQLite3)
    ->Name("Cache/Miss/sqlite3")->Apply(missArguments);

BENCHMARK_TEMPLATE(cacheMixed, Backend::FileSystem)
    ->Name("Cache/Mixed95Read5Write/filesystem")->Apply(mixedArguments);
BENCHMARK_TEMPLATE(cacheMixed, Backend::SQLite3)
    ->Name("Cache/Mixed95Read5Write/sqlite3")->Apply(mixedArguments);

BENCHMARK_TEMPLATE(sqliteMultiBinDurableWrite, true)
    ->Name("Cache/MultiBinDurableWrite/sqlite3/separate")->Apply(multiBinWriteArguments);
BENCHMARK_TEMPLATE(sqliteMultiBinDurableWrite, false)
    ->Name("Cache/MultiBinDurableWrite/sqlite3/shared")->Apply(multiBinWriteArguments);

BENCHMARK_TEMPLATE(sqliteMultiBinWarmRead, true)
    ->Name("Cache/MultiBinWarmRead/sqlite3/separate")->Apply(multiBinReadArguments);
BENCHMARK_TEMPLATE(sqliteMultiBinWarmRead, false)
    ->Name("Cache/MultiBinWarmRead/sqlite3/shared")->Apply(multiBinReadArguments);

BENCHMARK_TEMPLATE(sqliteMultiBinMixed, true)
    ->Name("Cache/MultiBinMixed95Read5Write/sqlite3/separate")->Apply(multiBinMixedArguments);
BENCHMARK_TEMPLATE(sqliteMultiBinMixed, false)
    ->Name("Cache/MultiBinMixed95Read5Write/sqlite3/shared")->Apply(multiBinMixedArguments);
