/* osgEarth
 * Copyright 2025 Pelican Mapping
 * MIT License
 */
#include "SQLite3Cache"
#include <osgEarth/Cache>
#include <osgEarth/StringUtils>
#include <osgEarth/Threading>
#include <osgEarth/URI>
#include <osgEarth/FileUtils>
#include <osgEarth/Registry>
#include <osgEarth/NetworkMonitor>
#include <osgEarth/Metrics>
#include <osgDB/FileUtils>
#include <osgDB/FileNameUtils>
#include <sqlite3.h>
#include <sstream>

using namespace osgEarth;
using namespace osgEarth::Drivers;

#undef  LC
#define LC "[SQLite3Cache] "

#define OSG_FORMAT "osgb"

namespace
{
    struct WriteCacheRecord {
        Config meta;
        osg::ref_ptr<const osg::Object> object;
    };
    typedef std::unordered_map<std::string, WriteCacheRecord> WriteCache;

    class SQLite3CacheBin; // forward

    // Per-thread set of prepared statements for a single cache bin.
    // Each thread lazily prepares its own statements so that no
    // application-level mutex is needed — SQLITE_OPEN_FULLMUTEX
    // handles connection-level thread safety internally.
    struct StmtSet {
        sqlite3_stmt* selectStmt = nullptr;
        sqlite3_stmt* insertStmt = nullptr;
        sqlite3_stmt* deleteStmt = nullptr;
        sqlite3_stmt* touchStmt = nullptr;
        sqlite3_stmt* existsStmt = nullptr;
        sqlite3_stmt* clearStmt = nullptr;
        sqlite3_stmt* sizeStmt = nullptr;
        sqlite3* preparedFor = nullptr;

        ~StmtSet() { finalize(); }

        void finalize() {
            if (selectStmt) { sqlite3_finalize(selectStmt); selectStmt = nullptr; }
            if (insertStmt) { sqlite3_finalize(insertStmt); insertStmt = nullptr; }
            if (deleteStmt) { sqlite3_finalize(deleteStmt); deleteStmt = nullptr; }
            if (touchStmt)  { sqlite3_finalize(touchStmt);  touchStmt = nullptr; }
            if (existsStmt) { sqlite3_finalize(existsStmt); existsStmt = nullptr; }
            if (clearStmt)  { sqlite3_finalize(clearStmt);  clearStmt = nullptr; }
            if (sizeStmt)   { sqlite3_finalize(sizeStmt);   sizeStmt = nullptr; }
            preparedFor = nullptr;
        }
    };

    // thread_local cache of prepared statements, keyed by bin pointer.
    static thread_local std::unordered_map<const SQLite3CacheBin*, StmtSet> t_stmtCache;

    // ---------------------------------------------------------------
    // SQLite3Cache
    // ---------------------------------------------------------------

    class SQLite3Cache : public Cache
    {
    public:
        SQLite3Cache() : _db(nullptr) { }
        SQLite3Cache( const SQLite3Cache& rhs, const osg::CopyOp& op ) : _db(nullptr) { }
        META_Object( osgEarth, SQLite3Cache );

        SQLite3Cache( const CacheOptions& options );
        ~SQLite3Cache() override;

    public: // Cache interface
        CacheBin* addBin( const std::string& binID ) override;
        CacheBin* getOrCreateDefaultBin() override;
        void setNumThreads(unsigned) override;
        bool compact() override;

    protected:
        friend class SQLite3CacheBin;
        std::string _rootPath;
        SQLite3CacheOptions _options;
        jobs::jobpool* _pool = nullptr;
        sqlite3* _db;             // shared db (null when separateBins is true)
        std::mutex _compactMutex; // only used for VACUUM on shared db
    };

    // ---------------------------------------------------------------
    // SQLite3CacheBin
    // ---------------------------------------------------------------

    class SQLite3CacheBin : public CacheBin
    {
    public:
        // Shared-database constructor (separateBins == false)
        SQLite3CacheBin(
            const std::string& binID,
            sqlite3* db,
            const SQLite3CacheOptions& options,
            jobs::jobpool* pool);

        // Per-bin database constructor (separateBins == true)
        SQLite3CacheBin(
            const std::string& binID,
            const std::string& rootPath,
            const SQLite3CacheOptions& options,
            jobs::jobpool* pool);

        ~SQLite3CacheBin() override;

    public: // CacheBin interface
        ReadResult readObject(const std::string& key, const osgDB::Options* dbo) override;
        ReadResult readImage(const std::string& key, const osgDB::Options* dbo) override;
        ReadResult readString(const std::string& key, const osgDB::Options* dbo) override;

        bool write(
            const std::string& key,
            const osg::Object* object,
            const Config& meta,
            const osgDB::Options* dbo) override;

        bool remove(const std::string& key) override;
        bool touch(const std::string& key) override;
        RecordStatus getRecordStatus(const std::string& key) override;
        bool clear() override;
        bool compact() override;
        unsigned getStorageSize() override;

    private:
        void initReaderWriter();
        void prepareStatements(StmtSet& stmts);
        ReadResult read(const std::string& key, const osgDB::Options* dbo, bool isImage);

        sqlite3*& db() { return _ownDb ? _ownDb : _db; }

        // Get (or lazily create) the calling thread's prepared statements.
        StmtSet& threadStatements() {
            StmtSet& s = t_stmtCache[this];
            if (s.preparedFor != db()) {
                s.finalize();
                prepareStatements(s);
                s.preparedFor = db();
            }
            return s;
        }

        sqlite3* _db;               // shared db (non-owning)
        sqlite3* _ownDb;            // per-bin db (owning, null in shared mode)

        bool _separate;             // true when using per-bin database
        jobs::jobpool* _pool;
        SQLite3CacheOptions _options;

        osg::ref_ptr<osgDB::ReaderWriter> _rw;
        osg::ref_ptr<osgDB::Options> _rwOptions;
        std::string _compressorName;

        WriteCache _writeCache;
        ReadWriteMutex _writeCacheRWM;
        std::mutex _compactMutex; // only for VACUUM

        bool _ok;
    };

    // ---------------------------------------------------------------
    // SQLite3Cache implementation
    // ---------------------------------------------------------------

    SQLite3Cache::SQLite3Cache( const CacheOptions& options ) :
        Cache(options),
        _options(options),
        _db(nullptr),
        _pool(nullptr)
    {
        if ( !_options.rootPath().isSet() )
        {
            const char* cachePath = ::getenv(OSGEARTH_ENV_CACHE_PATH);
            if ( cachePath )
                _options.rootPath() = cachePath;
        }

        if ( !_options.rootPath().isSet() )
        {
            _status.set(Status::ConfigurationError, "No cache path specified");
            return;
        }

        _rootPath = URI( *_options.rootPath(), options.referrer() ).full();

        if (osgDB::makeDirectory(_rootPath) == false)
        {
            _status.set(Status::ResourceUnavailable, Stringify()
                << "Failed to create or access folder \"" << _rootPath << "\"");
            return;
        }

        // In separate-bins mode each CacheBin opens its own database,
        // so we skip creating the shared database here.
        if (_options.separateBins() == true)
        {
            OE_INFO << LC << "SQLite3 cache (separate bins) at \"" << _rootPath << "\"" << std::endl;
            setNumThreads(_options.threads().get());
            return;
        }

        std::string dbPath = osgDB::concatPaths(_rootPath, "osgearth_cache.db");

        int flags = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_FULLMUTEX;
        int rc = sqlite3_open_v2(dbPath.c_str(), &_db, flags, nullptr);
        if (rc != SQLITE_OK)
        {
            _status.set(Status::ResourceUnavailable, Stringify()
                << "Failed to open SQLite3 database: " << sqlite3_errmsg(_db));
            sqlite3_close(_db);
            _db = nullptr;
            return;
        }

        // Configure for multi-process concurrent access.
        // Use the C API call for busy timeout — it installs a busy handler
        // that retries internally on every sqlite3_step, which is more
        // reliable than the PRAGMA when using prepared statements.
        sqlite3_busy_timeout(_db, (int)_options.busyTimeout().get());

        char* errMsg = nullptr;

        std::string pragmas =
            "PRAGMA journal_mode=WAL;"
            "PRAGMA synchronous=NORMAL;";

        rc = sqlite3_exec(_db, pragmas.c_str(), nullptr, nullptr, &errMsg);
        if (rc != SQLITE_OK)
        {
            OE_WARN << LC << "PRAGMA error: " << (errMsg ? errMsg : "unknown") << std::endl;
            sqlite3_free(errMsg);
        }

        // Create schema
        const char* schema =
            "CREATE TABLE IF NOT EXISTS cache ("
            "  bin_id    TEXT NOT NULL,"
            "  key       TEXT NOT NULL,"
            "  data      BLOB NOT NULL,"
            "  metadata  TEXT,"
            "  timestamp INTEGER NOT NULL DEFAULT (strftime('%s','now')),"
            "  PRIMARY KEY (bin_id, key)"
            ");"
            "CREATE INDEX IF NOT EXISTS idx_cache_timestamp ON cache(bin_id, timestamp);";

        rc = sqlite3_exec(_db, schema, nullptr, nullptr, &errMsg);
        if (rc != SQLITE_OK)
        {
            _status.set(Status::ResourceUnavailable, Stringify()
                << "Failed to create cache schema: " << (errMsg ? errMsg : "unknown"));
            sqlite3_free(errMsg);
            sqlite3_close(_db);
            _db = nullptr;
            return;
        }

        OE_INFO << LC << "Opened SQLite3 cache at \"" << dbPath << "\"" << std::endl;

        setNumThreads(_options.threads().get());
    }

    SQLite3Cache::~SQLite3Cache()
    {
        if (_db)
        {
            sqlite3_close_v2(_db);
            _db = nullptr;
        }
    }

    void
    SQLite3Cache::setNumThreads(unsigned num)
    {
        if (num > 0u)
        {
            _pool = jobs::get_pool("oe.sqlite3cache");
            _pool->set_can_steal_work(false);
            _pool->set_concurrency(osg::clampBetween(num, 1u, 8u));
        }
        else
        {
            _pool = nullptr;
        }
    }

    CacheBin*
    SQLite3Cache::addBin( const std::string& name )
    {
        if (getStatus().isError())
            return nullptr;

        if (_options.separateBins() == true)
        {
            return _bins.getOrCreate(name,
                new SQLite3CacheBin(name, _rootPath, _options, _pool));
        }
        else
        {
            if (!_db) return nullptr;
            return _bins.getOrCreate(name,
                new SQLite3CacheBin(name, _db, _options, _pool));
        }
    }

    CacheBin*
    SQLite3Cache::getOrCreateDefaultBin()
    {
        if (getStatus().isError())
            return nullptr;

        static Mutex s_defaultBinMutex;
        if ( !_defaultBin.valid() )
        {
            std::lock_guard<std::mutex> lock( s_defaultBinMutex );
            if ( !_defaultBin.valid() )
            {
                if (_options.separateBins() == true)
                    _defaultBin = new SQLite3CacheBin("__default", _rootPath, _options, _pool);
                else if (_db)
                    _defaultBin = new SQLite3CacheBin("__default", _db, _options, _pool);
            }
        }
        return _defaultBin.get();
    }

    bool
    SQLite3Cache::compact()
    {
        if (_options.separateBins() == true)
        {
            // Compact each bin's own database
            bool ok = true;
            _bins.forEach([&ok](const std::string&, osg::ref_ptr<CacheBin>& bin) {
                if (!bin->compact()) ok = false;
            });
            return ok;
        }

        if (!_db) return false;
        std::lock_guard<std::mutex> lock(_compactMutex);
        char* errMsg = nullptr;
        int rc = sqlite3_exec(_db, "PRAGMA optimize; VACUUM;", nullptr, nullptr, &errMsg);
        if (rc != SQLITE_OK)
        {
            OE_WARN << LC << "Compact failed: " << (errMsg ? errMsg : "unknown") << std::endl;
            sqlite3_free(errMsg);
            return false;
        }
        return true;
    }

    // ---------------------------------------------------------------
    // SQLite3CacheBin implementation
    // ---------------------------------------------------------------

    void
    SQLite3CacheBin::initReaderWriter()
    {
        _rw = osgDB::Registry::instance()->getReaderWriterForExtension(OSG_FORMAT);
        if (!_rw.valid())
        {
            OE_WARN << LC << "Failed to find ReaderWriter for \"" << OSG_FORMAT << "\"" << std::endl;
            _ok = false;
            return;
        }

        _rwOptions = Registry::instance()->cloneOrCreateOptions();

        if (::getenv(OSGEARTH_ENV_DEFAULT_COMPRESSOR) != nullptr)
        {
            _compressorName = ::getenv(OSGEARTH_ENV_DEFAULT_COMPRESSOR);
        }
        else
        {
            _compressorName = "zlib";
        }

        if (!_compressorName.empty())
        {
            _rwOptions->setPluginStringData("Compressor", _compressorName);
        }
    }

    void
    SQLite3CacheBin::prepareStatements(StmtSet& stmts)
    {
        if (_separate)
        {
            // Per-bin database: no bin_id column
            sqlite3_prepare_v2(db(),
                "SELECT data, metadata, timestamp FROM cache WHERE key=?",
                -1, &stmts.selectStmt, nullptr);

            sqlite3_prepare_v2(db(),
                "INSERT OR REPLACE INTO cache (key, data, metadata, timestamp) "
                "VALUES (?, ?, ?, strftime('%s','now'))",
                -1, &stmts.insertStmt, nullptr);

            sqlite3_prepare_v2(db(),
                "DELETE FROM cache WHERE key=?",
                -1, &stmts.deleteStmt, nullptr);

            sqlite3_prepare_v2(db(),
                "UPDATE cache SET timestamp=strftime('%s','now') WHERE key=?",
                -1, &stmts.touchStmt, nullptr);

            sqlite3_prepare_v2(db(),
                "SELECT 1 FROM cache WHERE key=?",
                -1, &stmts.existsStmt, nullptr);

            sqlite3_prepare_v2(db(),
                "DELETE FROM cache",
                -1, &stmts.clearStmt, nullptr);

            sqlite3_prepare_v2(db(),
                "SELECT SUM(LENGTH(data)) FROM cache",
                -1, &stmts.sizeStmt, nullptr);
        }
        else
        {
            // Shared database: filter by bin_id
            sqlite3_prepare_v2(db(),
                "SELECT data, metadata, timestamp FROM cache WHERE bin_id=? AND key=?",
                -1, &stmts.selectStmt, nullptr);

            sqlite3_prepare_v2(db(),
                "INSERT OR REPLACE INTO cache (bin_id, key, data, metadata, timestamp) "
                "VALUES (?, ?, ?, ?, strftime('%s','now'))",
                -1, &stmts.insertStmt, nullptr);

            sqlite3_prepare_v2(db(),
                "DELETE FROM cache WHERE bin_id=? AND key=?",
                -1, &stmts.deleteStmt, nullptr);

            sqlite3_prepare_v2(db(),
                "UPDATE cache SET timestamp=strftime('%s','now') WHERE bin_id=? AND key=?",
                -1, &stmts.touchStmt, nullptr);

            sqlite3_prepare_v2(db(),
                "SELECT 1 FROM cache WHERE bin_id=? AND key=?",
                -1, &stmts.existsStmt, nullptr);

            sqlite3_prepare_v2(db(),
                "DELETE FROM cache WHERE bin_id=?",
                -1, &stmts.clearStmt, nullptr);

            sqlite3_prepare_v2(db(),
                "SELECT SUM(LENGTH(data)) FROM cache WHERE bin_id=?",
                -1, &stmts.sizeStmt, nullptr);
        }
    }

    // Shared-database constructor
    SQLite3CacheBin::SQLite3CacheBin(
        const std::string& binID,
        sqlite3* db,
        const SQLite3CacheOptions& options,
        jobs::jobpool* pool) :

        CacheBin(binID, options.enableNodeCaching().get()),
        _db(db),
        _ownDb(nullptr),
        _separate(false),
        _pool(pool),
        _options(options),
        _ok(true)
    {
        initReaderWriter();
    }

    // Per-bin database constructor
    SQLite3CacheBin::SQLite3CacheBin(
        const std::string& binID,
        const std::string& rootPath,
        const SQLite3CacheOptions& options,
        jobs::jobpool* pool) :

        CacheBin(binID, options.enableNodeCaching().get()),
        _db(nullptr),
        _ownDb(nullptr),
        _separate(true),
        _pool(pool),
        _options(options),
        _ok(true)
    {
        initReaderWriter();
        if (!_ok) return;

        std::string dbPath = osgDB::concatPaths(rootPath, binID + ".db");

        int flags = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_FULLMUTEX;
        int rc = sqlite3_open_v2(dbPath.c_str(), &_ownDb, flags, nullptr);
        if (rc != SQLITE_OK)
        {
            OE_WARN << LC << "Failed to open per-bin database: " << sqlite3_errmsg(_ownDb) << std::endl;
            sqlite3_close(_ownDb);
            _ownDb = nullptr;
            _ok = false;
            return;
        }

        sqlite3_busy_timeout(_ownDb, (int)_options.busyTimeout().get());

        char* errMsg = nullptr;
        std::string pragmas =
            "PRAGMA journal_mode=WAL;"
            "PRAGMA synchronous=NORMAL;";

        rc = sqlite3_exec(_ownDb, pragmas.c_str(), nullptr, nullptr, &errMsg);
        if (rc != SQLITE_OK)
        {
            OE_WARN << LC << "PRAGMA error: " << (errMsg ? errMsg : "unknown") << std::endl;
            sqlite3_free(errMsg);
        }

        const char* schema =
            "CREATE TABLE IF NOT EXISTS cache ("
            "  key       TEXT NOT NULL PRIMARY KEY,"
            "  data      BLOB NOT NULL,"
            "  metadata  TEXT,"
            "  timestamp INTEGER NOT NULL DEFAULT (strftime('%s','now'))"
            ");"
            "CREATE INDEX IF NOT EXISTS idx_cache_timestamp ON cache(timestamp);";

        rc = sqlite3_exec(_ownDb, schema, nullptr, nullptr, &errMsg);
        if (rc != SQLITE_OK)
        {
            OE_WARN << LC << "Failed to create per-bin schema: " << (errMsg ? errMsg : "unknown") << std::endl;
            sqlite3_free(errMsg);
            sqlite3_close(_ownDb);
            _ownDb = nullptr;
            _ok = false;
            return;
        }

        OE_INFO << LC << "Opened per-bin database at \"" << dbPath << "\"" << std::endl;
    }

    SQLite3CacheBin::~SQLite3CacheBin()
    {
        _ok = false;

        // Clean up this thread's cached statements for this bin.
        // Statements on other threads will be cleaned up by
        // sqlite3_close_v2 deferring the actual close until all
        // statements are finalized (at thread exit).
        t_stmtCache.erase(this);

        if (_ownDb)
        {
            sqlite3_close_v2(_ownDb);
            _ownDb = nullptr;
        }
    }

    ReadResult
    SQLite3CacheBin::read(const std::string& key, const osgDB::Options* dbo, bool isImage)
    {
        if (!_ok)
            return ReadResult(ReadResult::RESULT_NOT_FOUND);

        // Check the async write cache first
        if (_pool)
        {
            ScopedReadLock lock(_writeCacheRWM);
            auto i = _writeCache.find(key);
            if (i != _writeCache.end())
            {
                if (isImage)
                {
                    ReadResult rr(
                        const_cast<osg::Image*>(dynamic_cast<const osg::Image*>(i->second.object.get())),
                        i->second.meta);
                    rr.setLastModifiedTime(DateTime().asTimeStamp());
                    return rr;
                }
                else
                {
                    ReadResult rr(
                        const_cast<osg::Object*>(i->second.object.get()),
                        i->second.meta);
                    rr.setLastModifiedTime(DateTime().asTimeStamp());
                    return rr;
                }
            }
        }

        // Read from database
        std::string data;
        std::string metaJSON;
        TimeStamp timestamp = 0;

        {
            StmtSet& stmts = threadStatements();
            sqlite3_reset(stmts.selectStmt);
            if (_separate)
            {
                sqlite3_bind_text(stmts.selectStmt, 1, key.c_str(), -1, SQLITE_TRANSIENT);
            }
            else
            {
                sqlite3_bind_text(stmts.selectStmt, 1, getID().c_str(), -1, SQLITE_TRANSIENT);
                sqlite3_bind_text(stmts.selectStmt, 2, key.c_str(), -1, SQLITE_TRANSIENT);
            }

            int rc = sqlite3_step(stmts.selectStmt);
            if (rc == SQLITE_ROW)
            {
                const void* blob = sqlite3_column_blob(stmts.selectStmt, 0);
                int blobSize = sqlite3_column_bytes(stmts.selectStmt, 0);
                if (blob && blobSize > 0)
                    data.assign(static_cast<const char*>(blob), blobSize);

                const char* meta = (const char*)sqlite3_column_text(stmts.selectStmt, 1);
                if (meta)
                    metaJSON = meta;

                timestamp = (TimeStamp)sqlite3_column_int64(stmts.selectStmt, 2);
            }
            sqlite3_reset(stmts.selectStmt);

            if (rc != SQLITE_ROW)
            {
                return ReadResult(ReadResult::RESULT_NOT_FOUND);
            }
        }

        if (data.empty())
            return ReadResult(ReadResult::RESULT_NOT_FOUND);

        // Deserialize the blob
        std::istringstream datastream(data);
        osgDB::ReaderWriter::ReadResult r;

        if (isImage)
            r = _rw->readImage(datastream, _rwOptions.get());
        else
            r = _rw->readObject(datastream, _rwOptions.get());

        if (!r.success())
        {
            OE_WARN << LC << "Failed to deserialize cached object for key \"" << key
                << "\" in bin [" << getID() << "]: " << r.message() << std::endl;
            return ReadResult(ReadResult::RESULT_READER_ERROR);
        }

        Config meta;
        if (!metaJSON.empty())
            meta.fromJSON(metaJSON);

        ReadResult rr(r.getObject(), meta);
        rr.setLastModifiedTime(timestamp);
        return rr;
    }

    ReadResult
    SQLite3CacheBin::readImage(const std::string& key, const osgDB::Options* dbo)
    {
        return read(key, dbo, true);
    }

    ReadResult
    SQLite3CacheBin::readObject(const std::string& key, const osgDB::Options* dbo)
    {
        return read(key, dbo, false);
    }

    ReadResult
    SQLite3CacheBin::readString(const std::string& key, const osgDB::Options* dbo)
    {
        ReadResult r = readObject(key, dbo);
        if (r.succeeded())
        {
            if (r.get<StringObject>())
                return r;
            else
                return ReadResult("Empty string");
        }
        return r;
    }

    bool
    SQLite3CacheBin::write(
        const std::string& key,
        const osg::Object* raw_object,
        const Config& meta,
        const osgDB::Options* raw_writeOptions)
    {
        if (!_ok || !raw_object)
            return false;

        bool isNode = dynamic_cast<const osg::Node*>(raw_object) != nullptr;
        if (isNode && _options.enableNodeCaching() == false)
            return true;

        osg::ref_ptr<const osg::Object> object(raw_object);
        osg::ref_ptr<const osgDB::Options> writeOptions(_rwOptions);
        Config metadata(meta);
        std::string binID = getID();

        auto write_op = [=]()
        {
            OE_PROFILING_ZONE_NAMED("OE SQLite3 Cache Write");

            // Serialize the object to a stringstream
            osgDB::ReaderWriter::WriteResult r;
            std::stringstream datastream;

            if (dynamic_cast<const osg::Image*>(object.get()))
            {
                r = _rw->writeImage(
                    *static_cast<const osg::Image*>(object.get()),
                    datastream, writeOptions.get());
            }
            else if (dynamic_cast<const osg::Node*>(object.get()))
            {
                r = _rw->writeNode(
                    *static_cast<const osg::Node*>(object.get()),
                    datastream, writeOptions.get());
            }
            else
            {
                r = _rw->writeObject(*object.get(), datastream, writeOptions.get());
            }

            if (!r.success())
            {
                OE_WARN << LC << "FAILED to serialize object for key \"" << key
                    << "\" in bin [" << binID << "]: " << r.message() << std::endl;
            }
            else
            {
                std::string data = datastream.str();
                std::string metaJSON = metadata.toJSON();

                // Retry on SQLITE_BUSY/SQLITE_LOCKED, which can occur under
                // multi-process contention even with a busy handler installed.
                StmtSet& stmts = threadStatements();
                int rc;
                int tries = 0;
                do {
                    sqlite3_reset(stmts.insertStmt);
                    if (_separate)
                    {
                        sqlite3_bind_text(stmts.insertStmt, 1, key.c_str(), -1, SQLITE_TRANSIENT);
                        sqlite3_bind_blob(stmts.insertStmt, 2, data.data(), (int)data.size(), SQLITE_TRANSIENT);
                        sqlite3_bind_text(stmts.insertStmt, 3, metaJSON.c_str(), -1, SQLITE_TRANSIENT);
                    }
                    else
                    {
                        sqlite3_bind_text(stmts.insertStmt, 1, binID.c_str(), -1, SQLITE_TRANSIENT);
                        sqlite3_bind_text(stmts.insertStmt, 2, key.c_str(), -1, SQLITE_TRANSIENT);
                        sqlite3_bind_blob(stmts.insertStmt, 3, data.data(), (int)data.size(), SQLITE_TRANSIENT);
                        sqlite3_bind_text(stmts.insertStmt, 4, metaJSON.c_str(), -1, SQLITE_TRANSIENT);
                    }
                    rc = sqlite3_step(stmts.insertStmt);
                    sqlite3_reset(stmts.insertStmt);
                }
                while (++tries < 100 && (rc == SQLITE_BUSY || rc == SQLITE_LOCKED));

                if (rc != SQLITE_DONE)
                {
                    OE_WARN << LC << "FAILED to write key \"" << key << "\" to bin ["
                        << binID << "]: " << sqlite3_errmsg(db()) << std::endl;
                }
            }

            // Remove from write cache
            {
                ScopedWriteLock lock(_writeCacheRWM);
                _writeCache.erase(key);
            }
        };

        if (_pool != nullptr)
        {
            // Store in write cache for read-back before async write completes
            _writeCacheRWM.lock();
            WriteCacheRecord& record = _writeCache[key];
            record.meta = meta;
            record.object = object;
            _writeCacheRWM.unlock();

            jobs::dispatch(write_op, jobs::context{ key, _pool });
        }
        else
        {
            write_op();
        }

        return true;
    }

    CacheBin::RecordStatus
    SQLite3CacheBin::getRecordStatus(const std::string& key)
    {
        if (!_ok) return STATUS_NOT_FOUND;

        // Check write cache first
        if (_pool)
        {
            ScopedReadLock lock(_writeCacheRWM);
            if (_writeCache.find(key) != _writeCache.end())
                return STATUS_OK;
        }

        StmtSet& stmts = threadStatements();
        sqlite3_reset(stmts.existsStmt);
        if (_separate)
        {
            sqlite3_bind_text(stmts.existsStmt, 1, key.c_str(), -1, SQLITE_TRANSIENT);
        }
        else
        {
            sqlite3_bind_text(stmts.existsStmt, 1, getID().c_str(), -1, SQLITE_TRANSIENT);
            sqlite3_bind_text(stmts.existsStmt, 2, key.c_str(), -1, SQLITE_TRANSIENT);
        }

        int rc = sqlite3_step(stmts.existsStmt);
        sqlite3_reset(stmts.existsStmt);
        return (rc == SQLITE_ROW) ? STATUS_OK : STATUS_NOT_FOUND;
    }

    bool
    SQLite3CacheBin::remove(const std::string& key)
    {
        if (!_ok) return false;

        // Remove from write cache
        {
            ScopedWriteLock lock(_writeCacheRWM);
            _writeCache.erase(key);
        }

        StmtSet& stmts = threadStatements();
        int rc, tries = 0;
        do {
            sqlite3_reset(stmts.deleteStmt);
            if (_separate)
            {
                sqlite3_bind_text(stmts.deleteStmt, 1, key.c_str(), -1, SQLITE_TRANSIENT);
            }
            else
            {
                sqlite3_bind_text(stmts.deleteStmt, 1, getID().c_str(), -1, SQLITE_TRANSIENT);
                sqlite3_bind_text(stmts.deleteStmt, 2, key.c_str(), -1, SQLITE_TRANSIENT);
            }
            rc = sqlite3_step(stmts.deleteStmt);
            sqlite3_reset(stmts.deleteStmt);
        }
        while (++tries < 100 && (rc == SQLITE_BUSY || rc == SQLITE_LOCKED));
        return rc == SQLITE_DONE;
    }

    bool
    SQLite3CacheBin::touch(const std::string& key)
    {
        if (!_ok) return false;

        StmtSet& stmts = threadStatements();
        int rc, tries = 0;
        do {
            sqlite3_reset(stmts.touchStmt);
            if (_separate)
            {
                sqlite3_bind_text(stmts.touchStmt, 1, key.c_str(), -1, SQLITE_TRANSIENT);
            }
            else
            {
                sqlite3_bind_text(stmts.touchStmt, 1, getID().c_str(), -1, SQLITE_TRANSIENT);
                sqlite3_bind_text(stmts.touchStmt, 2, key.c_str(), -1, SQLITE_TRANSIENT);
            }
            rc = sqlite3_step(stmts.touchStmt);
            sqlite3_reset(stmts.touchStmt);
        }
        while (++tries < 100 && (rc == SQLITE_BUSY || rc == SQLITE_LOCKED));
        return rc == SQLITE_DONE;
    }

    bool
    SQLite3CacheBin::clear()
    {
        if (!_ok) return false;

        // Clear write cache
        {
            ScopedWriteLock lock(_writeCacheRWM);
            _writeCache.clear();
        }

        StmtSet& stmts = threadStatements();
        int rc, tries = 0;
        do {
            sqlite3_reset(stmts.clearStmt);
            if (!_separate)
            {
                sqlite3_bind_text(stmts.clearStmt, 1, getID().c_str(), -1, SQLITE_TRANSIENT);
            }
            rc = sqlite3_step(stmts.clearStmt);
            sqlite3_reset(stmts.clearStmt);
        }
        while (++tries < 100 && (rc == SQLITE_BUSY || rc == SQLITE_LOCKED));
        return rc == SQLITE_DONE;
    }

    bool
    SQLite3CacheBin::compact()
    {
        if (!_ok) return false;

        std::lock_guard<std::mutex> lock(_compactMutex);
        char* errMsg = nullptr;
        int rc = sqlite3_exec(db(), "PRAGMA optimize; VACUUM;", nullptr, nullptr, &errMsg);
        if (rc != SQLITE_OK)
        {
            OE_WARN << LC << "Compact failed: " << (errMsg ? errMsg : "unknown") << std::endl;
            sqlite3_free(errMsg);
            return false;
        }
        return true;
    }

    unsigned
    SQLite3CacheBin::getStorageSize()
    {
        if (!_ok) return 0u;

        StmtSet& stmts = threadStatements();
        sqlite3_reset(stmts.sizeStmt);
        if (!_separate)
        {
            sqlite3_bind_text(stmts.sizeStmt, 1, getID().c_str(), -1, SQLITE_TRANSIENT);
        }

        int rc = sqlite3_step(stmts.sizeStmt);
        unsigned result = 0u;
        if (rc == SQLITE_ROW)
            result = (unsigned)sqlite3_column_int64(stmts.sizeStmt, 0);
        sqlite3_reset(stmts.sizeStmt);
        return result;
    }

} // anonymous namespace

// ---------------------------------------------------------------
// Driver registration
// ---------------------------------------------------------------

class SQLite3CacheDriver : public CacheDriver
{
public:
    SQLite3CacheDriver()
    {
        supportsExtension( "osgearth_cache_sqlite3", "SQLite3 cache for osgEarth" );
    }

    virtual const char* className() const
    {
        return "SQLite3 cache for osgEarth";
    }

    virtual ReadResult readObject(const std::string& file_name, const Options* options) const
    {
        if ( !acceptsExtension(osgDB::getLowerCaseFileExtension( file_name )))
            return ReadResult::FILE_NOT_HANDLED;

        return ReadResult( new SQLite3Cache( getCacheOptions(options) ) );
    }
};

REGISTER_OSGPLUGIN(osgearth_cache_sqlite3, SQLite3CacheDriver)
