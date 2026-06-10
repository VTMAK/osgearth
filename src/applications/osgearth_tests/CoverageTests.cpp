/* osgEarth
/* osgEarth
 * Copyright 2025 Pelican Mapping
 * MIT License
 */

#include <osgEarth/catch.hpp>
#include <osgEarth/Coverage>
#include <osgEarth/GeoData>

using namespace osgEarth;

namespace
{
    struct TestCoverageValue
    {
        bool _valid;
        int _value;

        TestCoverageValue() : _valid(false), _value(0) { }
        explicit TestCoverageValue(int value) : _valid(true), _value(value) { }

        explicit TestCoverageValue(const Config& conf) : _valid(false), _value(0)
        {
            conf.get("valid", _valid);
            conf.get("value", _value);
        }

        bool valid() const
        {
            return _valid;
        }

        bool operator < (const TestCoverageValue& rhs) const
        {
            if (_valid != rhs._valid)
                return _valid < rhs._valid;
            return _value < rhs._value;
        }

        bool operator == (const TestCoverageValue& rhs) const
        {
            return _valid == rhs._valid && _value == rhs._value;
        }

        Config getConfig() const
        {
            Config conf("value");
            conf.set("valid", _valid);
            conf.set("value", _value);
            return conf;
        }
    };
}

TEST_CASE("Coverage allocate/write/read and nodata tracking")
{
    Coverage<TestCoverageValue>::Ptr cov = Coverage<TestCoverageValue>::create();
    cov->allocate(4, 4);

    REQUIRE(cov->s() == 4u);
    REQUIRE(cov->t() == 4u);
    REQUIRE(cov->nodataCount() == 16u);
    REQUIRE(cov->empty());

    int k = cov->write(TestCoverageValue(7), 1, 2);
    REQUIRE(k > 0);
    REQUIRE(cov->hasDataAt(1, 2));
    REQUIRE(cov->nodataCount() == 15u);
    REQUIRE_FALSE(cov->empty());

    TestCoverageValue out;
    REQUIRE(cov->read(out, 1u, 2u));
    REQUIRE(out == TestCoverageValue(7));

    REQUIRE(cov->write(TestCoverageValue(), 1, 2) == 0);
    REQUIRE_FALSE(cov->hasDataAt(1, 2));
    REQUIRE(cov->nodataCount() == 16u);
    REQUIRE(cov->empty());
}

TEST_CASE("Coverage deduplicates equal values and supports write hint")
{
    Coverage<TestCoverageValue>::Ptr cov = Coverage<TestCoverageValue>::create();
    cov->allocate(4, 1);

    int k1 = cov->write(TestCoverageValue(5), 0, 0);
    int k2 = cov->write(TestCoverageValue(5), 1, 0);
    int k3 = cov->write(TestCoverageValue(9), 2, 0);
    int k4 = cov->write(TestCoverageValue(5), 3, 0, k1);

    REQUIRE(k1 > 0);
    REQUIRE(k2 == k1);
    REQUIRE(k3 != k1);
    REQUIRE(k4 == k1);

    TestCoverageValue v0, v1, v2, v3;
    REQUIRE(cov->read(v0, 0u, 0u));
    REQUIRE(cov->read(v1, 1u, 0u));
    REQUIRE(cov->read(v2, 2u, 0u));
    REQUIRE(cov->read(v3, 3u, 0u));

    REQUIRE(v0 == TestCoverageValue(5));
    REQUIRE(v1 == TestCoverageValue(5));
    REQUIRE(v2 == TestCoverageValue(9));
    REQUIRE(v3 == TestCoverageValue(5));
}

TEST_CASE("Coverage config round-trip preserves dimensions and pixels")
{
    Coverage<TestCoverageValue>::Ptr src = Coverage<TestCoverageValue>::create();
    src->allocate(3, 2);

    src->write(TestCoverageValue(10), 0, 0);
    src->write(TestCoverageValue(20), 1, 0);
    src->write(TestCoverageValue(10), 2, 0);
    src->write(TestCoverageValue(),    0, 1);
    src->write(TestCoverageValue(20), 1, 1);
    src->write(TestCoverageValue(30), 2, 1);

    Config conf = src->getConfig();

    Coverage<TestCoverageValue>::Ptr dst = Coverage<TestCoverageValue>::create();
    dst->setConfig(conf);

    REQUIRE(dst->s() == src->s());
    REQUIRE(dst->t() == src->t());
    REQUIRE(dst->nodataCount() == src->nodataCount());

    for (unsigned t = 0; t < src->t(); ++t)
    {
        for (unsigned s = 0; s < src->s(); ++s)
        {
            TestCoverageValue a, b;
            bool hasA = src->read(a, s, t);
            bool hasB = dst->read(b, s, t);

            REQUIRE(hasA == hasB);
            if (hasA)
                REQUIRE(a == b);
        }
    }
}

TEST_CASE("GeoCoverage readAtCoords maps world coords to raster values")
{
    Coverage<TestCoverageValue>::Ptr cov = Coverage<TestCoverageValue>::create();
    cov->allocate(4, 4);

    cov->write(TestCoverageValue(11), 1, 1);
    cov->write(TestCoverageValue(22), 2, 2);

    const SpatialReference* wgs84 = SpatialReference::get("wgs84");
    REQUIRE(wgs84 != 0);

    GeoExtent extent(wgs84, 0.0, 0.0, 40.0, 40.0);
    GeoCoverage<TestCoverageValue> geo(cov, extent);

    TestCoverageValue out;

    // (15,15) -> xs=floor((4-1)*(15/40))=1, yt=1
    REQUIRE(geo.readAtCoords(out, 15.0, 15.0));
    REQUIRE(out == TestCoverageValue(11));

    // (30,30) -> xs=floor((4-1)*(30/40))=2, yt=2
    REQUIRE(geo.readAtCoords(out, 30.0, 30.0));
    REQUIRE(out == TestCoverageValue(22));

    // Outside extent should fail.
    REQUIRE_FALSE(geo.readAtCoords(out, 50.0, 50.0));
}