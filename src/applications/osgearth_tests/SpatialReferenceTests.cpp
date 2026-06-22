/* osgEarth
* Copyright 2025 Pelican Mapping
* MIT License
*/

#include <osgEarth/catch.hpp>
#include <cmath>
#include <osgEarth/SpatialReference>
#include <osgEarth/Profile>

using namespace osgEarth;

namespace
{
    template<typename T>
    bool vec_eq(const T& a, const T& b) {
        for(int i=0; i<a.num_components; ++i)
            if (!osg::equivalent(a[i], b[i])) return false;
        return true;
    }

    bool near(double a, double b, double eps = 1e-3)
    {
        return std::abs(a - b) <= eps;
    }
                 
}

TEST_CASE("Parsing doubles") {
    REQUIRE(Strings::parseDouble("123.4") == 123.4);
    REQUIRE(Strings::parseDouble(" 123.4") == 123.4);
    REQUIRE(Strings::parseDouble(" 123.4 ") == 123.4);
    REQUIRE(Strings::parseDouble("123.4abc") == 123.4);
    REQUIRE(std::isnan(Strings::parseDouble("abc")));

    REQUIRE(Strings::parseDoubleAndIndex("123.4m") == std::make_pair(123.4, 5));
}

TEST_CASE( "SpatialReferences are cached" ) {
    osg::ref_ptr< const SpatialReference > srs1 = SpatialReference::create("spherical-mercator");
    REQUIRE(srs1.valid());

    osg::ref_ptr< const SpatialReference > srs2 = SpatialReference::create("spherical-mercator");
    REQUIRE(srs2.valid());

    REQUIRE(srs1.get() == srs2.get());
}

TEST_CASE( "Spherical Mercator SpatialReferences can be created" ) {
    // Create a spherical mercator SRS.
    osg::ref_ptr< const SpatialReference > mercSRS = SpatialReference::create("spherical-mercator");

    REQUIRE( mercSRS.valid() );
    REQUIRE( mercSRS->isSphericalMercator() );
    REQUIRE( mercSRS->isProjected() );
    REQUIRE( !mercSRS->isGeodetic() );
    REQUIRE( !mercSRS->isGeographic() );

    SECTION("epsg:900913 is equivalent") {
        osg::ref_ptr< const SpatialReference > epsg900913 = SpatialReference::create("epsg:900913");
        REQUIRE(epsg900913.valid());
        REQUIRE(epsg900913->isEquivalentTo(mercSRS.get()));
    }

    SECTION("epsg:3785 is equivalent") {
        osg::ref_ptr< const SpatialReference > epsg3785 = SpatialReference::create("epsg:3785");
        REQUIRE(epsg3785.valid());
        REQUIRE(epsg3785->isEquivalentTo(mercSRS.get()));
    }

    SECTION("epsg:102113 is equivalent") {
        osg::ref_ptr< const SpatialReference > epsg102113 = SpatialReference::create("epsg:102113");
        REQUIRE(epsg102113.valid());
        REQUIRE(epsg102113->isEquivalentTo(mercSRS.get()));
    }
}

TEST_CASE( "WGS84 SpatialReferences can be created" ) {
    // Create a wgs84 mercator SRS.
    osg::ref_ptr< const SpatialReference > wgs84 = SpatialReference::create("wgs84");
    REQUIRE( wgs84.valid() );

    REQUIRE(wgs84->isGeographic());
    REQUIRE(wgs84->isGeodetic());
    REQUIRE(!wgs84->isMercator());
    REQUIRE(!wgs84->isProjected());

    SECTION("epsg:4326 is equivalent") {
        osg::ref_ptr< const SpatialReference > epsg4326 = SpatialReference::create("epsg:4326");
        REQUIRE(epsg4326.valid());
        REQUIRE(epsg4326->isEquivalentTo(wgs84.get()));
    }

    SECTION("WGS84 and Mercator are not equivalent") {
        osg::ref_ptr<const SpatialReference> sm = SpatialReference::get("spherical-mercator");
        REQUIRE(!wgs84->isEquivalentTo(sm.get()));
    }
}

TEST_CASE("Plate Carree SpatialReferences can be created") {
    osg::ref_ptr< const SpatialReference > plateCarre = SpatialReference::create("plate-carree");
    REQUIRE(plateCarre.valid());
    REQUIRE(!plateCarre->isGeographic());
    REQUIRE(!plateCarre->isMercator());
    REQUIRE(!plateCarre->isGeodetic());
    REQUIRE(plateCarre->isProjected());
}

TEST_CASE("PC/WGS84 transform") {
    const SpatialReference* wgs84 = SpatialReference::get("wgs84");
    const SpatialReference* pceqc = SpatialReference::get("plate-carree");

    double pcMinX = -20037508.3427892476320267;
    double pcMaxX = 20037508.3427892476320267;
    double pcMinY = -10018754.1713946219533682;
    double pcMaxY = 10018754.1713946219533682;

    osg::Vec3d output;
    REQUIRE(wgs84->transform(osg::Vec3d(-180, -90, 0), pceqc, output));
    REQUIRE(vec_eq(output, osg::Vec3d(pcMinX, pcMinY, 0)));

    REQUIRE(wgs84->transform(osg::Vec3d(-180, +90, 0), pceqc, output));
    REQUIRE(vec_eq(output, osg::Vec3d(pcMinX, pcMaxY, 0)));

    REQUIRE(wgs84->transform(osg::Vec3d(180, -90, 0), pceqc, output));
    REQUIRE(vec_eq(output, osg::Vec3d(pcMaxX, pcMinY, 0)));

    REQUIRE(wgs84->transform(osg::Vec3d(180, 90, 0), pceqc, output));
    REQUIRE(vec_eq(output, osg::Vec3d(pcMaxX, pcMaxY, 0)));
}

TEST_CASE("SphMercator/WGS84 transform") {
    const SpatialReference* wgs84 = SpatialReference::get("wgs84");
    const SpatialReference* sm = SpatialReference::get("spherical-mercator");
    osg::Vec3d temp1, temp2;

    // illegal input:
    //REQUIRE(wgs84->transform(osg::Vec3d(-180, -90, 0), sm, temp1) == false);

    // valid input, there and back:
    REQUIRE(wgs84->transform(osg::Vec3d(-180, -85, 0), sm, temp1) == true);
    REQUIRE(sm->transform(temp1, wgs84, temp2) == true);
    REQUIRE(vec_eq(temp2, osg::Vec3d(-180, -85, 0)));
}

TEST_CASE("Vertical Datum Tests") {
    const SpatialReference* wgs84 = SpatialReference::get("wgs84");
    const SpatialReference* wgs84_egm96 = SpatialReference::get("wgs84", "egm96");

    double eps = 0.2;
    // Vertical datum tests.
    // Reference: http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/intpt.html
    osg::Vec3d output;

    REQUIRE(wgs84->transform(osg::Vec3d(0, 0, 17.16), wgs84_egm96, output));
    REQUIRE(osg::equivalent(output.z(), 0.0, eps));

    REQUIRE(wgs84->transform(osg::Vec3d(90, 0, -63.24), wgs84_egm96, output));
    REQUIRE(osg::equivalent(output.z(), 0.0, eps));

    REQUIRE(wgs84->transform(osg::Vec3d(180, 0, 21.15), wgs84_egm96, output));
    REQUIRE(osg::equivalent(output.z(), 0.0, eps));

    REQUIRE(wgs84->transform(osg::Vec3d(-90, 0, -4.29), wgs84_egm96, output));
    REQUIRE(osg::equivalent(output.z(), 0.0, eps));
}

TEST_CASE("getGeographicsSRS") {
    const SpatialReference* mercator = SpatialReference::get("spherical-mercator", "egm96");
    const SpatialReference* geo = mercator->getGeographicSRS();
    REQUIRE(geo->isGeographic());
    REQUIRE(geo->getVerticalDatum() != nullptr); // ditto
    REQUIRE(!geo->isGeodetic()); // vertical datum must be perserved
}

TEST_CASE("getGeodeticSRS") {
    const SpatialReference* mercator = SpatialReference::get("spherical-mercator", "egm96");
    const SpatialReference* geodetic = mercator->getGeodeticSRS();
    REQUIRE(geodetic->isGeographic());
    REQUIRE(geodetic->isGeodetic());
    REQUIRE(geodetic->getVerticalDatum() == nullptr);
}

TEST_CASE("getGeocentricSRS") {
    const SpatialReference* wgs84 = SpatialReference::get("wgs84");
    const SpatialReference* ecef = wgs84->getGeocentricSRS();

    osg::Vec3d np_wgs84, np_ecef, temp;

    np_wgs84.set(0.0, 90.0, 0.0);

    np_ecef.set(0.0, 0.0, wgs84->getEllipsoid().getRadiusPolar());

    REQUIRE(wgs84->transform(np_wgs84, ecef, temp));
    REQUIRE(vec_eq(temp, np_ecef));

    REQUIRE(ecef->transform(np_ecef, wgs84, temp));
    REQUIRE(vec_eq(temp, np_wgs84));
}

TEST_CASE("MFE SRS") {
    const SpatialReference* wgs84 = SpatialReference::get("wgs84");
    const SpatialReference* mfe = SpatialReference::get("+proj=eqc +lat_0=21n +lon_0=157w +lat_ts=21n +ellps=WGS84 +units=m");
    osg::Vec3d p_wgs84, p_mfe;
    p_mfe.set(1, 2, 3);
    p_wgs84.set(-157.0, 21.0, 0.0);

    wgs84->transform(p_wgs84, mfe, p_mfe);
    REQUIRE(p_mfe.x() == 0.0);
    REQUIRE(p_mfe.y() == 0.0);

    p_wgs84.set(0, 0, 0);
    mfe->transform(p_mfe, wgs84, p_wgs84);
    REQUIRE(p_wgs84.x() == -157.0);
    REQUIRE(p_wgs84.y() == 21.0);
}

TEST_CASE("transformExtentToMBR is a no-op for equivalent SRS") {
    const SpatialReference* wgs84 = SpatialReference::get("wgs84");

    double xmin = -10.0;
    double ymin = -20.0;
    double xmax = 30.0;
    double ymax = 40.0;

    REQUIRE(wgs84->transformExtentToMBR(wgs84, xmin, ymin, xmax, ymax));
    REQUIRE(xmin == -10.0);
    REQUIRE(ymin == -20.0);
    REQUIRE(xmax == 30.0);
    REQUIRE(ymax == 40.0);
}

TEST_CASE("transformExtentToMBR clamps WGS84 latitude when targeting spherical mercator") {
    const SpatialReference* wgs84 = SpatialReference::get("wgs84");
    const SpatialReference* sm = SpatialReference::get("spherical-mercator");
    osg::ref_ptr<const Profile> mercProfile = Profile::create(Profile::SPHERICAL_MERCATOR);

    double xmin = -5.0;
    double ymin = -90.0; // outside spherical-mercator legal lat range
    double xmax = 5.0;
    double ymax = 90.0;  // outside spherical-mercator legal lat range

    REQUIRE(wgs84->transformExtentToMBR(sm, xmin, ymin, xmax, ymax));
    REQUIRE(std::isfinite(xmin));
    REQUIRE(std::isfinite(ymin));
    REQUIRE(std::isfinite(xmax));
    REQUIRE(std::isfinite(ymax));

    osg::Vec3d pMin, pMax;
    const double latMin = mercProfile->getLatLongExtent().yMin();
    const double latMax = mercProfile->getLatLongExtent().yMax();

    REQUIRE(wgs84->transform(osg::Vec3d(0.0, latMin, 0.0), sm, pMin));
    REQUIRE(wgs84->transform(osg::Vec3d(0.0, latMax, 0.0), sm, pMax));

    REQUIRE(near(ymin, pMin.y(), 1e-2));
    REQUIRE(near(ymax, pMax.y(), 1e-2));
}

TEST_CASE("transformExtentToMBR rejects inverted extents") {
    const SpatialReference* wgs84 = SpatialReference::get("wgs84");
    const SpatialReference* sm = SpatialReference::get("spherical-mercator");

    double xmin = 10.0;
    double ymin = 5.0;
    double xmax = -10.0; // inverted X range
    double ymax = 15.0;

    REQUIRE_FALSE(wgs84->transformExtentToMBR(sm, xmin, ymin, xmax, ymax));
}
