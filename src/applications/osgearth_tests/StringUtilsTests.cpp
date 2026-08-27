/* osgEarth
 * Copyright 2025 Pelican Mapping
 * MIT License
 */

#include <osgEarth/catch.hpp>
#include <osgEarth/StringUtils>

#include <cmath>
#include <limits>
#include <sstream>

using namespace osgEarth;
using namespace osgEarth::Util;

namespace
{
    bool near(double lhs, double rhs, double epsilon = 1e-6)
    {
        return std::fabs(lhs - rhs) <= epsilon;
    }
}

TEST_CASE("StringTokenizer handles delimiters, quotes, trimming, and empty tokens")
{
    SECTION("standard quotes protect delimiters and remain in tokens")
    {
        const auto tokens = StringTokenizer()
            .delim(",")
            .standardQuotes()
            .tokenize("  alpha , ,\"beta,gamma\",'delta,epsilon', ");

        REQUIRE((tokens == std::vector<std::string>{
            "alpha", "", "\"beta,gamma\"", "'delta,epsilon'" }));
    }

    SECTION("quotes and asymmetric quote pairs can be removed")
    {
        const auto quoted = StringTokenizer()
            .delim(",")
            .quote('"', false)
            .tokenize("alpha,\"beta,gamma\",delta");
        REQUIRE((quoted == std::vector<std::string>{ "alpha", "beta,gamma", "delta" }));

        const auto paired = StringTokenizer()
            .delim(",")
            .quotePair('{', '}', false)
            .tokenize("alpha,{beta,gamma},delta");
        REQUIRE((paired == std::vector<std::string>{ "alpha", "beta,gamma", "delta" }));
    }

    SECTION("multi-character delimiters may be retained")
    {
        const auto tokens = StringTokenizer()
            .delim("::", true)
            .tokenize("alpha::beta::::gamma");
        REQUIRE((tokens == std::vector<std::string>{
            "alpha", "::", "beta", "::", "", "::", "gamma" }));
    }

    SECTION("empty-token and trimming policies are configurable")
    {
        const auto withoutEmpties = StringTokenizer()
            .delim(",")
            .keepEmpties(false)
            .tokenize("alpha,,,beta,");
        REQUIRE((withoutEmpties == std::vector<std::string>{ "alpha", "beta" }));

        const auto untrimmed = StringTokenizer()
            .delim(",")
            .trimTokens(false)
            .tokenize(" alpha , beta ");
        REQUIRE((untrimmed == std::vector<std::string>{ " alpha ", " beta " }));

        const auto empties = StringTokenizer().delim(",").tokenize(",,");
        REQUIRE((empties == std::vector<std::string>{ "", "" }));
    }

    SECTION("standard whitespace delimiters split mixed whitespace")
    {
        const auto tokens = StringTokenizer()
            .whitespaceDelims()
            .keepEmpties(false)
            .tokenize("alpha\tbeta\n gamma\rdelta");
        REQUIRE((tokens == std::vector<std::string>{ "alpha", "beta", "gamma", "delta" }));
    }

    SECTION("empty and delimiter-free input are well-defined")
    {
        REQUIRE(StringTokenizer().tokenize("").empty());
        REQUIRE((StringTokenizer().tokenize("  value  ") == std::vector<std::string>{ "value" }));
    }
}

TEST_CASE("StringTokenizer reports and optionally ignores dangling quotes")
{
    bool error = false;
    auto tokens = StringTokenizer()
        .delim(",")
        .quote('"', false)
        ("alpha,\"beta", &error);

    REQUIRE(error);
    REQUIRE((tokens == std::vector<std::string>{ "alpha", "beta" }));

    error = true;
    tokens = StringTokenizer()
        .delim(",")
        .quote('"', false)
        .ignoreDanglingQuotes()
        ("alpha,\"beta", &error);

    REQUIRE_FALSE(error);
    REQUIRE((tokens == std::vector<std::string>{ "alpha", "beta" }));
}

TEST_CASE("String integral conversion handles bases, signs, limits, and failures")
{
    REQUIRE(as<int>("42", -1) == 42);
    REQUIRE(as<int>("  +42suffix", -1) == 42);
    REQUIRE(as<int>("0x20", -1) == 32);
    REQUIRE(as<int>("010", -1) == 8);
    REQUIRE(as<int>(std::to_string((std::numeric_limits<int>::min)()), 0) ==
        (std::numeric_limits<int>::min)());
    REQUIRE(as<int>(std::to_string((std::numeric_limits<int>::max)()), 0) ==
        (std::numeric_limits<int>::max)());

    REQUIRE(as<int>("invalid", -7) == -7);
    REQUIRE(as<int>("", -7) == -7);
    REQUIRE(as<int>("999999999999999999999", -7) == -7);
    REQUIRE(as<unsigned>("-1", 7u) == 7u);

    REQUIRE(static_cast<int>(as<signed char>("-128", 0)) == -128);
    REQUIRE(static_cast<int>(as<signed char>("127", 0)) == 127);
    REQUIRE(static_cast<int>(as<signed char>("128", 9)) == 9);
    REQUIRE(static_cast<unsigned>(as<unsigned char>("255", 0u)) == 255u);
    REQUIRE(static_cast<unsigned>(as<unsigned char>("256", 9u)) == 9u);

    REQUIRE(as<short>("32768", 11) == 11);
    REQUIRE(as<unsigned short>("65536", 12u) == 12u);
    REQUIRE(as<long>(std::to_string((std::numeric_limits<long>::max)()), 0L) ==
        (std::numeric_limits<long>::max)());
    REQUIRE(as<unsigned long>(std::to_string((std::numeric_limits<unsigned long>::max)()), 0UL) ==
        (std::numeric_limits<unsigned long>::max)());
}

TEST_CASE("String floating-point and boolean conversion preserves established forms")
{
    REQUIRE(near(as<float>(" 12.5suffix", -1.0f), 12.5));
    REQUIRE(near(as<double>("1.25e3", -1.0), 1250.0));
    REQUIRE(near(static_cast<double>(as<long double>("-2.5e-3", -1.0L)), -0.0025));
    REQUIRE(as<double>("invalid", -3.0) == -3.0);
    REQUIRE(as<double>("1e9999", -3.0) == -3.0);

    REQUIRE(as<bool>("true", false));
    REQUIRE(as<bool>("YES", false));
    REQUIRE(as<bool>("On", false));
    REQUIRE_FALSE(as<bool>("false", true));
    REQUIRE_FALSE(as<bool>("NO", true));
    REQUIRE_FALSE(as<bool>("off", true));
    REQUIRE(as<bool>("unknown", true));
    REQUIRE_FALSE(as<bool>(" unknown ", false));

    REQUIRE(as<std::string>("value", "fallback") == "value");

    const osg::Vec3f fallback(9.0f, 9.0f, 9.0f);
    const osg::Vec3f vector = as<osg::Vec3f>("1 2 3", fallback);
    REQUIRE(near(vector.x(), 1.0));
    REQUIRE(near(vector.y(), 2.0));
    REQUIRE(near(vector.z(), 3.0));
}

TEST_CASE("StringCursor reads whitespace-separated primitive and generic values")
{
    const std::string input = "  -12\t255  1.25e2 word 1 \n";
    StringCursor cursor(input);

    int integer = 0;
    unsigned char byte = 0u;
    double floatingPoint = 0.0;
    std::string word;
    bool boolean = false;

    REQUIRE_FALSE(cursor.atEnd());
    REQUIRE(cursor.read(integer));
    REQUIRE(integer == -12);
    REQUIRE(cursor.read(byte));
    REQUIRE(static_cast<unsigned>(byte) == 255u);
    REQUIRE(cursor.read(floatingPoint));
    REQUIRE(floatingPoint == 125.0);
    REQUIRE(cursor.read(word));
    REQUIRE(word == "word");
    REQUIRE(cursor.read(boolean));
    REQUIRE(boolean);
    REQUIRE(cursor.atEnd());
    REQUIRE_FALSE(cursor.read(integer));

    const std::string boundaries =
        std::to_string((std::numeric_limits<int>::min)()) + " " +
        std::to_string((std::numeric_limits<unsigned int>::max)());
    StringCursor boundaryCursor(boundaries);
    unsigned int unsignedInteger = 0u;
    REQUIRE(boundaryCursor.read(integer));
    REQUIRE(integer == (std::numeric_limits<int>::min)());
    REQUIRE(boundaryCursor.read(unsignedInteger));
    REQUIRE(unsignedInteger == (std::numeric_limits<unsigned int>::max)());
    REQUIRE(boundaryCursor.atEnd());
}

TEST_CASE("StringCursor failures do not modify outputs or advance past bad tokens")
{
    int integer = 7;
    const std::string overflow = "999999999999999999999 3";
    StringCursor overflowCursor(overflow);
    REQUIRE_FALSE(overflowCursor.read(integer));
    REQUIRE(integer == 7);
    REQUIRE_FALSE(overflowCursor.read(integer));
    REQUIRE(integer == 7);

    unsigned unsignedInteger = 9u;
    const std::string negative = "-1";
    StringCursor negativeCursor(negative);
    REQUIRE_FALSE(negativeCursor.read(unsignedInteger));
    REQUIRE(unsignedInteger == 9u);

    double floatingPoint = 4.0;
    const std::string floatingOverflow = "1e9999";
    StringCursor floatingOverflowCursor(floatingOverflow);
    REQUIRE_FALSE(floatingOverflowCursor.read(floatingPoint));
    REQUIRE(floatingPoint == 4.0);

    const std::string whitespace = " \t\r\n";
    StringCursor whitespaceCursor(whitespace);
    REQUIRE(whitespaceCursor.atEnd());
    REQUIRE_FALSE(whitespaceCursor.read(integer));
}

TEST_CASE("String replacement utilities replace all matches safely")
{
    std::string value = "one fish, two fish, fish";
    REQUIRE(&replaceIn(value, "fish", "cat") == &value);
    REQUIRE(value == "one cat, two cat, cat");

    value = "aaaa";
    replaceIn(value, "aa", "b");
    REQUIRE(value == "bb");

    const std::string unchanged = value;
    replaceIn(value, "", "ignored");
    REQUIRE(value == unchanged);

    value = "Alpha aLPHa ALPHABET";
    REQUIRE(&ciReplaceIn(value, "alpha", "x") == &value);
    REQUIRE(value == "x x xBET");

    ciReplaceIn(value, "", "ignored");
    REQUIRE(value == "x x xBET");
}

TEST_CASE("String whitespace and case helpers cover empty and mixed input")
{
    REQUIRE(trim(" \t\r\n value \f\v") == "value");
    REQUIRE(trim(" \t\r\n").empty());
    REQUIRE(trim("").empty());

    std::string inPlace = " \t value \n";
    trim2(inPlace);
    REQUIRE(inPlace == "value");
    inPlace = " \t\n";
    trim2(inPlace);
    REQUIRE(inPlace.empty());

    REQUIRE(trimAndCompress(" \t alpha \n beta\r\n gamma ") == "alpha beta gamma");
    REQUIRE(trimAndCompress(" \t\r\n").empty());

    REQUIRE(toLower('A') == 'a');
    REQUIRE(toLower('7') == '7');
    REQUIRE(toLower("AbC-123") == "abc-123");
    std::string lowerInPlace = "OsG_EaRtH";
    REQUIRE(&toLowerInPlace(lowerInPlace) == &lowerInPlace);
    REQUIRE(lowerInPlace == "osg_earth");
}

TEST_CASE("String prefix, suffix, containment, and equality helpers honor case policy")
{
    REQUIRE(startsWith("osgEarth", "osg"));
    REQUIRE_FALSE(startsWith("osgEarth", "OSG"));
    REQUIRE(startsWith("osgEarth", ""));
    REQUIRE_FALSE(startsWith("osg", "osgEarth"));

    REQUIRE(ciStartsWith("osgEarth", "OSG"));
    REQUIRE(endsWith("map.earth", ".earth"));
    REQUIRE_FALSE(endsWith("map.EARTH", ".earth"));
    REQUIRE(ciEndsWith("map.EARTH", ".earth"));

    REQUIRE(contains("osgEarth", "Earth"));
    REQUIRE_FALSE(contains("osgEarth", "earth"));
    REQUIRE(ciContains("osgEarth", "EARTH"));
    REQUIRE(ciContains("osgEarth", ""));
    REQUIRE_FALSE(ciContains("short", "longer"));

    REQUIRE(ci_equals("MiXeD", "mixed"));
    REQUIRE(ciEquals("MiXeD", "MIXED"));
    REQUIRE_FALSE(ci_equals("mixed", "mismatch"));
}

TEST_CASE("String number validation and parsing report consumed input")
{
    const auto parsed = parseDoubleAndIndex("  -2.5suffix");
    REQUIRE(near(parsed.first, -2.5));
    REQUIRE(parsed.second == 6);

    const auto invalid = parseDoubleAndIndex("invalid");
    REQUIRE(std::isnan(invalid.first));
    REQUIRE(invalid.second == 0);
    REQUIRE(std::isnan(parseDouble("")));
    REQUIRE(std::isnan(parseDouble("1e9999")));
    REQUIRE(near(parseDouble("1.25e2suffix"), 125.0));

    const auto validNumber = isValidNumber("-1.25e2");
    REQUIRE(validNumber.first);
    REQUIRE(near(validNumber.second, -125.0));
    REQUIRE_FALSE(isValidNumber("12suffix").first);
    REQUIRE_FALSE(isValidNumber("").first);
}

TEST_CASE("String joining, token extraction, and unquoting handle edge cases")
{
    REQUIRE(joinStrings({}, ',').empty());
    REQUIRE(joinStrings({ "one" }, ',') == "one");
    REQUIRE(joinStrings({ "one", "", "three" }, ',') == "one,,three");

    REQUIRE(getToken("alpha,\"beta,gamma\",delta", 0u, ',') == "alpha");
    REQUIRE(getToken("alpha,\"beta,gamma\",delta", 1u, ',') == "\"beta,gamma\"");
    REQUIRE(getToken("alpha,\"beta,gamma\",delta", 3u, ',').empty());

    REQUIRE(unquote("  'value'  ") == "value");
    REQUIRE(unquote("\"'nested'\"") == "nested");
    REQUIRE(unquote("'mismatched\"") == "'mismatched\"");
    REQUIRE(unquote("x") == "x");
}

TEST_CASE("String color and vector conversions preserve values and fallbacks")
{
    const osg::Vec4ub fallbackColor(1u, 2u, 3u, 4u);
    const osg::Vec4ub color = stringToColor("10 20 30 40", fallbackColor);
    REQUIRE(static_cast<unsigned>(color.r()) == 10u);
    REQUIRE(static_cast<unsigned>(color.g()) == 20u);
    REQUIRE(static_cast<unsigned>(color.b()) == 30u);
    REQUIRE(static_cast<unsigned>(color.a()) == 40u);
    REQUIRE(colorToString(color) == "10 20 30 40");

    const osg::Vec4ub invalidColor = stringToColor("10 20 30", fallbackColor);
    REQUIRE(invalidColor == fallbackColor);

    const osg::Vec3f fallbackVector(9.0f, 8.0f, 7.0f);
    const osg::Vec3f vector = stringToVec3f("1 2 3", fallbackVector);
    REQUIRE(near(vector.x(), 1.0));
    REQUIRE(near(vector.y(), 2.0));
    REQUIRE(near(vector.z(), 3.0));

    const osg::Vec3f replicated = stringToVec3f("4", fallbackVector);
    REQUIRE(near(replicated.x(), 4.0));
    REQUIRE(near(replicated.y(), 4.0));
    REQUIRE(near(replicated.z(), 4.0));
    REQUIRE(stringToVec3f("1 2", fallbackVector) == fallbackVector);
    REQUIRE(stringToVec3f("invalid", fallbackVector) == fallbackVector);
    REQUIRE(vec3fToString(osg::Vec3f(1.0f, 2.0f, 3.0f)) == "1 2 3\n");
}

TEST_CASE("String HTML color conversion handles RGB and RGBA")
{
    const osg::Vec4f red = htmlColorToVec4f("#ff0000");
    REQUIRE(near(red.r(), 1.0));
    REQUIRE(near(red.g(), 0.0));
    REQUIRE(near(red.b(), 0.0));
    REQUIRE(near(red.a(), 1.0));

    const osg::Vec4f rgba = htmlColorToVec4f("#33669980");
    REQUIRE(near(rgba.r(), 0x33 / 255.0));
    REQUIRE(near(rgba.g(), 0x66 / 255.0));
    REQUIRE(near(rgba.b(), 0x99 / 255.0));
    REQUIRE(near(rgba.a(), 0x80 / 255.0));
    REQUIRE(vec4fToHtmlColor(rgba) == "#33669980");
    REQUIRE(vec4fToHtmlColor(osg::Vec4f(1.0f, 0.0f, 0.0f, 1.0f)) == "#ff0000");

    const osg::Vec4f invalid = htmlColorToVec4f("invalid");
    REQUIRE(near(invalid.r(), 0.0));
    REQUIRE(near(invalid.g(), 0.0));
    REQUIRE(near(invalid.b(), 0.0));
    REQUIRE(near(invalid.a(), 1.0));
}

TEST_CASE("String filename, hash, and display helpers are stable")
{
    REQUIRE(toLegalFileName("https://host/path name") == "host-2f-path-20-name");
    REQUIRE(toLegalFileName("https://host/path name", true) == "host/path-20-name");
    REQUIRE(toLegalFileName("https://host/path name", false, "_") == "host_path_name");
    REQUIRE(toLegalFileName("alpha_beta.txt") == "alpha_beta.txt");

    REQUIRE(hashString("") == 0u);
    REQUIRE(hashString("osgEarth") == 0x10215ec679fdcc76ULL);
    REQUIRE(hashString("abcdefghijklmnop") == 0x87f9ea36917d5c39ULL);
    REQUIRE(hashToString("") == "00000000");
    REQUIRE(hashToString("osgEarth") == "10215ec679fdcc76");

    REQUIRE(prettyPrintTime(3661.5) == "1:1:1.5");
    REQUIRE(prettyPrintSize(512.0) == "512 MB");
    REQUIRE(prettyPrintSize(2048.0) == "2 GB");
    REQUIRE(prettyPrintSize(2.0 * 1024.0 * 1024.0) == "2 TB");
}

TEST_CASE("Stringify and toString format supported scalar and vector types")
{
    REQUIRE(toString(true) == "true");
    REQUIRE(toString(false) == "false");

    const std::string assembled = Stringify() << "value=" << 7 << ", enabled=" << true;
    REQUIRE(assembled == "value=7, enabled=true");

    const std::string vec3f = Stringify() << osg::Vec3f(1.0f, 2.0f, 3.0f);
    const std::string vec3d = Stringify() << osg::Vec3d(4.0, 5.0, 6.0);
    const std::string vec4f = Stringify() << osg::Vec4f(0.1f, 0.2f, 0.3f, 0.4f);
    REQUIRE(vec3f == "1 2 3");
    REQUIRE(vec3d == "4 5 6");
    REQUIRE(vec4f == "0.1 0.2 0.3 0.4");
    REQUIRE(toString(osg::Vec3f(1.0f, 2.0f, 3.0f)) == "1 2 3\n");

    Stringify nested;
    nested << "inner";
    const std::string outer = Stringify() << "[" << nested << "]";
    REQUIRE(outer == "[inner]");
}
