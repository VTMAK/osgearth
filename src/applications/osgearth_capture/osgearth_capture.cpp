/* osgEarth
 * Copyright 2026 Pelican Mapping
 * MIT License
 */
#include <osgEarth/MapNode>
#include <osgEarth/EarthManipulator>
#include <osgEarth/ExampleResources>
#include <osgEarth/ElevationPool>
#include <osgEarth/GLUtils>
#include <osgEarth/Registry>
#include <osgEarth/AutoClipPlaneHandler>
#include <osgViewer/Viewer>
#include <osgDB/WriteFile>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <thread>

using namespace osgEarth;

namespace
{
    //! Prints command-line usage for both modes.
    int usage(int code)
    {
        std::cout <<
            "osgearth_capture file.earth --vp lon lat height heading pitch range --out image.png\n"
            "    [--size 1600 1000] [--frames 900] [--utc-hours 12] [--prestige-sky] [--nvgl] [--samples 4]\n"
            "  Renders the map offscreen from one viewpoint and saves the frame; nothing is drawn to the desktop.\n"
            "  --frames sets how many frames (10 ms apart) terrain and model paging get to settle first.\n"
            "\n"
            "osgearth_capture file.earth --sample input.csv output.csv [--resolution 1]\n"
            "  Reads id,longitude,latitude rows (WGS84) and writes id,longitude,latitude,base_height,final_height\n"
            "  in WGS84 ellipsoid meters: base from the first open elevation layer, final from the complete\n"
            "  composed stack (grading layers included). Missing data is written as nan.\n";
        return code;
    }

    //! Samples one elevation working set at map coordinates and converts the results to WGS84 ellipsoid heights.
    //! An empty layer list selects every open elevation layer. Points without data come back as NO_DATA_VALUE.
    //! The query is split into contiguous chunks (queries are written in spatial order, so a chunk stays local)
    //! sampled concurrently, each with its own working set; the elevation pool and its tile cache are shared.
    std::vector<double> ellipsoidHeights(const Map* map, const std::vector<osg::Vec3d>& mapPoints,
        const ElevationLayerVector& layers, double resolution)
    {
        std::vector<osg::Vec3d> points = mapPoints;
        const std::size_t threads = std::max(1u, std::thread::hardware_concurrency());
        const std::size_t chunk = std::max<std::size_t>(512u, (points.size() + threads*4 - 1) / (threads*4));
        std::atomic<std::size_t> next{ 0 };
        auto worker = [&]()
        {
            ElevationPool::WorkingSet workingSet;
            if (!layers.empty()) workingSet.setElevationLayers(layers);
            for (std::size_t start = next.fetch_add(chunk); start < points.size(); start = next.fetch_add(chunk))
            {
                const auto end = points.begin() + std::min(points.size(), start + chunk);
                map->getElevationPool()->sampleMapCoords(points.begin() + start, end,
                    Distance(resolution, Units::METERS), &workingSet, nullptr);
            }
        };
        std::vector<std::thread> pool;
        for (std::size_t i = 1; i < threads; ++i) pool.emplace_back(worker);
        worker();
        for (auto& thread : pool) thread.join();
        std::vector<bool> valid(points.size());
        for (std::size_t i = 0; i < points.size(); ++i)
        {
            valid[i] = points[i].z() != NO_DATA_VALUE;
            if (!valid[i]) points[i].z() = 0.0;
        }
        // Map heights follow the map SRS vertical datum; report ellipsoid heights.
        map->getSRS()->transform(points, SpatialReference::get("wgs84"));
        std::vector<double> heights(points.size());
        for (std::size_t i = 0; i < points.size(); ++i) heights[i] = valid[i] ? points[i].z() : NO_DATA_VALUE;
        return heights;
    }

    //! Writes one height column; NO_DATA becomes nan so numeric readers keep the row count.
    void writeHeight(std::ostream& out, double height)
    {
        if (height == NO_DATA_VALUE) out << "nan";
        else out << height;
    }

    //! Sample mode: reads query rows, samples the base layer and the composed stack, and writes the result CSV.
    //! Returns nonzero when the earth file or input cannot be read or the map SRS transform fails.
    int sample(osg::ArgumentParser& args, const std::string& input, const std::string& output, double resolution)
    {
        osgEarth::initialize(args);
        osg::ref_ptr<MapNode> mapNode = MapNode::load(args);
        if (!mapNode || !mapNode->open())
        {
            std::cerr << "Cannot load the earth file\n";
            return 1;
        }
        std::ifstream in(input);
        if (!in)
        {
            std::cerr << "Cannot read " << input << '\n';
            return 1;
        }
        std::string line;
        std::getline(in, line);
        std::vector<std::string> ids;
        std::vector<osg::Vec3d> points;
        while (std::getline(in, line))
        {
            std::istringstream row(line);
            std::string id, lon, lat;
            if (!std::getline(row, id, ',') || !std::getline(row, lon, ',') || !std::getline(row, lat, ',')) continue;
            ids.push_back(id);
            points.emplace_back(std::stod(lon), std::stod(lat), 0.0);
        }
        const auto* map = mapNode->getMap();
        std::vector<osg::Vec3d> mapPoints = points;
        if (!SpatialReference::get("wgs84")->transform(mapPoints, map->getSRS()))
        {
            std::cerr << "Cannot transform the query points into the map SRS\n";
            return 1;
        }
        ElevationLayerVector layers, base;
        map->getLayers(layers);
        for (const auto& layer : layers)
            if (layer->isOpen()) { base.push_back(layer); break; }
        const auto finals = ellipsoidHeights(map, mapPoints, {}, resolution);
        const auto bases = base.empty() ? finals : ellipsoidHeights(map, mapPoints, base, resolution);
        std::ofstream out(output);
        out << "id,longitude,latitude,base_height,final_height\n" << std::setprecision(15);
        for (std::size_t i = 0; i < ids.size(); ++i)
        {
            out << ids[i] << ',' << points[i].x() << ',' << points[i].y() << ',';
            writeHeight(out, bases[i]);
            out << ',';
            writeHeight(out, finals[i]);
            out << '\n';
        }
        std::cout << "Sampled " << ids.size() << " points\n";
        return out ? 0 : 1;
    }
}

//! Offscreen capture and composed-elevation sampling for osgEarth scenes. Capture mode renders into a pbuffer,
//! lets paging settle for a fixed number of frames, then saves the frame; sample mode answers height queries
//! from the same composed elevation the renderer uses (see usage()).
int main(int argc, char** argv)
{
    osg::ArgumentParser args(&argc, argv);
    if (args.read("--help") || args.read("-h")) return usage(0);
    std::string sampleIn, sampleOut;
    double resolution = 1.0;
    args.read("--resolution", resolution);
    if (args.read("--sample", sampleIn, sampleOut)) return sample(args, sampleIn, sampleOut, resolution);

    std::string out = "capture.png";
    unsigned width = 1600u, height = 1000u, frames = 900u;
    double lon = 0, lat = 0, alt = 0, heading = 0, pitch = -30, range = 500, hours = -1.0;
    args.read("--out", out);
    args.read("--size", width, height);
    args.read("--frames", frames);
    const bool viewpoint = args.read("--vp", lon, lat, alt, heading, pitch, range);
    args.read("--utc-hours", hours);
    if (!viewpoint || width < 16u || height < 16u) return usage(1);

    osgEarth::initialize(args);
    osg::DisplaySettings::instance()->readCommandLine(args);
    if (args.read("--nvgl")) GLUtils::useNVGL(true);
    osg::ref_ptr<MapNode> mapNode = MapNode::load(args);
    if (!mapNode || !mapNode->open())
    {
        std::cerr << "Cannot load the earth file\n";
        return 1;
    }
    osgViewer::Viewer viewer;
    viewer.setThreadingModel(osgViewer::Viewer::SingleThreaded);
    viewer.setRealizeOperation(new GL3RealizeOperation());
    auto manipulator = new Util::EarthManipulator();
    // A fixed capture view must not be pushed outward by coarse terrain during paging.
    manipulator->getSettings()->setTerrainAvoidanceEnabled(false);
    manipulator->getSettings()->setMinMaxDistance(0.5, 1e8);
    viewer.setCameraManipulator(manipulator);

    auto traits = new osg::GraphicsContext::Traits(osg::DisplaySettings::instance());
    traits->readDISPLAY();
    traits->setUndefinedScreenDetailsToDefaultScreen();
    traits->width = width;
    traits->height = height;
    traits->pbuffer = true;
    traits->doubleBuffer = false;
    osg::ref_ptr<osg::GraphicsContext> context = osg::GraphicsContext::createGraphicsContext(traits);
    if (!context)
    {
        std::cerr << "Cannot create an offscreen graphics context\n";
        return 1;
    }
    context->getState()->setUseVertexAttributeAliasing(true);
    context->getState()->setUseModelViewAndProjectionUniforms(true);
    viewer.getCamera()->setGraphicsContext(context);
    viewer.getCamera()->setViewport(0, 0, width, height);
    viewer.getCamera()->setDrawBuffer(GL_FRONT);
    viewer.getCamera()->setReadBuffer(GL_FRONT);
    viewer.getCamera()->setProjectionMatrixAsPerspective(35.0, double(width)/height, 0.5, 1e7);
    viewer.getCamera()->setSmallFeatureCullingPixelSize(-1.0f);

    auto root = new osg::Group();
    root->addChild(mapNode);
    Util::MapNodeHelper().parse(mapNode, args, &viewer, root);
    GLUtils::setGlobalDefaults(viewer.getCamera()->getOrCreateStateSet());
    auto clip = new Util::AutoClipPlaneCullCallback(mapNode);
    clip->setMinNearFarRatio(0.00001);
    viewer.getCamera()->addCullCallback(clip);
    if (hours >= 0.0 && hours < 24.0)
    {
        const auto date = Registry::instance()->getDateTime();
        Registry::instance()->setDateTime(DateTime(date.year(), date.month(), date.day(), hours));
    }
    viewer.setSceneData(root);
    Viewpoint vp;
    vp.setFocalPoint(GeoPoint(SpatialReference::get("wgs84"), lon, lat, alt, ALTMODE_ABSOLUTE));
    vp.setHeading(Angle(heading, Units::DEGREES));
    vp.setPitch(Angle(pitch, Units::DEGREES));
    vp.setRange(Distance(range, Units::METERS));
    manipulator->setViewpoint(vp, 0.0);
    viewer.realize();
    if (!viewer.isRealized()) return 1;
    // The earth file's viewpoint extension may set its own home view on the first frame; override it after.
    viewer.frame();
    for (unsigned frame = 0; frame < frames; ++frame)
    {
        manipulator->setViewpoint(vp, 0.0);
        viewer.frame();
        // Give asynchronous terrain and model paging wall-clock time; this is a visual fixture, not a benchmark.
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    if (!context->makeCurrent()) return 1;
    osg::ref_ptr<osg::Image> image = new osg::Image();
    image->readPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE);
    const bool saved = osgDB::writeImageFile(*image, out);
    context->releaseContext();
    std::cout << (saved ? "Captured " : "Failed to save ") << out << '\n';
    return saved ? 0 : 1;
}
