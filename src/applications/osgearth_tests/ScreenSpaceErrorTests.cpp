/* osgEarth
 * Copyright 2026 Pelican Mapping
 * MIT License
 */
#include <osgEarth/catch.hpp>
#include <osgEarth/MapNode>
#include <osgEarth/PagedNode>
#include <osgEarth/TerrainEngineNode>
#include <osgEarthImGui/ScreenSpaceErrorControls>
#include <osg/Texture2D>
#include <osgUtil/CullVisitor>
#include <osgUtil/RenderStage>
#include <chrono>
#include <condition_variable>
#include <limits>
#include <thread>

using namespace osgEarth;
using namespace osgEarth::Util;

namespace
{
    struct CullContext
    {
        osg::ref_ptr<osg::Camera> camera = new osg::Camera;
        osg::ref_ptr<osgUtil::CullVisitor> visitor = new osgUtil::CullVisitor;
        CullContext()
        {
            visitor->setStateGraph(new osgUtil::StateGraph);
            auto* stage = new osgUtil::RenderStage;
            stage->setCamera(camera);
            visitor->setRenderStage(stage);
            visitor->setFrameStamp(new osg::FrameStamp);
            visitor->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);
            visitor->pushViewport(new osg::Viewport(0, 0, 800, 800));
            visitor->pushProjectionMatrix(new osg::RefMatrix(osg::Matrix::ortho(-10, 10, -10, 10, -100, 100)));
            visitor->pushModelViewMatrix(new osg::RefMatrix, osg::Transform::ABSOLUTE_RF);
        }
        ~CullContext()
        {
            visitor->popModelViewMatrix();
            visitor->popProjectionMatrix();
            visitor->popViewport();
        }
    };

    struct Probe : osg::Node
    {
        std::function<void(osg::NodeVisitor&)> inspect;
        Probe(const std::function<void(osg::NodeVisitor&)>& fn) : inspect(fn) { setCullingActive(false); }
        void traverse(osg::NodeVisitor& nv) override { inspect(nv); }
    };

    osg::ref_ptr<MapNode> makeMap(float sse)
    {
        MapNode::Options options;
        options.screenSpaceError() = sse;
        osg::ref_ptr<MapNode> map = new MapNode(options);
        map->open();
        // Exercise MapNode's real cull publication without rendering terrain or requiring GL.
        for (unsigned i = 0; i < map->getNumChildren(); ++i)
            map->getChild(i)->setNodeMask(0u);
        map->setCullingActive(false);
        return map;
    }

    float gpuSSE(osgUtil::StateGraph* graph)
    {
        for (; graph; graph = graph->_parent)
        {
            auto* ss = graph->getStateSet();
            auto* uniform = ss ? ss->getUniform("oe_sse") : nullptr;
            float value;
            if (uniform && uniform->get(value)) return value;
        }
        return -1.0f;
    }

    struct TestEngine : TerrainEngineNode
    {
        explicit TestEngine(float sse) { _optionsConcrete.screenSpaceError() = sse; }
        TerrainEngineRequirements requirements;
        UID getUID() const override { return 0; }
        const TerrainEngineRequirements& getRequirements() const override { return requirements; }
        osg::StateSet* getTerrainStateSet() override { return getOrCreateStateSet(); }
        unsigned getNumResidentTiles() const override { return 0; }
    };
}

TEST_CASE("Terrain panel discovers views and preserves their existing SSE resolver", "[sse][terrain-gui]")
{
    osg::ref_ptr<TestEngine> engine = new TestEngine(10.0f);
    CullContext a, b;
    a.camera->setName("Main");
    b.camera->setName("Main"); // duplicate names still have distinct row identities
    float base = 20.0f;
    engine->setScreenSpaceErrorResolver([&](osgUtil::CullVisitor& cv, float) {
        return cv.getCurrentCamera() == a.camera.get() ? base : base * 2.0f;
    });
    auto controls = ImGuiDetail::ScreenSpaceErrorControls::install(*engine);
    CHECK(controls->rows(a.visitor->getFrameStamp()).empty());
    CHECK(engine->getScreenSpaceError(*a.visitor) == 20.0f);
    CHECK(engine->getScreenSpaceError(*b.visitor) == 40.0f);
    auto rows = controls->rows(a.visitor->getFrameStamp());
    REQUIRE(rows.size() == 2u);
    CHECK(rows[0].id != rows[1].id);
    CHECK(rows[0].label != rows[1].label);
    CHECK(rows[0].label.find("Main") != std::string::npos);

    controls->set(rows[0].id, 100.0f);
    base = 30.0f;
    CHECK(engine->getScreenSpaceError(*a.visitor) == 100.0f);
    CHECK(engine->getScreenSpaceError(*b.visitor) == 60.0f);
    auto updated = controls->rows(a.visitor->getFrameStamp());
    CHECK(updated[0].id == rows[0].id);
    CHECK(updated[0].sse == 100.0f);
    CHECK(updated[1].sse == 60.0f);
    controls->reset(rows[0].id);
    CHECK(engine->getScreenSpaceError(*a.visitor) == 30.0f);

    base = std::numeric_limits<float>::quiet_NaN();
    CHECK(engine->getScreenSpaceError(*b.visitor) == 10.0f);
    CHECK(controls->rows(b.visitor->getFrameStamp())[1].sse == 10.0f);
}

TEST_CASE("Terrain panel excludes auxiliary cameras while allowing color texture views", "[sse][terrain-gui]")
{
    osg::ref_ptr<TestEngine> engine = new TestEngine(10.0f);
    engine->setScreenSpaceErrorResolver([](osgUtil::CullVisitor&, float configured) { return configured * 2.0f; });
    auto controls = ImGuiDetail::ScreenSpaceErrorControls::install(*engine);
    CullContext color, shadow, depth, pick;
    for (auto* context : { &color, &shadow, &depth, &pick })
        context->camera->setRenderTargetImplementation(osg::Camera::FRAME_BUFFER_OBJECT);
    color.camera->attach(osg::Camera::COLOR_BUFFER0, new osg::Texture2D);
    for (auto* context : { &shadow, &depth, &pick })
        context->camera->attach(osg::Camera::DEPTH_BUFFER, new osg::Texture2D);
    CameraUtils::setIsShadowCamera(shadow.camera);
    CameraUtils::setIsDepthCamera(depth.camera);
    CameraUtils::setIsPickCamera(pick.camera);
    for (auto* context : { &shadow, &depth, &pick })
        CHECK(engine->getScreenSpaceError(*context->visitor) == 20.0f);
    CHECK(controls->rows(color.visitor->getFrameStamp()).empty());

    CHECK(engine->getScreenSpaceError(*color.visitor) == 20.0f);
    auto rows = controls->rows(color.visitor->getFrameStamp());
    REQUIRE(rows.size() == 1u);
    controls->set(rows[0].id, 50.0f);
    CHECK(engine->getScreenSpaceError(*color.visitor) == 50.0f);
    for (auto* context : { &shadow, &depth, &pick })
        CHECK(engine->getScreenSpaceError(*context->visitor) == 20.0f);
    CHECK(controls->rows(color.visitor->getFrameStamp()).size() == 1u);
}

TEST_CASE("SSE panels retain the original callback instance", "[sse][terrain-gui][rendering-gui]")
{
    CullContext context;
    auto exercise = [&](auto& target) {
        target.setScreenSpaceErrorResolver([calls = 0u](osgUtil::CullVisitor&, float configured) mutable {
            return configured + static_cast<float>(++calls);
        });
        CHECK(target.getScreenSpaceError(*context.visitor) == 11.0f);
        auto controls = ImGuiDetail::ScreenSpaceErrorControls::install(target);
        CHECK(target.getScreenSpaceError(*context.visitor) == 12.0f);
        CHECK(target.getScreenSpaceError(*context.visitor) == 13.0f);
        auto rows = controls->rows(context.visitor->getFrameStamp());
        REQUIRE(rows.size() == 1u);
        CHECK(rows[0].sse == 13.0f);
    };
    SECTION("Terrain")
    {
        osg::ref_ptr<TestEngine> engine = new TestEngine(10.0f);
        exercise(*engine);
    }
    SECTION("Map")
    {
        auto map = makeMap(10.0f);
        exercise(*map);
    }
}

TEST_CASE("Terrain panel expires inactive views without retaining their cameras", "[sse][terrain-gui]")
{
    osg::ref_ptr<TestEngine> engine = new TestEngine(10.0f);
    auto controls = ImGuiDetail::ScreenSpaceErrorControls::install(*engine);
    CullContext active, removed;
    active.visitor->getFrameStamp()->setFrameNumber(10u);
    removed.visitor->getFrameStamp()->setFrameNumber(10u);
    engine->getScreenSpaceError(*active.visitor);
    engine->getScreenSpaceError(*removed.visitor);
    auto rows = controls->rows(active.visitor->getFrameStamp());
    REQUIRE(rows.size() == 2u);
    const auto activeID = rows[0].id;
    const auto removedID = rows[1].id;
    controls->set(removedID, 200.0f);

    // Keep the removed camera alive, as a viewer/cache might, but stop its culls.
    osg::ref_ptr<osg::FrameStamp> drawFrame = new osg::FrameStamp;
    drawFrame->setFrameNumber(12u);
    CHECK(controls->rows(drawFrame).size() == 2u);
    active.visitor->getFrameStamp()->setFrameNumber(13u);
    engine->getScreenSpaceError(*active.visitor);
    drawFrame->setFrameNumber(13u);
    rows = controls->rows(drawFrame);
    REQUIRE(rows.size() == 1u);
    CHECK(rows[0].id == activeID);

    // Returning cameras are rediscovered without a stale override or widget ID.
    removed.visitor->getFrameStamp()->setFrameNumber(13u);
    CHECK(engine->getScreenSpaceError(*removed.visitor) == 10.0f);
    rows = controls->rows(drawFrame);
    REQUIRE(rows.size() == 2u);
    CHECK(rows[1].id != removedID);
    controls->set(removedID, 300.0f); // a stale UI snapshot cannot edit the new row
    CHECK(engine->getScreenSpaceError(*removed.visitor) == 10.0f);

    osg::observer_ptr<osg::Camera> closed;
    {
        CullContext temporary;
        temporary.visitor->getFrameStamp()->setFrameNumber(13u);
        closed = temporary.camera;
        engine->getScreenSpaceError(*temporary.visitor);
        CHECK(controls->rows(drawFrame).size() == 3u);
    }
    CHECK_FALSE(closed.valid());
    CHECK(controls->rows(drawFrame).size() == 2u);

    // Draw can lag a cull by a frame; do not expire rows with a newer stamp.
    drawFrame->setFrameNumber(12u);
    CHECK(controls->rows(drawFrame).size() == 2u);
    drawFrame->setFrameNumber(16u);
    CHECK(controls->rows(drawFrame).empty());
}

TEST_CASE("Terrain panel shares overrides safely with concurrent culls", "[sse][terrain-gui]")
{
    osg::ref_ptr<TestEngine> engine = new TestEngine(10.0f);
    auto controls = ImGuiDetail::ScreenSpaceErrorControls::install(*engine);
    CullContext a, b;
    engine->getScreenSpaceError(*a.visitor);
    engine->getScreenSpaceError(*b.visitor);
    auto rows = controls->rows(a.visitor->getFrameStamp());
    REQUIRE(rows.size() == 2u);
    controls->set(rows[0].id, 100.0f);
    controls->set(rows[1].id, 200.0f);
    bool correctA = true, correctB = true;
    std::thread cullA([&] {
        for (unsigned i = 0; i < 1000u; ++i)
        {
            float value = engine->getScreenSpaceError(*a.visitor);
            correctA = correctA && (value == 100.0f || value == 150.0f);
        }
    });
    std::thread cullB([&] {
        for (unsigned i = 0; i < 1000u; ++i)
            correctB = correctB && engine->getScreenSpaceError(*b.visitor) == 200.0f;
    });
    for (unsigned i = 0; i < 1000u; ++i)
    {
        controls->set(rows[0].id, i % 2u == 0u ? 150.0f : 100.0f);
        controls->rows(a.visitor->getFrameStamp());
    }
    cullA.join();
    cullB.join();
    CHECK(correctA);
    CHECK(correctB);

    // Closing the panel leaves the current settings safe to use until the
    // engine replaces the resolver. The controls do not capture the panel.
    std::weak_ptr<ImGuiDetail::ScreenSpaceErrorControls> retained = controls;
    controls.reset();
    CHECK_FALSE(retained.expired());
    CHECK(engine->getScreenSpaceError(*b.visitor) == 200.0f);
    engine->setScreenSpaceErrorResolver({});
    CHECK(retained.expired());
    CHECK(engine->getScreenSpaceError(*b.visitor) == 10.0f);
}

TEST_CASE("Rendering and terrain panels independently control each view's cull SSE", "[sse][rendering-gui][terrain-gui]")
{
    auto map = makeMap(25.0f);
    osg::ref_ptr<TestEngine> engine = new TestEngine(10.0f);
    CullContext a, b, shadow;
    shadow.camera->setRenderTargetImplementation(osg::Camera::FRAME_BUFFER_OBJECT);
    shadow.camera->attach(osg::Camera::DEPTH_BUFFER, new osg::Texture2D);
    CameraUtils::setIsShadowCamera(shadow.camera);

    float mapMultiplier = 2.0f;
    map->setScreenSpaceErrorResolver([&](osgUtil::CullVisitor& cv, float configured) {
        return cv.getCurrentCamera() == b.camera.get() ? configured * mapMultiplier : configured;
    });
    engine->setScreenSpaceErrorResolver([&](osgUtil::CullVisitor& cv, float configured) {
        return cv.getCurrentCamera() == b.camera.get() ? configured * 3.0f : configured;
    });
    auto mapControls = ImGuiDetail::ScreenSpaceErrorControls::install(*map);
    auto terrainControls = ImGuiDetail::ScreenSpaceErrorControls::install(*engine);

    float mapSSE = -1.0f, gpuMapSSE = -1.0f, terrainSSE = -1.0f;
    map->addChild(new Probe([&](osg::NodeVisitor& nv) {
        auto& cv = static_cast<osgUtil::CullVisitor&>(nv);
        REQUIRE(nv.getUserValue("oe_sse", mapSSE));
        gpuMapSSE = gpuSSE(cv.getCurrentStateGraph());
        terrainSSE = engine->getScreenSpaceError(cv);
    }));
    auto cull = [&](CullContext& context, float expectedMap, float expectedTerrain) {
        map->traverse(*context.visitor);
        CHECK(mapSSE == expectedMap);
        CHECK(gpuMapSSE == expectedMap);
        CHECK(terrainSSE == expectedTerrain);
    };

    // Opening either panel preserves the application's per-camera callbacks.
    cull(a, 25.0f, 10.0f);
    cull(b, 50.0f, 30.0f);
    auto mapRows = mapControls->rows(a.visitor->getFrameStamp());
    auto terrainRows = terrainControls->rows(a.visitor->getFrameStamp());
    REQUIRE(mapRows.size() == 2u);
    REQUIRE(terrainRows.size() == 2u);
    mapControls->set(mapRows[0].id, 100.0f);
    terrainControls->set(terrainRows[1].id, 200.0f);
    cull(a, 100.0f, 10.0f);
    cull(b, 50.0f, 200.0f);
    cull(shadow, 25.0f, 10.0f);
    CHECK(mapControls->rows(a.visitor->getFrameStamp()).size() == 2u);
    CHECK(terrainControls->rows(a.visitor->getFrameStamp()).size() == 2u);

    mapMultiplier = 3.0f;
    cull(b, 75.0f, 200.0f); // unedited views still follow the original callback
    mapControls->reset(mapRows[0].id);
    cull(a, 25.0f, 10.0f);
    cull(b, 75.0f, 200.0f); // resetting map SSE does not reset terrain SSE
    terrainControls->reset(terrainRows[1].id);
    cull(b, 75.0f, 30.0f);
    CHECK(map->getScreenSpaceError() == 25.0f);
    CHECK(engine->getOptions().getScreenSpaceError() == 10.0f);

    // A camera retained by the viewer still disappears when its culls stop.
    a.visitor->getFrameStamp()->setFrameNumber(3u);
    cull(a, 25.0f, 10.0f);
    auto remaining = mapControls->rows(a.visitor->getFrameStamp());
    REQUIRE(remaining.size() == 1u);
    CHECK(remaining[0].id == mapRows[0].id);
    CHECK(terrainControls->rows(a.visitor->getFrameStamp()).size() == 1u);

    b.visitor->getFrameStamp()->setFrameNumber(3u);
    cull(b, 75.0f, 30.0f);
    auto rediscovered = mapControls->rows(b.visitor->getFrameStamp());
    REQUIRE(rediscovered.size() == 2u);
    CHECK(rediscovered[1].id != mapRows[1].id);
}

TEST_CASE("Effective SSE getters resolve per-camera values and configured fallbacks", "[sse]")
{
    auto map = makeMap(25.0f);
    osg::ref_ptr<TestEngine> engine = new TestEngine(10.0f);
    const MapNode& mapAPI = *map;
    const TerrainEngine& terrainAPI = *engine;
    CullContext a, b;
    a.visitor->setUserValue("oe_sse", 77.0f);

    CHECK(mapAPI.getScreenSpaceError(*a.visitor) == 25.0f);
    CHECK(terrainAPI.getScreenSpaceError(*a.visitor) == 10.0f);

    unsigned mapCalls = 0, terrainCalls = 0;
    map->setScreenSpaceErrorResolver([&](osgUtil::CullVisitor& cv, float configured) {
        ++mapCalls;
        return cv.getCurrentCamera() == a.camera.get() ? 0.0f : configured * 2.0f;
    });
    engine->setScreenSpaceErrorResolver([&](osgUtil::CullVisitor& cv, float configured) {
        ++terrainCalls;
        return cv.getCurrentCamera() == a.camera.get() ? 0.0f : configured * 2.0f;
    });
    CHECK(mapAPI.getScreenSpaceError(*a.visitor) == 0.0f);
    CHECK(mapAPI.getScreenSpaceError(*b.visitor) == 50.0f);
    CHECK(mapAPI.getScreenSpaceError(*b.visitor) == 50.0f);
    CHECK(terrainAPI.getScreenSpaceError(*a.visitor) == 0.0f);
    CHECK(terrainAPI.getScreenSpaceError(*b.visitor) == 20.0f);
    CHECK(terrainAPI.getScreenSpaceError(*b.visitor) == 20.0f);
    CHECK(mapCalls == 3u);
    CHECK(terrainCalls == 3u);
    CHECK(mapAPI.getScreenSpaceError() == 25.0f);
    CHECK(engine->getOptions().getScreenSpaceError() == 10.0f);

    for (float invalid : { -1.0f, std::numeric_limits<float>::infinity(),
        -std::numeric_limits<float>::infinity(), std::numeric_limits<float>::quiet_NaN() })
    {
        auto resolver = [invalid](osgUtil::CullVisitor&, float) { return invalid; };
        map->setScreenSpaceErrorResolver(resolver);
        engine->setScreenSpaceErrorResolver(resolver);
        CHECK(mapAPI.getScreenSpaceError(*a.visitor) == 25.0f);
        CHECK(terrainAPI.getScreenSpaceError(*a.visitor) == 10.0f);
    }

    map->setScreenSpaceErrorResolver({});
    engine->setScreenSpaceErrorResolver({});
    CHECK(mapAPI.getScreenSpaceError(*b.visitor) == 25.0f);
    CHECK(terrainAPI.getScreenSpaceError(*b.visitor) == 10.0f);

    float unchanged;
    REQUIRE(a.visitor->getUserValue("oe_sse", unchanged));
    CHECK(unchanged == 77.0f); // queries do not publish a cull's traversal state
    CHECK_FALSE(b.visitor->getUserValue("oe_sse", unchanged));
}

TEST_CASE("Map SSE is scoped to nested culls and retained for queued draws", "[sse]")
{
    auto outer = makeMap(25.0f);
    auto inner = makeMap(50.0f);
    CullContext context;
    std::vector<float> cpu, gpu;
    std::vector<osg::ref_ptr<osgUtil::StateGraph>> retained;
    auto inspect = [&](osg::NodeVisitor& nv) {
        float value = -1.0f;
        nv.getUserValue("oe_sse", value);
        cpu.push_back(value);
        auto* graph = static_cast<osgUtil::CullVisitor&>(nv).getCurrentStateGraph();
        gpu.push_back(gpuSSE(graph));
        retained.push_back(graph);
    };
    outer->addChild(new Probe(inspect));
    inner->addChild(new Probe(inspect));
    outer->addChild(inner);
    outer->addChild(new Probe(inspect));

    unsigned calls = 0;
    outer->setScreenSpaceErrorResolver([&](osgUtil::CullVisitor& cv, float configured) {
        ++calls;
        CHECK(cv.getCurrentCamera() == context.camera.get());
        CHECK(configured == 25.0f);
        return 100.0f;
    });
    inner->setScreenSpaceErrorResolver([](osgUtil::CullVisitor&, float) { return 200.0f; });
    outer->traverse(*context.visitor);
    CHECK(calls == 1u);
    CHECK((cpu == std::vector<float>{100.0f, 200.0f, 100.0f}));
    CHECK(gpu == cpu);
    float restored;
    CHECK_FALSE(context.visitor->getUserValue("oe_sse", restored));
    CHECK(outer->getScreenSpaceError() == 25.0f);

    // Clear the callback without altering configured defaults or queued draw state.
    outer->setScreenSpaceErrorResolver({});
    context.visitor->setUserValue("oe_sse", 77.0f);
    cpu.clear(); gpu.clear();
    outer->traverse(*context.visitor);
    CHECK((cpu == std::vector<float>{25.0f, 200.0f, 25.0f}));
    CHECK(gpu == cpu);
    REQUIRE(context.visitor->getUserValue("oe_sse", restored));
    CHECK(restored == 77.0f);
    CHECK(gpuSSE(retained[0]) == 100.0f);
    CHECK(gpuSSE(retained[1]) == 200.0f);
    CHECK(gpuSSE(retained[2]) == 100.0f);

    // A constant callback replaces the former value setter.
    const float overrideSSE = 30.0f;
    outer->setScreenSpaceErrorResolver([overrideSSE](osgUtil::CullVisitor&, float) { return overrideSSE; });
    cpu.clear(); gpu.clear();
    outer->traverse(*context.visitor);
    CHECK((cpu == std::vector<float>{30.0f, 200.0f, 30.0f}));
    CHECK(gpu == cpu);
    CHECK(outer->getScreenSpaceError() == 25.0f);
    CHECK(gpuSSE(retained[0]) == 100.0f);

    outer->setScreenSpaceErrorResolver([](osgUtil::CullVisitor&, float) {
        return std::numeric_limits<float>::quiet_NaN();
    });
    cpu.clear(); gpu.clear();
    outer->traverse(*context.visitor);
    CHECK(cpu.front() == 25.0f);
    CHECK(gpu.front() == 25.0f);
}

TEST_CASE("Independent SSE resolvers support simultaneous camera culls", "[sse]")
{
    auto map = makeMap(25.0f);
    osg::ref_ptr<TestEngine> engine = new TestEngine(10.0f);
    CullContext a, b;
    std::mutex mutex;
    std::condition_variable condition;
    unsigned entered = 0, terrainEntered = 0;
    std::atomic_uint mapCalls{0}, terrainCalls{0}, failures{0};
    map->setScreenSpaceErrorResolver([&](osgUtil::CullVisitor& cv, float configured) {
        ++mapCalls;
        // Both callbacks must enter together; a shared firing guard would fail this.
        std::unique_lock<std::mutex> lock(mutex);
        ++entered;
        condition.notify_all();
        if (!condition.wait_for(lock, std::chrono::seconds(5), [&]() { return entered == 2; })) ++failures;
        return cv.getCurrentCamera() == a.camera.get() ? configured : 100.0f;
    });
    engine->setScreenSpaceErrorResolver([&](osgUtil::CullVisitor& cv, float configured) {
        ++terrainCalls;
        std::unique_lock<std::mutex> lock(mutex);
        ++terrainEntered;
        condition.notify_all();
        if (!condition.wait_for(lock, std::chrono::seconds(5), [&]() { return terrainEntered == 2; })) ++failures;
        return cv.getCurrentCamera() == a.camera.get() ? 512.0f : configured;
    });
    map->addChild(new Probe([&](osg::NodeVisitor& nv) {
        auto& cv = static_cast<osgUtil::CullVisitor&>(nv);
        float mapSSE = -1.0f;
        nv.getUserValue("oe_sse", mapSSE);
        float terrainSSE = engine->getScreenSpaceError(cv);
        bool first = cv.getCurrentCamera() == a.camera.get();
        if (mapSSE != (first ? 25.0f : 100.0f) || terrainSSE != (first ? 512.0f : 10.0f) ||
            gpuSSE(cv.getCurrentStateGraph()) != mapSSE) ++failures;
    }));
    // Warm shared bounds on the setup thread, as a viewer update traversal does.
    map->getBound();
    std::thread t1([&]() { map->traverse(*a.visitor); });
    std::thread t2([&]() { map->traverse(*b.visitor); });
    t1.join(); t2.join();
    CHECK(failures.load() == 0u);
    CHECK(mapCalls.load() == 2u);
    CHECK(terrainCalls.load() == 2u);
    CHECK(map->getScreenSpaceError() == 25.0f);
    CHECK(engine->getOptions().getScreenSpaceError() == 10.0f);
    engine->setScreenSpaceErrorResolver({});
    CHECK(engine->getScreenSpaceError(*a.visitor) == 10.0f);
    engine->setScreenSpaceErrorResolver([](osgUtil::CullVisitor&, float) { return -1.0f; });
    CHECK(engine->getScreenSpaceError(*a.visitor) == 10.0f);
}

TEST_CASE("SSE resolvers can be replaced and cleared during an in-flight cull", "[sse]")
{
    auto exercise = [](auto install, auto evaluate, float configured) {
        std::mutex mutex;
        std::condition_variable condition;
        bool entered = false, release = false;
        auto captured = std::make_shared<float>(100.0f);
        std::weak_ptr<float> lifetime = captured;
        install([&, captured](osgUtil::CullVisitor&, float) {
            std::unique_lock<std::mutex> lock(mutex);
            entered = true;
            condition.notify_all();
            condition.wait_for(lock, std::chrono::seconds(5), [&]() { return release; });
            return *captured;
        });
        captured.reset();

        float inFlight = -1.0f;
        std::thread worker([&]() { inFlight = evaluate(); });
        bool started;
        {
            std::unique_lock<std::mutex> lock(mutex);
            started = condition.wait_for(lock, std::chrono::seconds(5), [&]() { return entered; });
        }
        if (started)
        {
            install([](osgUtil::CullVisitor&, float) { return 200.0f; });
            CHECK(evaluate() == 200.0f);
            install(MapNode::ScreenSpaceErrorResolver{});
            CHECK(evaluate() == configured);
            CHECK_FALSE(lifetime.expired());
        }
        {
            std::lock_guard<std::mutex> lock(mutex);
            release = true;
        }
        condition.notify_all();
        worker.join();
        CHECK(started);
        CHECK(inFlight == 100.0f);
        CHECK(lifetime.expired());
    };

    SECTION("Map resolver")
    {
        auto map = makeMap(25.0f);
        map->addChild(new Probe([](osg::NodeVisitor& nv) {
            float value = -1.0f;
            nv.getUserValue("oe_sse", value);
            nv.setUserValue("test.resolvedSSE", value);
        }));
        map->getBound();
        exercise([&](const MapNode::ScreenSpaceErrorResolver& resolver) {
            map->setScreenSpaceErrorResolver(resolver);
        }, [&]() {
            CullContext context;
            map->traverse(*context.visitor);
            float value = -1.0f;
            context.visitor->getUserValue("test.resolvedSSE", value);
            return value;
        }, 25.0f);
    }

    SECTION("Terrain resolver")
    {
        osg::ref_ptr<TestEngine> engine = new TestEngine(10.0f);
        exercise([&](const TerrainEngine::ScreenSpaceErrorResolver& resolver) {
            engine->setScreenSpaceErrorResolver(resolver);
        }, [&]() {
            CullContext context;
            return engine->getScreenSpaceError(*context.visitor);
        }, 10.0f);
    }
}

TEST_CASE("Paged models select using the current cull SSE", "[sse][pagednode]")
{
    CullContext context;
    auto map = makeMap(0.0f);
    float effective = 0.0f;
    map->setScreenSpaceErrorResolver([&](osgUtil::CullVisitor&, float) { return effective; });
    osg::ref_ptr<PagingManager> manager = new PagingManager("sse-test-pages");
    map->addChild(manager);
    osg::ref_ptr<PagedNode2> page = new PagedNode2;
    page->setCenter(osg::Vec3());
    page->setRadius(1.0f);
    page->setCullingActive(false);
    page->setLODMethod(LODMethod::SCREEN_SPACE);
    page->setMinPixels(1.0f);
    page->setMaxPixels(10000.0f);
    page->setAutoUnload(false);
    page->setPreCompileGLObjects(false);
    page->setRefinePolicy(REFINE_REPLACE);
    unsigned selected = 0;
    osg::ref_ptr<Probe> payload = new Probe([&](osg::NodeVisitor&) { ++selected; });
    payload->setInitialBound(osg::BoundingSphere(osg::Vec3(), 1.0f));
    page->setLoadFunction([payload](Cancelable*) { return osg::ref_ptr<osg::Node>(payload); });
    manager->addChild(page);
    osg::NodeVisitor bind(osg::NodeVisitor::TRAVERSE_ALL_CHILDREN);
    manager->accept(bind);

    map->traverse(*context.visitor);
    auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(5);
    while (!page->isLoadComplete() && std::chrono::steady_clock::now() < deadline)
    {
        manager->update();
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    REQUIRE(page->isLoadComplete());
    selected = 0;
    map->traverse(*context.visitor);
    CHECK(selected == 1u);
    effective = 10000.0f;
    map->traverse(*context.visitor);
    CHECK(selected == 1u); // loaded payload is excluded by this view's effective SSE
    effective = 0.0f;
    map->traverse(*context.visitor);
    CHECK(selected == 2u);
    CHECK(manager->sse() == 0.0f); // no update-time default mutation
}

TEST_CASE("Map and terrain overrides compose and return to configured defaults", "[sse]")
{
    auto map = makeMap(25.0f);
    osg::ref_ptr<TestEngine> engine = new TestEngine(10.0f);
    CullContext a, b;
    bool mapOverride = false, terrainOverride = false, swapped = false;
    unsigned mapCalls = 0, terrainCalls = 0;
    map->setScreenSpaceErrorResolver([&](osgUtil::CullVisitor& cv, float configured) {
        ++mapCalls;
        return mapOverride ? ((cv.getCurrentCamera() == a.camera.get()) != swapped ? 0.0f : 100.0f) : configured;
    });
    engine->setScreenSpaceErrorResolver([&](osgUtil::CullVisitor& cv, float configured) {
        ++terrainCalls;
        return terrainOverride ? ((cv.getCurrentCamera() == a.camera.get()) != swapped ? 0.0f : 512.0f) : configured;
    });
    map->addChild(new Probe([&](osg::NodeVisitor& nv) {
        auto& cv = static_cast<osgUtil::CullVisitor&>(nv);
        bool first = (cv.getCurrentCamera() == a.camera.get()) != swapped;
        float mapSSE = -1.0f;
        REQUIRE(nv.getUserValue("oe_sse", mapSSE));
        CHECK(mapSSE == (mapOverride ? (first ? 0.0f : 100.0f) : 25.0f));
        CHECK(gpuSSE(cv.getCurrentStateGraph()) == mapSSE);
        CHECK(engine->getScreenSpaceError(cv) == (terrainOverride ? (first ? 0.0f : 512.0f) : 10.0f));
    }));
    for (int mode : {0, 1, 2, 3, 0})
    {
        mapOverride = (mode & 1) != 0;
        terrainOverride = (mode & 2) != 0;
        for (bool swap : {false, true})
        {
            swapped = swap;
            map->traverse(*a.visitor);
            map->traverse(*b.visitor);
        }
    }
    CHECK(mapCalls == 20u);
    CHECK(terrainCalls == 20u);
    CHECK(map->getScreenSpaceError() == 25.0f);
    CHECK(engine->getOptions().getScreenSpaceError() == 10.0f);
}
