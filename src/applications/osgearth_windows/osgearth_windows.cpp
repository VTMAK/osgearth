/* osgEarth
* Copyright 2025 Pelican Mapping
* MIT License
*/
#include <osgEarthImGui/ImGuiApp>

#include <osgViewer/CompositeViewer>
#include <osgEarth/EarthManipulator>
#include <osgEarth/ExampleResources>
#include <osgEarth/MapNode>
#include <osgEarth/Registry>
#include <osgEarth/TerrainEngineNode>
#include <osg/PolygonMode>
#include <osg/Stats>
#include <osgUtil/CullVisitor>
#include <algorithm>
#include <mutex>
#include <iostream>

#define LC "[osgearth_windows] "

using namespace osgEarth;
using namespace osgEarth::Util;

int
usage(const char* name)
{
    OE_NOTICE
        << "\nUsage: " << name << " file.earth"
        << "\n          --views [num] : Number of windows to open"
        << "\n          --shared      : Share a single OSG graphics context across all windows (forces SingleThreaded mode)"
        << "\n          --updates [num] : Number of update traversals"
        << "\n          --sse-demo    : Two linked views, screen-space LOD, wireframe, terrain SSE 0/512, GUI"
        << "\n                          Uses concurrent culls when supported by this build; --shared uses SingleThreaded"
        << "\n          --sse-frames [num] : Exit after this many frames and print per-view cull statistics"
        << "\n"
        << MapNodeHelper().usage() << std::endl;

    return 0;
}

struct SSEControls
{
    struct View
    {
        osg::observer_ptr<osg::Camera> camera;
        osg::ref_ptr<osg::Stats> stats;
        std::string name;
        bool mapOverride = false, terrainOverride = false;
        float mapSSE = 25.0f, terrainSSE = 0.0f;
        float resolvedMap = -1.0f, resolvedTerrain = -1.0f;
    };
    std::mutex mutex;
    std::vector<View> views;
    bool linked = false, wireframe = false, screenSpace = false;
    bool dirty = false;

    float resolve(osgUtil::CullVisitor& cv, float configured, bool terrain)
    {
        std::lock_guard<std::mutex> lock(mutex);
        for (auto& view : views)
        {
            if (view.camera.get() == cv.getCurrentCamera())
            {
                float value = terrain ? (view.terrainOverride ? view.terrainSSE : configured) :
                    (view.mapOverride ? view.mapSSE : configured);
                (terrain ? view.resolvedTerrain : view.resolvedMap) = value;
                return value;
            }
        }
        return configured; // auxiliary cameras retain the configured defaults
    }
};

struct App
{
    osgViewer::CompositeViewer _viewer;
    bool _sharedGC;
    int _size;
    osg::ref_ptr<osg::Node> _node;
    osg::ref_ptr<MapNode> _mapNode;
    std::shared_ptr<SSEControls> _sse = std::make_shared<SSEControls>();
    std::mutex _commandMutex;
    std::vector<std::function<void()>> _commands;
    bool _linked = false;

    App(osg::ArgumentParser& args) :
        _viewer(args),
        _size(800)
    {
        _sharedGC = args.read("--shared");
    }

    ~App() { _viewer.stopThreading(); }

    // GUI runs on a draw thread in the concurrent demo. Scene/view changes are
    // applied by the main loop with rendering threads stopped.
    void enqueue(const std::function<void()>& command)
    {
        std::lock_guard<std::mutex> lock(_commandMutex);
        _commands.push_back(command);
    }

    void installResolvers()
    {
        if (!_mapNode || !_mapNode->open() || !_mapNode->getTerrainEngine()) return;
        auto controls = _sse;
        _mapNode->setScreenSpaceErrorResolver([controls](osgUtil::CullVisitor& cv, float configured) {
            return controls->resolve(cv, configured, false);
        });
        _mapNode->getTerrainEngine()->setScreenSpaceErrorResolver([controls](osgUtil::CullVisitor& cv, float configured) {
            return controls->resolve(cv, configured, true);
        });
        _sse->screenSpace = _mapNode->getTerrainOptions().getLODMethod() == LODMethod::SCREEN_SPACE;
    }

    void applyControls()
    {
        std::vector<std::function<void()>> commands;
        {
            std::lock_guard<std::mutex> lock(_commandMutex);
            commands.swap(_commands);
        }
        bool dirty, wireframe, screenSpace;
        {
            std::lock_guard<std::mutex> lock(_sse->mutex);
            _linked = _sse->linked;
            dirty = _sse->dirty;
            _sse->dirty = false;
            wireframe = _sse->wireframe;
            screenSpace = _sse->screenSpace;
        }
        if (commands.empty() && !dirty) return;

        bool running = _viewer.areThreadsRunning();
        _viewer.stopThreading();
        for (auto& command : commands) command();
        if (_mapNode && dirty)
        {
            auto options = _mapNode->getTerrainOptions();
            auto method = screenSpace ? LODMethod::SCREEN_SPACE : LODMethod::CAMERA_DISTANCE;
            if (options.getLODMethod() != method)
            {
                options.setLODMethod(method);
                _mapNode->getTerrainEngine()->dirtyTerrainOptions();
            }
        }
        // Replace camera state only with draw readers stopped.
        for (unsigned i = 0; i < _viewer.getNumViews(); ++i)
        {
            auto* ss = _viewer.getView(i)->getCamera()->getOrCreateStateSet();
            auto* polygon = new osg::PolygonMode(osg::PolygonMode::FRONT_AND_BACK,
                wireframe ? osg::PolygonMode::LINE : osg::PolygonMode::FILL);
            ss->setAttributeAndModes(polygon, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
        }
        if (running) _viewer.startThreading();
    }

    void linkCameras()
    {
        if (!_linked || _viewer.getNumViews() < 2) return;
        auto* source = _viewer.getView(0);
        auto* camera = source->getCamera();
        int width = camera->getGraphicsContext()->getTraits()->width;
        int height = camera->getGraphicsContext()->getTraits()->height;
        for (unsigned i = 1; i < _viewer.getNumViews(); ++i)
        {
            auto* traits = _viewer.getView(i)->getCamera()->getGraphicsContext()->getTraits();
            width = std::min(width, traits->width);
            height = std::min(height, traits->height);
        }
        for (unsigned i = 0; i < _viewer.getNumViews(); ++i)
        {
            auto* target = _viewer.getView(i)->getCamera();
            auto* traits = target->getGraphicsContext()->getTraits();
            target->setViewport((traits->width-width)/2, (traits->height-height)/2, width, height);
        }
        for (unsigned i = 1; i < _viewer.getNumViews(); ++i)
        {
            auto* view = _viewer.getView(i);
            view->getCameraManipulator()->setByMatrix(source->getCameraManipulator()->getMatrix());
            view->getCamera()->setViewMatrix(camera->getViewMatrix());
            view->getCamera()->setProjectionMatrix(camera->getProjectionMatrix());
        }
    }

    void addView(const std::string& name)
    {
        int i = _viewer.getNumViews();

        int x = 10 + i*(_size + 20);
        osg::GraphicsContext* gc_to_share = _sharedGC && i > 0 ? _viewer.getView(0)->getCamera()->getGraphicsContext() : nullptr;

        osgViewer::View* view = createView(name, x, 10, _size, _size, gc_to_share);
        if (!view) return;

        view->setCameraManipulator(new EarthManipulator());
        view->setSceneData(_node.get());
        MapNodeHelper().configureView(view);

        _viewer.addView(view);
        SSEControls::View settings;
        settings.name = name;
        settings.camera = view->getCamera();
        settings.stats = settings.camera->getStats();
        if (!settings.stats)
        {
            settings.stats = new osg::Stats(name);
            settings.camera->setStats(settings.stats);
        }
        settings.stats->collectStats("osgEarth.SSE", true);
        if (_mapNode)
        {
            settings.mapSSE = _mapNode->getScreenSpaceError();
            settings.terrainSSE = _mapNode->getTerrainOptions().getScreenSpaceError();
        }
        std::lock_guard<std::mutex> lock(_sse->mutex);
        _sse->views.push_back(settings);
    }

    osgViewer::View* createView(const std::string& name, int x, int y, int width, int height, osg::GraphicsContext* sharedGC)
    {
        osg::ref_ptr<osg::DisplaySettings>& ds = osg::DisplaySettings::instance();
        osg::ref_ptr<osg::GraphicsContext::Traits> traits = new osg::GraphicsContext::Traits(ds.get());
        traits->readDISPLAY();
        if (traits->displayNum < 0) traits->displayNum = 0;
        traits->screenNum = 0;
        traits->x = x;
        traits->y = y;
        traits->width = width;
        traits->height = height;
        traits->windowDecoration = true;
        traits->windowName = name;
        traits->doubleBuffer = true;
        traits->sharedContext = sharedGC;

        osg::GraphicsContext* gc = osg::GraphicsContext::createGraphicsContext(traits.get());
        if (!gc)
        {
            OE_WARN << "Cannot create graphics context for " << name << std::endl;
            return nullptr;
        }
        gc->setName(name);

        osgViewer::View* view = new osgViewer::View();
        view->getCamera()->setGraphicsContext(gc);

        view->getCamera()->setViewport(0, 0, width, height);
        view->getCamera()->setProjectionMatrixAsPerspective(45, 1, 1, 10);

        GLenum buffer = traits->doubleBuffer ? GL_BACK : GL_FRONT;
        view->getCamera()->setDrawBuffer(buffer);
        view->getCamera()->setReadBuffer(buffer);
        
        return view;
    }

    void releaseGLObjects(osg::State* state)
    {
        osgViewer::ViewerBase::Cameras cameras;
        _viewer.getCameras(cameras);
        for (auto c : cameras)
        {
            c->releaseGLObjects(state);
        }
    }
};

struct GCPanel : public ImGuiPanel
{
    App& _app;
    GCPanel(App& app) : ImGuiPanel("GCs"), _app(app) { }

    void draw(osg::RenderInfo& ri) override
    {
        if (!isVisible()) return;
        ImGui::Begin(name(), visible());
        auto gcs = osg::GraphicsContext::getAllRegisteredGraphicsContexts();
        for (auto gc : gcs)
        {
            if (gc->getState() != nullptr)
            {
                ImGui::Text("OSG GC ID = %d", gc->getState()->getContextID());
                ImGui::Text("OSG GC Addr: %0x", (std::uintptr_t)gc);
                ImGui::Text("Name = %s", gc->getName().c_str());
                ImGui::Text("Operations = %d", gc->getGraphicsThread() && gc->getGraphicsThread()->getOperationQueue() ? gc->getGraphicsThread()->getOperationQueue()->getNumOperationsInQueue() : 0);
                ImGui::Text("Size = %d x %d", gc->getTraits() ? gc->getTraits()->width : -1, gc->getTraits() ? gc->getTraits()->height : -1);
                if (ImGui::Button("release GL objects"))
                {
                    osg::ref_ptr<osg::State> state = gc->getState();
                    _app.enqueue([this, state]() { _app.releaseGLObjects(state); });
                }
                ImGui::Unindent();
            }
        }

        if (ImGui::Button("New window"))
        {
            std::string name = Stringify() << "Window " << _app._viewer.getNumViews();
            _app.enqueue([this, name]() { _app.addView(name); });
        }

        ImGui::End();
    }
};

struct ViewerPanel : public ImGuiPanel
{
    App& _app;
    ViewerPanel(App& app) : ImGuiPanel("Views"), _app(app) { }

    void draw(osg::RenderInfo& ri) override
    {
        if (!isVisible()) return;
        ImGui::Begin(name(), visible());

        if (ImGui::Button("New view"))
        {
            std::string name = Stringify() << "View " << _app._viewer.getNumViews();
            _app.enqueue([this, name]() { _app.addView(name); });
        }

        if (ImGui::Button("Release GL Objects"))
        {
            _app.enqueue([this]() { _app.releaseGLObjects(nullptr); });
        }

        ImGui::Text("Active Views");
        osgViewer::ViewerBase::Views views;
        _app._viewer.getViews(views);
        int ptr = 0;
        for (auto view : views)
        {
            ImGui::PushID(view);
            ImGui::Text("View #%d", ptr++);
            ImGui::Indent();

            if (view->getCamera() && view->getCamera()->getGraphicsContext() && view->getCamera()->getGraphicsContext()->getState())
            {
                auto state = view->getCamera()->getGraphicsContext()->getState();

                ImGui::Text("OSG GC ID = %d", GLUtils::getSharedContextID(*state));

                if (_app._sharedGC)
                {
                    ImGui::Text("Unique State ID = %d", GLUtils::getUniqueStateID(*state));
                }
            }
            else
            {
                ImGui::Text("**Invalid**");
            }

            if (ImGui::Button("close"))
            {
                OE_WARN << "Closing a view" << std::endl;
                osg::ref_ptr<osgViewer::View> keep = view;
                _app.enqueue([this, keep]() {
                    _app._viewer.removeView(keep);
                    keep->setSceneData(nullptr);
                    std::lock_guard<std::mutex> lock(_app._sse->mutex);
                    auto& views = _app._sse->views;
                    views.erase(std::remove_if(views.begin(), views.end(), [&keep](const SSEControls::View& item) {
                        return item.camera.get() == keep->getCamera();
                    }), views.end());
                });
            }

            ImGui::Unindent();
            ImGui::PopID();
        }            

        ImGui::End();
    }
};

struct SSEPanel : public ImGuiPanel
{
    App& _app;
    SSEPanel(App& app) : ImGuiPanel("Per-view SSE"), _app(app) { }

    void draw(osg::RenderInfo& ri) override
    {
        if (!isVisible()) return;
        ImGui::Begin(name(), visible());
        auto& controls = *_app._sse;
        std::lock_guard<std::mutex> lock(controls.mutex);
        ImGui::Checkbox("Link cameras (navigate first view)", &controls.linked);
        controls.dirty |= ImGui::Checkbox("Screen-space terrain LOD", &controls.screenSpace);
        controls.dirty |= ImGui::Checkbox("Wireframe", &controls.wireframe);
        if (ImGui::Button("Compare terrain 0 / 512") && controls.views.size() >= 2)
        {
            controls.linked = controls.screenSpace = controls.wireframe = controls.dirty = true;
            for (unsigned i = 0; i < controls.views.size(); ++i)
            {
                controls.views[i].mapOverride = false;
                controls.views[i].terrainOverride = true;
                controls.views[i].terrainSSE = i == 0 ? 0.0f : 512.0f;
            }
        }
        if (ImGui::Button("Use configured defaults"))
            for (auto& view : controls.views) view.mapOverride = view.terrainOverride = false;
        if (ImGui::Button("Swap first two overrides") && controls.views.size() >= 2)
        {
            auto& a = controls.views[0]; auto& b = controls.views[1];
            std::swap(a.mapOverride, b.mapOverride); std::swap(a.mapSSE, b.mapSSE);
            std::swap(a.terrainOverride, b.terrainOverride); std::swap(a.terrainSSE, b.terrainSSE);
        }
        ImGui::TextUnformatted("Terrain threshold = tile pixels + terrain SSE + map SSE");
        ImGui::TextUnformatted("Resources are shared; tile counts are per-view cull submissions.");
        auto threading = _app._viewer.getThreadingModel();
        ImGui::Text("Cull threading: %s", threading == osgViewer::ViewerBase::CullDrawThreadPerContext ?
            "CullDrawThreadPerContext" : threading == osgViewer::ViewerBase::SingleThreaded ?
            "SingleThreaded" : "viewer-selected");
#if defined(OSGEARTH_SINGLE_GL_CONTEXT) || defined(OSGEARTH_SINGLE_THREADED_OSG)
        ImGui::TextUnformatted("This build requires single-threaded rendering.");
#endif
        unsigned frame = ri.getState()->getFrameStamp()->getFrameNumber();
        for (auto& view : controls.views)
        {
            ImGui::PushID(view.camera.get());
            ImGui::Separator();
            ImGui::TextUnformatted(view.name.c_str());
            ImGui::Checkbox("Override Map SSE", &view.mapOverride);
            ImGui::SliderFloat("Map SSE", &view.mapSSE, 0.0f, 2048.0f);
            ImGui::Checkbox("Override Terrain SSE", &view.terrainOverride);
            ImGui::SliderFloat("Terrain SSE", &view.terrainSSE, 0.0f, 2048.0f);
            ImGui::Text("Last resolved cull: map %.1f, terrain %.1f", view.resolvedMap, view.resolvedTerrain);
            double tiles;
            if (frame > 0 && view.stats->getAttribute(frame-1, "osgEarth.SurfaceTiles", tiles))
                ImGui::Text("Surface tiles submitted (previous frame): %.0f", tiles);
            ImGui::PopID();
        }
        ImGui::End();
    }
};

int
main(int argc, char** argv)
{
    osgEarth::initialize();

    osgEarth::Registry::instance()->unRefImageDataAfterApply() = false;

    osg::ArgumentParser arguments(&argc,argv);

    // help?
    if ( arguments.read("--help") )
        return usage(argv[0]);

    int numViews = 1;
    arguments.read("--views", numViews);
    bool sseDemo = arguments.read("--sse-demo");
    bool gui = arguments.read("--gui") || sseDemo;
    if (sseDemo) numViews = std::max(2, numViews);
    numViews = std::max(1, numViews);
    unsigned sseFrames = 0;
    arguments.read("--sse-frames", sseFrames);

    int numUpdates = 1;
    arguments.read("--updates", numUpdates);


    // create a viewer:
    App app(arguments);

#ifdef OSGEARTH_SINGLE_GL_CONTEXT
    // This build maps all GL resource caches to context zero. Independent contexts
    // would reuse incompatible GL objects, and concurrent draws would race on the
    // same VirtualProgram scratch storage. Share resources and serialize rendering.
    if (!app._sharedGC)
    {
        app._sharedGC = true;
        OE_NOTICE << "Enabling shared graphics contexts because this build assumes a single GL context"
            << std::endl;
    }
#endif

#ifdef OSGEARTH_SINGLE_THREADED_OSG
    app._viewer.setThreadingModel(osgViewer::ViewerBase::SingleThreaded);
    OE_NOTICE << "Using SingleThreaded rendering because this build assumes single-threaded OSG"
        << std::endl;
#else
    // Shared OSG context IDs also share mutable rendering caches.
    if (app._sharedGC)
    {
        app._viewer.setThreadingModel(osgViewer::ViewerBase::SingleThreaded);
        OE_NOTICE << "Using SingleThreaded rendering because graphics contexts share resources" << std::endl;
    }
    else if (sseDemo)
        app._viewer.setThreadingModel(osgViewer::ViewerBase::CullDrawThreadPerContext);
#endif

    app._node = MapNodeHelper().load(arguments, &app._viewer);
    if (!app._node.get())
        return usage(argv[0]);

    app._mapNode = MapNode::get(app._node.get());
    app.installResolvers();
    if (sseDemo && (!app._mapNode || !app._mapNode->getTerrainEngine()))
    {
        OE_WARN << "--sse-demo requires a map with a terrain engine" << std::endl;
        return 1;
    }

    for(int i=0; i<numViews; ++i)
    {
        app.addView(Stringify() << "View " << i);
    }
    if (app._viewer.getNumViews() != static_cast<unsigned>(numViews))
        return 1;

    auto view = app._viewer.getView(0);

    if (sseDemo)
    {
        app._sse->linked = app._sse->screenSpace = app._sse->wireframe = app._sse->dirty = true;
        for (unsigned i = 0; i < app._sse->views.size(); ++i)
        {
            app._sse->views[i].terrainOverride = true;
            app._sse->views[i].terrainSSE = i == 0 ? 0.0f : 512.0f;
        }
        app.applyControls();
    }

    if (gui)
    {
        // install the Gui.
        auto* engine = new ImGuiAppEngine(arguments);
        engine->add(new ViewerPanel(app), !sseDemo);
        engine->add(new GCPanel(app), !sseDemo);
        if (app._mapNode) engine->add(new SSEPanel(app), true);
        view->getEventHandlers().push_front(engine);
    }

    OE_NOTICE << "Press 'n' to create a new view" << std::endl;
    EventRouter::get(view).onKeyPress(EventRouter::KEY_N, [&](...) { 
        std::cout << "Creating new view" << std::endl;
        app.enqueue([&app]() { app.addView(Stringify()<<"View " << app._viewer.getNumViews()); }); });

    OE_NOTICE << "Press 'r' to call releaseGLObjects" << std::endl;
    EventRouter::get(view).onKeyPress(EventRouter::KEY_R, [&](...) { 
        app.enqueue([&app]() { app.releaseGLObjects(nullptr); });
    });

    app._viewer.realize();
    if (sseDemo)
        OE_NOTICE << "SSE demo threading model: " <<
            (app._viewer.getThreadingModel() == osgViewer::ViewerBase::CullDrawThreadPerContext ?
                "CullDrawThreadPerContext" : "SingleThreaded/shared or viewer-selected fallback") << std::endl;

    while (!app._viewer.done())
    {
        app.applyControls();

        app._viewer.advance();
        app._viewer.eventTraversal();
        for (int i = 0; i < numUpdates; ++i)
            app._viewer.updateTraversal();
        app.linkCameras();
        app._viewer.renderingTraversals();
        if (sseFrames > 0 && app._viewer.getFrameStamp()->getFrameNumber() >= sseFrames)
            app._viewer.setDone(true);
    }
    app._viewer.stopThreading();
    if (sseFrames > 0)
    {
        unsigned frame = app._viewer.getFrameStamp()->getFrameNumber() - 1;
        for (const auto& item : app._sse->views)
        {
            double tiles = -1, mapSSE = -1, terrainSSE = -1;
            item.stats->getAttribute(frame, "osgEarth.SurfaceTiles", tiles);
            item.stats->getAttribute(frame, "osgEarth.MapSSE", mapSSE);
            item.stats->getAttribute(frame, "osgEarth.TerrainSSE", terrainSSE);
            OE_NOTICE << item.name << ": map SSE=" << mapSSE << ", terrain SSE=" << terrainSSE
                << ", surface tiles=" << tiles << std::endl;
        }
    }
    return 0;
}
