/* osgEarth
* Copyright 2025 Pelican Mapping
* MIT License
*/
#include <osgEarth/MapNode>
#include <osgEarth/CascadeDrapingDecorator>
#include <osgEarth/ClampingTechnique>
#include <osgEarth/CullingUtils>
#include <osgEarth/DrapingTechnique>
#include <osgEarth/NodeUtils>
#include <osgEarth/Registry>
#include <osgEarth/Lighting>
#include <osgEarth/GLUtils>
#include <osgEarth/HorizonClipPlane>
#include <osgEarth/SceneGraphCallback>
#include <osgEarth/MapNodeObserver>
#include <osgEarth/Utils>
#include <osgEarth/Shaders>
#include <osgEarth/Chonk>
#include <osgEarth/HTTPClient>

#include <osgDB/DatabasePager>
#include <osgUtil/Optimizer>

#ifdef OSGEARTH_HAVE_SUPERLUMINALAPI
#include <Superluminal/PerformanceAPI.h>
#endif

using namespace osgEarth;
using namespace osgEarth::Util;
using namespace osgEarth::Contrib;

#define LC "[MapNode] "

namespace
{
    /**
     * A group that the osgUtil::Optimizer won't remove even when it's empty.
     */
    struct StickyGroup : public osg::Group { };

    // adapter that lets MapNode listen to Map events
    struct MapNodeMapCallbackProxy : public MapCallback
    {
        MapNodeMapCallbackProxy(MapNode* node) : _node(node) { }

        void onLayerAdded(Layer* layer, unsigned index) {
            _node->onLayerAdded(layer, index);
        }
        void onLayerRemoved(Layer* layer, unsigned index) {
            _node->onLayerRemoved(layer, index);
        }
        void onLayerMoved(Layer* layer, unsigned oldIndex, unsigned newIndex) {
            _node->onLayerMoved(layer, oldIndex, newIndex);
        }
        void onLayerOpened(Layer* layer) {
            _node->onLayerAdded(layer, _node->getMap()->getIndexOfLayer(layer));
        }
        void onLayerClosed(Layer* layer) {
            _node->onLayerRemoved(layer, _node->getMap()->getIndexOfLayer(layer));
        }

        osg::observer_ptr<MapNode> _node;
    };

    // callback that will run the MapNode installer on model layers so that
    // MapNodeObservers can have MapNode access
    struct MapNodeObserverInstaller : public SceneGraphCallback
    {
        MapNodeObserverInstaller( MapNode* mapNode ) : _mapNode( mapNode ) { }

        virtual void onPostMergeNode(osg::Node* node, osg::Object* sender)
        {
            if ( _mapNode.valid() && node )
            {
                MapNodeReplacer replacer( _mapNode.get() );
                node->accept( replacer );
            }
        }

        osg::observer_ptr<MapNode> _mapNode;
    };

    // proxy cull callback for layers that have their own cull callback
    struct LayerCullCallbackDispatch : public osg::NodeCallback
    {
        Layer* _layer;

        LayerCullCallbackDispatch(Layer* layer) : _layer(layer) { }

        virtual void operator()(osg::Node* node, osg::NodeVisitor* nv)
        {
#ifdef OSGEARTH_HAVE_SUPERLUMINALAPI
            // Track which layer is being culled
            PERFORMANCEAPI_INSTRUMENT_FUNCTION();
            PERFORMANCEAPI_INSTRUMENT_DATA("layer", _layer->getName().c_str());
#endif

            if (_layer->getNode())
            {
                _layer->apply(_layer->getNode(), nv);
            }
            else
            {
                traverse(node, nv);
            }
        };
    };

    typedef std::vector< osg::ref_ptr<Extension> > Extensions;

    class RemoveBlacklistedFilenamesVisitor : public osg::NodeVisitor
    {
    public:
        RemoveBlacklistedFilenamesVisitor():
          osg::NodeVisitor(osg::NodeVisitor::TRAVERSE_ALL_CHILDREN),
              _numRemoved(0)
          {
          }

          virtual void apply(osg::PagedLOD& node)
          {
              //The PagedLOD node will contain two filenames, the first is empty and is the actual geometry of the
              //tile and the second is the filename of the next tile.
              if (node.getNumFileNames() > 1)
              {
                  //Get the child filename
                  const std::string &filename = node.getFileName(1);
                  if (osgEarth::Registry::instance()->isBlacklisted(filename))
                  {
                      //If the tile is blacklisted, we set the actual geometry, child 0, to always display
                      //and the second child to never display
                      node.setRange(0, 0, FLT_MAX);
                      node.setRange(1, FLT_MAX, FLT_MAX);
                  }
              }
              traverse(node);
          }

          unsigned int _numRemoved;
    };
}

//---------------------------------------------------------------------------

MapNode*
MapNode::load(osg::ArgumentParser& args)
{
    for( int i=1; i<args.argc(); ++i )
    {
        if ( args[i] ) //&& (endsWith(args[i], ".earth") || endsWith(args[i], ".earth.template")) )
        {
            ReadResult r = URI(args[i]).readNode();
            if ( r.succeeded() )
            {
                return r.release<MapNode>();
            }
        }
    }
    return 0L;
}

MapNode*
MapNode::load(osg::ArgumentParser& args, const MapNode::Options& defaults)
{
    for( int i=1; i<args.argc(); ++i )
    {
        if ( args[i] ) //&& (endsWith(args[i], ".earth") || endsWith(args[i], ".earth.template")))
        {
            osg::ref_ptr<osgDB::Options> dbo = new osgDB::Options();
            std::string optionsJSON = defaults.getConfig().toJSON();
            dbo->setPluginStringData( "osgEarth.defaultOptions", optionsJSON );
            ReadResult r = URI(args[i]).readNode( dbo.get() );
            if ( r.succeeded() )
            {
                return r.release<MapNode>();
            }
        }
    }
    return 0L;
}

//....................................................................

Config
MapNode::Options::getConfig() const
{
    Config conf("options");
    conf.set( "proxy",                    proxySettings() );
    conf.set( "lighting",                 enableLighting() );
    conf.set( "overlay_blending",         overlayBlending() );
    conf.set( "overlay_blending_source",  overlayBlendingSource());
    conf.set( "overlay_texture_size",     overlayTextureSize() );
    conf.set( "overlay_mipmapping",       overlayMipMapping() );
    conf.set( "overlay_resolution_ratio", overlayResolutionRatio() );
    conf.set( "cascade_draping",          useCascadeDraping() );
    conf.set( "draping_render_bin_number",drapingRenderBinNumber() );
    conf.set("screen_space_error", screenSpaceError());

    if (terrain().isSet() && !terrain()->empty())
        conf.set( "terrain", terrain()->getConfig() );

    return conf;
}

void
MapNode::Options::fromConfig(const Config& conf)
{
    conf.get( "proxy",                    proxySettings() );
    conf.get( "lighting",                 enableLighting() );
    conf.get( "overlay_blending",         overlayBlending() );
    conf.get( "overlay_blending_source",  overlayBlendingSource());
    conf.get( "overlay_texture_size",     overlayTextureSize() );
    conf.get( "overlay_mipmapping",       overlayMipMapping() );
    conf.get( "overlay_resolution_ratio", overlayResolutionRatio() );
    conf.get( "cascade_draping",          useCascadeDraping() );
    conf.get( "draping_render_bin_number",drapingRenderBinNumber() );
    conf.get("screen_space_error", screenSpaceError());

    if ( conf.hasChild( "terrain" ) )
        terrain() = TerrainOptions( conf.child("terrain") );
}

//....................................................................

MapNode::MapNode() :
_map( new Map() )
{
    init();
}

MapNode::MapNode( Map* map ) :
_map( map ? map : new Map() )
{
    init();
}

MapNode::MapNode(const MapNode::Options& options ) :
_map( new Map() ),
_optionsConcrete(options)
{
    init();
}

MapNode::MapNode(Map* map, const MapNode::Options& options) :
_map( map? map : new Map() ),
_optionsConcrete(options)
{
    init();
}

// Common CTOR method
void
MapNode::init()
{
    // Protect the MapNode from the Optimizer
    setDataVariance(osg::Object::DYNAMIC);

    // Protect the MapNode from the ShaderGenerator
    ShaderGenerator::setIgnoreHint(this, true);

    // initialize 0Ls
    _terrainEngine = 0L;
    _terrainGroup = 0L;
    _layerNodes = 0L;
    _lastNumBlacklistedFilenames = 0;
    _isOpen = false;

    setName( "osgEarth::MapNode" );

    // Construct the container for the terrain engine, which also hold the options.
    _terrainGroup = new StickyGroup();
    this->addChild(_terrainGroup);

    // make a group for the model layers.  This node is a PagingManager instead of a regular Group to allow PagedNode's to be used within the layers.
    auto pagingManager = new PagingManager();
    pagingManager->setName("osgEarth::MapNode.layerNodes");
    pagingManager->setMaxMergesPerFrame(options().terrain()->mergesPerFrame().value());
    _layerNodes = pagingManager;

    this->addChild( _layerNodes );

    // Make sure the Registry is not destroyed until we are done using
    // it (in ~MapNode).
    //_registry = Registry::instance();

    // the default SSE for all supporting geometries
    _sseU = new osg::Uniform("oe_sse", options().screenSpaceError().get());

    _readyForUpdate = true;
}

bool
MapNode::open()
{
    if (_isOpen)
        return _isOpen;

    _isOpen = true;

    // Take a reference to this object so that it doesn't get inadvertently
    // deleting during startup. It is possible that during startup, a driver
    // will load that will take a reference to the MapNode and we don't want
    // that deleting the MapNode out from under us.
    // This is paired by an unref_nodelete() at the end of this method.
    this->ref();

    // Set the global proxy settings
    // TODO: this should probably happen elsewhere, like in the registry?
    if ( options().proxySettings().isSet() )
    {
        HTTPClient::setProxySettings( options().proxySettings().get() );
    }

    // Install a callback that lets this MapNode know about any changes
    // to the map, and invoke it manually now so they start out in sync.
    _mapCallback = new MapNodeMapCallbackProxy(this);
    _map->addMapCallback( _mapCallback.get() );
    _mapCallback->invokeOnLayerAdded(_map.get());

    // load and attach the terrain engine.
    _terrainEngine = TerrainEngineNode::create(options().terrain().get());
    if ( _terrainEngine )
    {
        _terrainEngine->setMap(_map.get(), options().terrain().get());
    }
    else
    {
        OE_WARN << "FAILED to create a terrain engine for this map" << std::endl;
    }

    // now that the terrain engine exists, we can invoke the rendering callback
    LayerVector layers;
    _map->getLayers(layers, [&](const Layer* layer) { return layer->isOpen(); });
    for (auto& layer : layers)
        layer->invoke_prepareForRendering(getTerrainEngine());

    // initialize terrain-level lighting:
    if ( options().terrain()->enableLighting().isSet() )
    {
        GLUtils::setLighting(
            _terrainGroup->getOrCreateStateSet(),
            options().terrain()->enableLighting().get() ? 1 : 0 );
    }

    if (_terrainEngine.valid())
    {
        // a decorator for draped geometry (projected texturing)
        OverlayDecorator* overlayDecorator = new OverlayDecorator();
        _terrainGroup->addChild(overlayDecorator);

        // install the Clamping technique for overlays:
        ClampingTechnique* clamping = new ClampingTechnique();
        overlayDecorator->addTechnique(clamping);
        _clampingManager = &clamping->getClampingManager();

        bool envUseCascadedDraping = (::getenv("OSGEARTH_USE_CASCADE_DRAPING") != 0L);
        if (envUseCascadedDraping || options().useCascadeDraping() == true)
        {
            CascadeDrapingDecorator* cascadeDrapingDecorator = new CascadeDrapingDecorator(getMapSRS(), _terrainEngine->getResources());
            overlayDecorator->addChild(cascadeDrapingDecorator);
            _drapingManager = cascadeDrapingDecorator->getDrapingManager();
            cascadeDrapingDecorator->addChild(_terrainEngine);
        }

        else
        {
            // simple draping - faster but less accurate

            DrapingTechnique* draping = new DrapingTechnique();

            const char* envOverlayTextureSize = ::getenv("OSGEARTH_OVERLAY_TEXTURE_SIZE");

            if (options().overlayBlending().isSet())
                draping->setOverlayBlending(options().overlayBlending().get());

            if (options().overlayBlendingSource().isSetTo("color"))
                draping->setOverlayBlendingParams(GL_SRC_COLOR, GL_ONE_MINUS_SRC_COLOR);
            else if (options().overlayBlendingSource().isSetTo("alpha"))
                draping->setOverlayBlendingParams(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

            if (envOverlayTextureSize)
                draping->setTextureSize(as<int>(envOverlayTextureSize, 1024));

            else if (options().overlayTextureSize().isSet())
                draping->setTextureSize(options().overlayTextureSize().get());

            if (options().overlayMipMapping().isSet())
                draping->setMipMapping(options().overlayMipMapping().get());

            if (options().overlayResolutionRatio().isSet())
                draping->setResolutionRatio(options().overlayResolutionRatio().get());

            draping->reestablish(_terrainEngine);
            overlayDecorator->addTechnique(draping);
            _drapingManager = draping->getDrapingManager();

            if (options().drapingRenderBinNumber().isSet())
                _drapingManager->setRenderBinNumber(options().drapingRenderBinNumber().get());

            overlayDecorator->addChild(_terrainEngine);
        }

        overlayDecorator->setTerrainEngine(_terrainEngine);
    }

    osg::StateSet* stateset = getOrCreateStateSet();

    stateset->addUniform(_sseU.get());

    if ( options().enableLighting().isSet() )
    {
        setEnableLighting(options().enableLighting().get());
    }

    // A shader define indicating that this is a geocentric display
    if ( _map->getSRS()->isGeographic() )
    {
        stateset->setDefine("OE_IS_GEOCENTRIC");
    }

    // VRV_PATCH, don't use OSG's GL material system
    // install a default material for everything in the map
    osg::Material* defaultMaterial = new osg::Material();
    defaultMaterial->setDiffuse(defaultMaterial->FRONT, osg::Vec4(1,1,1,1));
    defaultMaterial->setAmbient(defaultMaterial->FRONT, osg::Vec4(1,1,1,1));
    stateset->setAttributeAndModes(defaultMaterial, 1);

    // Note: without this, osgEarth's Sky won't work -gw
    //MaterialCallback().operator()(defaultMaterial, 0L);

    // activate PBR support.
    VirtualProgram* vp = VirtualProgram::getOrCreate(stateset);
    vp->setName(className());
    Shaders shaders;

    shaders.load(vp, shaders.PBR);
    stateset->setDefine("OE_USE_PBR");

    shaders.load(vp, shaders.HexTilingLib);
    stateset->setDefine("OE_HAVE_HEX_TILING");

    dirtyBound();

    // install a callback that sets the viewport size uniform:
    this->addCullCallback(new InstallCameraUniform());

    // install a callback that updates a horizon object and installs a clipping plane
    if (getMapSRS()->isGeographic())
    {
        this->addCullCallback(new HorizonClipPlane(getMapSRS()->getEllipsoid()));
    }

    // connect any extensions that have already been added.
    for(Extensions::iterator i = _extensions.begin(); i != _extensions.end(); ++i)
    {
        ExtensionInterface<MapNode>* extensionIF = ExtensionInterface<MapNode>::get(i->get());
        if (extensionIF)
        {
            extensionIF->connect(this);
        }
    }

    // register for event traversals so we can deal with blacklisted filenames
    ADJUST_EVENT_TRAV_COUNT( this, 1 );

    // remove the temporary reference.
    this->unref_nodelete();

    return true;
}

MapNode::~MapNode()
{
    releaseGLObjects(nullptr);

    if (_terrainEngine)
        _terrainEngine->shutdown();

    if (_mapCallback.valid())
    {
        // Remove this node's map callback first:
        _map->removeMapCallback( _mapCallback.get() );

        // Then invoke "removed" on all the layers in a batch.
        _mapCallback->invokeOnLayerRemoved(_map.get());
    }

    this->clearExtensions();

    osg::observer_ptr<TerrainEngineNode> te = _terrainEngine;
    removeChildren(0, getNumChildren());
}

void
MapNode::setEnableLighting(const bool& value)
{
    options().enableLighting() = value;
    
    GLUtils::setLighting(
        getOrCreateStateSet(),
        options().enableLighting().value() ? 1 : 0);
}

const bool&
MapNode::getEnableLighting() const
{
    return options().enableLighting().get();
}

Config
MapNode::getConfig() const
{
    Config mapConf("map");
    mapConf.set("version", "3");

    // the map and node options:
    const Map* cmap = _map.get();
    Config optionsConf = cmap->options().getConfig();
    optionsConf.merge( options().getConfig() );
    if (!optionsConf.empty())
    {
        mapConf.add("options", optionsConf);
    }

    // the layers
    LayerVector layers;
    _map->getLayers(layers);
    for(auto& layer : layers)
    {
        Config layerConf = layer->getConfig();
        if (!layerConf.empty() && !layerConf.key().empty())
        {
            mapConf.add(layerConf);
        }
    }

    typedef std::vector< osg::ref_ptr<Extension> > Extensions;
    for(Extensions::const_iterator i = getExtensions().begin(); i != getExtensions().end(); ++i)
    {
        Extension* e = i->get();
        Config conf = e->getConfigOptions().getConfig();
        conf.key() = e->getConfigKey();
        mapConf.add( conf );
    }

    Config ext = externalConfig();
    if ( !ext.empty() )
    {
        ext.key() = "external";
        mapConf.add( ext );
    }

    return mapConf;
}

osg::BoundingSphere
MapNode::computeBound() const
{
    osg::BoundingSphere bs;
    if ( getTerrainEngine() )
    {
        bs.expandBy( getTerrainEngine()->getNode()->getBound() );
    }

    if (_layerNodes)
    {
        bs.expandBy( _layerNodes->getBound() );
    }

    return bs;
}

Map*
MapNode::getMap()
{
    return _map.get();
}

const Map*
MapNode::getMap() const
{
    return _map.get();
}

const SpatialReference*
MapNode::getMapSRS() const
{
    return getMap() && getMap()->getProfile() ?
        getMap()->getProfile()->getSRS() :
        0L;
}

TerrainOptionsAPI
MapNode::getTerrainOptions()
{
    if (_terrainEngine.valid())
        return _terrainEngine->getOptions();
    else
        return TerrainOptionsAPI(&_optionsConcrete.terrain().mutable_value());
}

Terrain*
MapNode::getTerrain()
{
    if (getTerrainEngine() == NULL)
          return NULL;
    return getTerrainEngine()->getTerrain();
}

const Terrain*
MapNode::getTerrain() const
{
    if (getTerrainEngine() == NULL)
        return NULL;
    return getTerrainEngine()->getTerrain();
}

TerrainEngine*
MapNode::getTerrainEngine() const
{
    return _terrainEngine;
}

void
MapNode::setScreenSpaceError(float value)
{
    // global option:
    options().screenSpaceError() = value;

    // update the corresponding terrain option:
    getTerrainOptions().setScreenSpaceError(value);

    // update the uniform:
    _sseU->set(value);
}

float
MapNode::getScreenSpaceError() const
{
    return options().screenSpaceError().get();
}

void
MapNode::addExtension(Extension* extension, const osgDB::Options* options)
{
    if ( extension )
    {
        _extensions.push_back( extension );

        // set the IO options is they were provided:
        if ( options )
            extension->setDBOptions( options );

        else if ( getMap()->getReadOptions() )
            extension->setDBOptions( getMap()->getReadOptions() );

        // start it if the map is open; otherwise this will happen during MapNode::open.
        if (_isOpen)
        {
            ExtensionInterface<MapNode>* extensionIF = ExtensionInterface<MapNode>::get(extension);
            if ( extensionIF )
            {
                extensionIF->connect( this );
            }
        }

        //OE_INFO << LC << "Added extension \"" << extension->getName() << "\"\n";
    }
}

void
MapNode::removeExtension(Extension* extension)
{
    Extensions::iterator i = std::find(_extensions.begin(), _extensions.end(), extension);
    if ( i != _extensions.end() )
    {
        ExtensionInterface<MapNode>* extensionIF = ExtensionInterface<MapNode>::get( i->get() );
        if ( extensionIF )
        {
            extensionIF->disconnect( this );
        }
        _extensions.erase( i );
    }
}

void
MapNode::clearExtensions()
{
    for(Extensions::iterator i = _extensions.begin(); i != _extensions.end(); ++i)
    {
        ExtensionInterface<MapNode>* extensionIF = ExtensionInterface<MapNode>::get( i->get() );
        if ( extensionIF )
        {
            extensionIF->disconnect( this );
        }
    }

    _extensions.clear();
}

osg::Group*
MapNode::getLayerNodeGroup() const
{
    return _layerNodes;
}

MapNode*
MapNode::findMapNode( osg::Node* graph, unsigned travmask )
{
    return findRelativeNodeOfType<MapNode>( graph, travmask );
}

bool
MapNode::isGeocentric() const
{
    return _map->getSRS()? _map->getSRS()->isGeographic() : true;
}

namespace
{
    void rebuildLayerNodes(const Map* map, osg::Group* layerNodes)
    {
        layerNodes->removeChildren(0, layerNodes->getNumChildren());

        LayerVector layers;
        map->getLayers(layers);
        for (LayerVector::iterator i = layers.begin(); i != layers.end(); ++i)
        {
            Layer* layer = i->get();

            if (layer->isOpen())
            {
                osg::Node* node = layer->getNode();
                if (node)
                {
                    osg::Group* container = new osg::Group();
                    container->setName(layer->getName());
                    container->addChild(node);
                    container->setStateSet(layer->getOrCreateStateSet());
                    container->setCullCallback(new LayerCullCallbackDispatch(layer));
                    layerNodes->addChild(container);
                }
            }
        }
    }
}

void
MapNode::onLayerAdded(Layer* layer, unsigned index)
{
    if (!layer || !layer->isOpen())
        return;

    if (getTerrainEngine())
        layer->invoke_prepareForRendering(getTerrainEngine());

    // Create the layer's node, if it has one:
    osg::Node* node = layer->getNode();
    if (node)
    {
        // notify before adding it to the graph:
        layer->getSceneGraphCallbacks()->firePreMergeNode(node);

        // update the layer-to-node table (adds the node to the graph)
        rebuildLayerNodes(_map.get(), _layerNodes);

        // after putting it in the graph:
        layer->getSceneGraphCallbacks()->firePostMergeNode(node);
    }
}

void
MapNode::onLayerRemoved(Layer* layer, unsigned index)
{
    if (layer)
    {
        osg::Node* node = layer->getNode();
        if (node)
        {
            layer->getSceneGraphCallbacks()->fireRemoveNode(node);
        }
        rebuildLayerNodes(_map.get(), _layerNodes);
    }
}

void
MapNode::onLayerMoved(Layer* layer, unsigned oldIndex, unsigned newIndex)
{
    if (layer)
    {
        rebuildLayerNodes(_map.get(), _layerNodes);
    }
}

void
MapNode::openMapLayers()
{
    LayerVector layers;
    _map->getLayers(layers);

    for (LayerVector::const_iterator i = layers.begin(); i != layers.end(); ++i)
    {
        Layer* layer = i->get();

        if (!layer->getStatus().isError())
        {
            const Status& status = layer->open();
            if (status.isError())
            {
                OE_WARN << LC << "Failed to open layer \"" << layer->getName() << "\" ... " << status.message() << std::endl;
            }
        }
    }
}

void
MapNode::traverse( osg::NodeVisitor& nv )
{
    // ensure the map node is open during the very first "rendering" traversal
    if ( !_isOpen )
    {
        if (nv.getVisitorType() == nv.EVENT_VISITOR ||
            nv.getVisitorType() == nv.CULL_VISITOR ||
            nv.getVisitorType() == nv.UPDATE_VISITOR)
        {
            static std::mutex s_openMutex;
            std::lock_guard<std::mutex> lock(s_openMutex);
            if (!_isOpen)
            {
                _isOpen = open();
            }
        }
    }

    if ( nv.getVisitorType() == nv.EVENT_VISITOR )
    {
        unsigned int numBlacklist = Registry::instance()->getNumBlacklistedFilenames();
        if (numBlacklist != _lastNumBlacklistedFilenames)
        {
            //Only remove the blacklisted filenames if new filenames have been added since last time.
            _lastNumBlacklistedFilenames = numBlacklist;
            RemoveBlacklistedFilenamesVisitor v;
            _terrainEngine->accept( v );

        }

        // traverse:
        for (auto& child : _children)
            child->accept(nv);
    }

    else if ( nv.getVisitorType() == nv.CULL_VISITOR)
    {
        osgUtil::CullVisitor* cv = Culling::asCullVisitor(nv);

        // find a database pager with an ICO
        if (cv && cv->getCurrentCamera())
        {
            // Store this MapNode itself
            ObjectStorage::set(&nv, this);

            osgDB::DatabasePager* pager = dynamic_cast<osgDB::DatabasePager*>(
                nv.getDatabaseRequestHandler());

            if (pager)
                ObjectStorage::set(&nv, pager->getIncrementalCompileOperation());

            if (_drapingManager != nullptr)
                ObjectStorage::set(&nv, _drapingManager);
        }


        LayerVector layers;
        getMap()->getLayers(layers);

        int count = 0;
        for (auto& layer : layers)
        {
            if (layer->isOpen())
            {
                osg::StateSet* ss = layer->getSharedStateSet(&nv);
                if (ss)
                {
                    cv->pushStateSet(ss);
                    ++count;
                }
            }
        }

        // traverse:
        for (auto& child : _children)
            child->accept(nv);

        for(int i=0; i<count; ++i)
            cv->popStateSet();

        //Config c = CullDebugger().dumpRenderBin(cv->getCurrentRenderBin());
        //OE_INFO << c.toJSON(true) << std::endl;
        //exit(0);

        // after any cull, allow an update traversal.
        _readyForUpdate.exchange(true);
    }

    else if (nv.getVisitorType() == nv.UPDATE_VISITOR)
    {
        ObjectStorage::set(&nv, this);

        // Ensures only one update will happen per frame loop
        if (_readyForUpdate.exchange(false))
        {
            // re-enable if we decide to use it.
            //JobArena::get(JobArena::UPDATE_TRAVERSAL)->runJobs();
        }

        // include these in the above condition as well??
        osg::Group::traverse(nv);
    }

    else
    {
        if (dynamic_cast<osgUtil::BaseOptimizerVisitor*>(&nv) == 0L)
        {
            // Store this MapNode itself
            ObjectStorage::set(&nv, this);

            osg::Group::traverse(nv);
        }
    }
}

void
MapNode::resizeGLObjectBuffers(unsigned maxSize)
{
    LayerVector layers;
    getMap()->getLayers(layers);
    for (LayerVector::const_iterator i = layers.begin(); i != layers.end(); ++i)
    {
        i->get()->resizeGLObjectBuffers(maxSize);
    }

    osg::Group::resizeGLObjectBuffers(maxSize);
}

void
MapNode::releaseGLObjects(osg::State* state) const
{
    LayerVector layers;
    getMap()->getLayers(layers);
    for (LayerVector::const_iterator i = layers.begin(); i != layers.end(); ++i)
    {
        i->get()->releaseGLObjects(state);
    }

    // osg::Node doesn't release nested callbacks, oops
    for(const osg::Callback* cc = getCullCallback(); cc; cc = cc->getNestedCallback())
        cc->releaseGLObjects(state);
    for(const osg::Callback* uc = getUpdateCallback(); uc; uc = uc->getNestedCallback())
        uc->releaseGLObjects(state);
    for(const osg::Callback* ec = getEventCallback(); ec; ec = ec->getNestedCallback())
        ec->releaseGLObjects(state);

    // run this prior to the static bins and pools.
    osg::Group::releaseGLObjects(state);

    // indirect rendering objects
    ChonkRenderBin::releaseSharedGLObjects(state);

    // managed GL objects
    GLObjectPool::releaseGLObjects(state);
}

ClampingManager*
MapNode::getClampingManager()
{
    return _clampingManager;
}

bool
MapNode::getGeoPointUnderMouse(
    osg::View* view,
    float mx, float my,
    GeoPoint& output) const
{
    osg::Vec3d world;
    if (getTerrain()->getWorldCoordsUnderMouse(view, mx, my, world))
    {
        return output.fromWorld(getMapSRS(), world);
    }
    return false;
}

GeoPoint
MapNode::getGeoPointUnderMouse(
    osg::View* view,
    float mx, float my) const
{
    GeoPoint p;
    getGeoPointUnderMouse(view, mx, my, p);
    return p;
}
