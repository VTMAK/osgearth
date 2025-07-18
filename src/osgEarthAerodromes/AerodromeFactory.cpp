/* -*-c++-*- */
/* osgEarth - Dynamic map generation toolkit for OpenSceneGraph
* Copyright 2008-2014 Pelican Mapping
* http://osgearth.org
*
* osgEarth is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>
*/

#include "AerodromeFactory"
#include "Common"
#include "AerodromeNode"
#include "AerodromeCatalog"
#include "BoundaryNode"
#include "LightBeaconNode"
#include "LightIndicatorNode"
#include "LinearFeatureNode"
#include "PavementNode"
#include "RunwayNode"
#include "RunwayThresholdNode"
#include "StartupLocationNode"
#include "StopwayNode"
#include "TaxiwayNode"
#include "TaxiwaySignNode"
#include "TerminalNode"
#include "WindsockNode"
#include "AerodromeRenderer"

#include <osgEarth/Feature>
#include <osgEarth/FeatureSource>
#include <osgEarth/FeatureCursor>
#include <osgEarth/Registry>
#include <osgEarth/HTM>
#include <osgEarth/ElevationRanges>
#include <osgEarth/Utils>

#include <osgDB/FileNameUtils>
#include <osgDB/Registry>
#include <osgUtil/Optimizer>

#define USE_PAGING_MANAGER
#ifdef USE_PAGING_MANAGER
#include <osgEarth/PagedNode>
#else
#include <osg/PagedLOD>
#endif

using namespace osgEarth;
using namespace osgEarth::Aerodrome;
using namespace osgEarth::Contrib;
using namespace osgEarth::Util;

#define LC "[AerodromeFactory] "


namespace
{
    template <typename T, typename Y> class FeatureNodeFinder : public osg::NodeVisitor
    {
    public:
        FeatureNodeFinder(const std::string& attr, const std::string& value)
          : _attr(attr), _value(value), osg::NodeVisitor(osg::NodeVisitor::TRAVERSE_ALL_CHILDREN)
        {
        }

        void apply(T& node)
        {
            if (node.getFeature()->getString(_attr) == _value)
                _found = &node;
        }

        void apply(osg::Group& node)
        {
            if (dynamic_cast<AerodromeNode*>(&node) || dynamic_cast<Y*>(&node))
                traverse(node);
            else if (dynamic_cast<T*>(&node))
                apply(static_cast<T&>(node));
        }

    public:
        T* foundNode() { return _found.get(); }

    private:
        std::string _attr;
        std::string _value;
        osg::ref_ptr<T> _found;
    };
}


// -----------------------------------------------------------------------------
// pseudo-loader for paging in aerodromes.

namespace
{
    UID                               _uid         = 0;
    Threading::ReadWriteMutex         _amfMutex;
    typedef std::map<UID, osg::observer_ptr<AerodromeFactory> > AMFRegistry;
    AMFRegistry _amfRegistry;

    static std::string s_makeURI( UID uid, const std::string& icao ) 
    {
        std::stringstream buf;
        buf << uid << "." << icao << ".osgearth_pseudo_amf";
        std::string str;
        str = buf.str();
        return str;
    }

    static std::string s_makeQuery(const std::string& field, const std::string& icao ) 
    {
        std::stringstream buf;
        buf << field << "='" << icao << "'";
        std::string str;
        str = buf.str();
        return str;
    }
}

#ifndef USE_PAGING_MANAGER
/**
 * A pseudo-loader for paged feature tiles.
 */
struct osgEarthAerodromeModelPseudoLoader : public osgDB::ReaderWriter
{
    osgEarthAerodromeModelPseudoLoader()
    {
        supportsExtension( "osgearth_pseudo_amf", "Aerodrome model pseudo-loader" );
    }

    const char* className() const
    { // override
        return "osgEarth Aerodrome Model Pseudo-Loader";
    }

    ReadResult readNode(const std::string& uri, const Options* options) const
    {
        if ( !acceptsExtension( osgDB::getLowerCaseFileExtension(uri) ) )
            return ReadResult::FILE_NOT_HANDLED;

        UID uid;
        char icao[11];
        sscanf( uri.c_str(), "%d.%10[^'.'].%*s", &uid, icao );

        osg::ref_ptr<AerodromeFactory> factory = getFactory(uid);
        if ( factory.valid() )
        {
            Registry::instance()->startActivity(uri);

            AerodromeNode* node = factory->getAerodromeNode(std::string(icao));

            Registry::instance()->endActivity(uri);

            if (node)
            {
                if (factory->getSceneGraphCallbacks())
                    factory->getSceneGraphCallbacks()->firePreMergeNode(node);

                return ReadResult(node);
            }
        }

        return ReadResult::ERROR_IN_READING_FILE;
    }

    static UID registerFactory( AerodromeFactory* factory )
    {
        Threading::ScopedWriteLock lock( _amfMutex );
        UID key = ++_uid;
        _amfRegistry[key] = factory;
        OE_DEBUG << "Registered AMF " << key << std::endl;
        return key;
    }

    static void unregisterFactory( UID uid )
    {
        Threading::ScopedWriteLock lock( _amfMutex );
        _amfRegistry.erase( uid );
        OE_DEBUG << "Unregistered AMF " << uid << std::endl;
    }

    static AerodromeFactory* getFactory( UID uid ) 
    {
        Threading::ScopedReadLock lock( _amfMutex );
        AMFRegistry::const_iterator i = _amfRegistry.find( uid );
        return i != _amfRegistry.end() ? i->second.get() : 0L;
    }
};

REGISTER_OSGPLUGIN(osgearth_pseudo_amf, osgEarthAerodromeModelPseudoLoader);
#endif


AerodromeFactory::RendererFactory AerodromeFactory::s_renderer_factory = nullptr;

AerodromeFactory::AerodromeFactory(
    const Map* map,
    AerodromeCatalog* catalog,
    SceneGraphCallbacks* callbacks) :
    
    _map(map),
    _catalog(catalog),
    _sceneGraphCallbacks(callbacks),
    _residentTiles(std::make_shared<std::atomic_int>(0))
{
    //nop
}

void
AerodromeFactory::init(const osgDB::Options* options)
{
#ifndef USE_PAGING_MANAGER
    _uid = osgEarthAerodromeModelPseudoLoader::registerFactory( this );
#endif

    _dbOptions = options;

    // setup the PagedLODs
    seedAerodromes(_catalog.get(), _dbOptions.get());
}

AerodromeFactory::~AerodromeFactory()
{
#ifndef USE_PAGING_MANAGER
    osgEarthAerodromeModelPseudoLoader::unregisterFactory( _uid );
#endif
}

template <typename T, typename Y, typename P>
void AerodromeFactory::createFeatureNodes(P featureOpts, AerodromeNode* aerodrome, const osgDB::Options* options, void (*processor)(T* node, AerodromeNode* aerodrome))
{
    Status s = featureOpts.featureSource().open(options);
    if (s.isError() || !featureOpts.featureSource().getLayer())
    {
        OE_WARN << LC << "Skipping feature source; open failed: " << s.message() << std::endl;
        return;
    }
    FeatureSource* featureSource = featureOpts.featureSource().getLayer();

    Y* parentGroup = new Y();
    aerodrome->addChild(parentGroup);

    OE_DEBUG << LC << "Reading features...\n";

    int featureCount = 0;

    Query query;
    query.expression() = s_makeQuery(featureOpts.icaoAttr().value(), aerodrome->icao());

    osg::ref_ptr<FeatureCursor> cursor = featureSource->createFeatureCursor(query);
    while ( cursor.valid() && cursor->hasMore() )
    {
        Feature* f = cursor->nextFeature();

        /* **************************************** */
        /* Necessary but not sure why               */

        const SpatialReference* ecefSRS = f->getSRS()->getGeocentricSRS();

        /* **************************************** */

        OE_DEBUG << LC << "Adding feature to aerodrome: " << aerodrome->icao() << std::endl;

        // create new node
        //T* tNode = new T(featureOpts, aerodrome->icao(), featureSource, f->getFID());
        T* tNode = new T(featureOpts, aerodrome->icao(), f);

        // if a processor function is passed in, call it
        if (processor)
            (*processor)(tNode, aerodrome);

        // add the new node to the parent AerodromeNode
        parentGroup->addChild(tNode);
        featureCount++;
    }

    OE_DEBUG << LC << "Added " << featureCount << " feature nodes to aerodrome " << aerodrome->icao() << std::endl;
}

template <typename T, typename Y, typename P>
void AerodromeFactory::createMergedFeatureNodes(P featureOpts, AerodromeNode* aerodrome, const osgDB::Options* options, void (*processor)(T* node, AerodromeNode* aerodrome))
{
    Status s = featureOpts.featureSource().open(options);
    if (s.isError() || featureOpts.featureSource().getLayer() == 0L)
    {
        OE_WARN << LC << "Skipping boundary source; open failed: " << s.message() << std::endl;
        return;
    }

    FeatureSource* featureSource = featureOpts.featureSource().getLayer();

    Y* parentGroup = new Y();
    aerodrome->addChild(parentGroup);

    OE_DEBUG << LC << "Reading features...\n";

    int featureCount = 0;

    Query query;
    query.expression() = s_makeQuery(featureOpts.icaoAttr().value(), aerodrome->icao());


    osg::ref_ptr<MultiGeometry> mg = new MultiGeometry();
    const SpatialReference* srs = 0L;

    osg::ref_ptr<Feature> newFeature = 0L;

    osg::ref_ptr<FeatureCursor> cursor = featureSource->createFeatureCursor(query);
    while ( cursor.valid() && cursor->hasMore() )
    {
        Feature* f = cursor->nextFeature();

        /* **************************************** */
        /* Necessary but not sure why               */

        const SpatialReference* ecefSRS = f->getSRS()->getGeocentricSRS();

        /* **************************************** */

        OE_DEBUG << LC << "Adding feature to aerodrome: " << aerodrome->icao() << std::endl;

        if ( !newFeature.valid() )
        {
            newFeature = new Feature( *f );
            newFeature->setGeometry( mg.get() );
        }

        mg->add( f->getGeometry() );
    }

    if ( newFeature.valid() )
    {
        T* tNode = new T(featureOpts, aerodrome->icao(), newFeature.get());

        // if a processor function is passed in, call it
        if (processor)
            (*processor)(tNode, aerodrome);

        // add the new node to the parent AerodromeNode
        parentGroup->addChild(tNode);
        featureCount++;
    }

    OE_DEBUG << LC << "Added " << featureCount << " feature nodes to aerodrome " << aerodrome->icao() << std::endl;
}



void AerodromeFactory::createBoundaryNodes(BoundaryFeatureOptions boundaryOpts, AerodromeNode* aerodrome, const osgDB::Options* options)
{
    if (!boundaryOpts.featureSource().getLayer())
    {
        OE_WARN << LC << "Cannot create boundary features: feature source is not set." << std::endl;
        return;
    }

    if (!aerodrome)
    {
        OE_WARN << LC << "Cannot create boundary features: AerodromeNode is not set." << std::endl;
        return;
    }

    Status status = boundaryOpts.featureSource().open(options);
    if (status.isError() || boundaryOpts.featureSource().getLayer() == 0L)
    {
        OE_WARN << LC << "No feature data: " << status.message() << std::endl;
        return;
    }
    //osg::ref_ptr<FeatureSource> featureSource = FeatureSource::create(boundaryOpts.featureOptions().value());
    //featureSource->setReadOptions(options);
    //const Status& status = featureSource->open();
    FeatureSource* featureSource = boundaryOpts.featureSource().getLayer();
    
    Query query;
    query.expression() = s_makeQuery(boundaryOpts.icaoAttr().value(), aerodrome->icao());

    osg::ref_ptr<FeatureCursor> cursor = featureSource->createFeatureCursor(query);
    while ( cursor.valid() && cursor->hasMore() )
    {
        Feature* f = cursor->nextFeature();

        /* **************************************** */
        /* Necessary but not sure why               */

        const SpatialReference* ecefSRS = f->getSRS()->getGeocentricSRS();

        /* **************************************** */

        OE_DEBUG << LC << "Adding boundary to aerodrome: " << aerodrome->icao() << std::endl;

        // create new node and add to parent AerodromeNode
        aerodrome->setBoundary(new BoundaryNode(boundaryOpts, aerodrome->icao(), f));
        break;
    }
}

void AerodromeFactory::processStopwayNode(StopwayNode* stopway, AerodromeNode* aerodrome)
{
    if (stopway)
    {
        std::string rwyNum = stopway->getFeature()->getString("rwy_num");
            
        FeatureNodeFinder<RunwayNode, RunwayGroup> finder("rwy_num1", rwyNum);
        aerodrome->accept(finder);

        osg::ref_ptr<RunwayNode> runway = finder.foundNode();
        if (!runway.valid())
        {
            FeatureNodeFinder<RunwayNode, RunwayGroup> finder2("rwy_num2", rwyNum);
            aerodrome->accept(finder2);
            runway = finder2.foundNode();
        }
            
        if (runway.valid())
            stopway->setReferencePoint(runway->getFeature()->getGeometry()->getBounds().center());
        else
            OE_WARN << LC << "Could not find runway " << rwyNum << " for stopway." << std::endl;
    }
}

AerodromeNode*
AerodromeFactory::createAerodrome(
    AerodromeCatalog* catalog, 
    const std::string& icao, 
    bool buildTerminals,
    const osgDB::Options* options)
{
    osg::ref_ptr<AerodromeNode> aerodrome = new AerodromeNode(icao);

    for(BoundaryOptionsSet::const_iterator i = catalog->boundaryOptions().begin(); i != catalog->boundaryOptions().end(); ++i)
        AerodromeFactory::createBoundaryNodes(*i, aerodrome.get(), options);
    
    for(AerodromeOptionsSet::const_iterator i = catalog->pavementOptions().begin(); i != catalog->pavementOptions().end(); ++i)
        //AerodromeFactory::createMergedFeatureNodes<PavementNode, PavementGroup, AerodromeFeatureOptions>(*i, aerodrome.get(), options);
        AerodromeFactory::createFeatureNodes<PavementNode, PavementGroup, AerodromeFeatureOptions>(*i, aerodrome.get(), options);
    
    for(AerodromeOptionsSet::const_iterator i = catalog->taxiwayOptions().begin(); i != catalog->taxiwayOptions().end(); ++i)
        AerodromeFactory::createMergedFeatureNodes<TaxiwayNode, TaxiwayGroup, AerodromeFeatureOptions>(*i, aerodrome.get(), options);

    for(AerodromeOptionsSet::const_iterator i = catalog->runwayOptions().begin(); i != catalog->runwayOptions().end(); ++i)
        AerodromeFactory::createFeatureNodes<RunwayNode, RunwayGroup, AerodromeFeatureOptions>(*i, aerodrome.get(), options);

    for(AerodromeOptionsSet::const_iterator i = catalog->runwayThresholdOptions().begin(); i != catalog->runwayThresholdOptions().end(); ++i)
        AerodromeFactory::createFeatureNodes<RunwayThresholdNode, RunwayThresholdGroup, AerodromeFeatureOptions>(*i, aerodrome.get(), options);

    for(AerodromeOptionsSet::const_iterator i = catalog->stopwayOptions().begin(); i != catalog->stopwayOptions().end(); ++i)
        AerodromeFactory::createFeatureNodes<StopwayNode, StopwayGroup, AerodromeFeatureOptions>(*i, aerodrome.get(), options, &AerodromeFactory::processStopwayNode);

    for(AerodromeOptionsSet::const_iterator i = catalog->linearFeatureOptions().begin(); i != catalog->linearFeatureOptions().end(); ++i)
        AerodromeFactory::createMergedFeatureNodes<LinearFeatureNode, LinearFeatureGroup, AerodromeFeatureOptions>(*i, aerodrome.get(), options);

    for(AerodromeOptionsSet::const_iterator i = catalog->startupLocationOptions().begin(); i != catalog->startupLocationOptions().end(); ++i)
        AerodromeFactory::createFeatureNodes<StartupLocationNode, StartupLocationGroup, AerodromeFeatureOptions>(*i, aerodrome.get(), options);

    for(AerodromeOptionsSet::const_iterator i = catalog->lightBeaconOptions().begin(); i != catalog->lightBeaconOptions().end(); ++i)
        AerodromeFactory::createFeatureNodes<LightBeaconNode, LightBeaconGroup, AerodromeFeatureOptions>(*i, aerodrome.get(), options);

    for(AerodromeOptionsSet::const_iterator i = catalog->lightIndicatorOptions().begin(); i != catalog->lightIndicatorOptions().end(); ++i)
         AerodromeFactory::createFeatureNodes<LightIndicatorNode, LightIndicatorGroup, AerodromeFeatureOptions>(*i, aerodrome.get(), options);

    for (AerodromeOptionsSet::const_iterator i = catalog->taxiWayLightOptions().begin(); i != catalog->taxiWayLightOptions().end(); ++i)
        AerodromeFactory::createFeatureNodes<LightIndicatorNode, LightIndicatorGroup, AerodromeFeatureOptions>(*i, aerodrome.get(), options);

    for(AerodromeOptionsSet::const_iterator i = catalog->runwayAccessLightOptions().begin(); i != catalog->runwayAccessLightOptions().end(); ++i)
        AerodromeFactory::createFeatureNodes<LightIndicatorNode, LightIndicatorGroup, AerodromeFeatureOptions>(*i, aerodrome.get(), options);

    for(AerodromeOptionsSet::const_iterator i = catalog->taxiwaySignOptions().begin(); i != catalog->taxiwaySignOptions().end(); ++i)
        AerodromeFactory::createFeatureNodes<TaxiwaySignNode, TaxiwaySignGroup, AerodromeFeatureOptions>(*i, aerodrome.get(), options);

    for(AerodromeOptionsSet::const_iterator i = catalog->windsockOptions().begin(); i != catalog->windsockOptions().end(); ++i)
        AerodromeFactory::createFeatureNodes<WindsockNode, WindsockGroup, AerodromeFeatureOptions>(*i, aerodrome.get(), options);

    if (buildTerminals)
    {
        for (TerminalOptionsSet::const_iterator i = catalog->terminalOptions().begin(); i != catalog->terminalOptions().end(); ++i)
            AerodromeFactory::createFeatureNodes<TerminalNode, TerminalGroup, TerminalFeatureOptions>(*i, aerodrome.get(), options);
    }

    return aerodrome.release();
}

AerodromeNode* AerodromeFactory::getAerodromeNode(const std::string& icao) const
{
    // create and initialize a renderer
    osg::ref_ptr<AerodromeRenderer> renderer = s_renderer_factory ? s_renderer_factory() : new AerodromeRenderer();
    if (!renderer.valid())
        return nullptr;

    renderer->setStyleSheet(_catalog->styleSheet());
    renderer->setClampToTerrain(_clampToTerrain);
    renderer->initialize(_map.get(), _dbOptions.get());

    Threading::scoped_lock_if(_serialMutex, _serialLoading);
    //Threading::ScopedWriteLock lock(_mutex);

    OE_START_TIMER(getAerodromeNode);

    // create AerodromeNode
    OE_START_TIMER(create);

    osg::ref_ptr<AerodromeNode> node = createAerodrome(_catalog.get(), icao, _buildTerminals, _dbOptions.get());

    node->accept(*renderer.get());

    // Generate shaders (using a cache)
    osg::ref_ptr<osgEarth::StateSetCache> cache = new osgEarth::StateSetCache();
    osgEarth::ShaderGenerator gen;
    gen.run(node.get(), "Aerodrome", cache.get());
    //osgEarth::Registry::shaderGenerator().run( node.get(), "Aerodrome", cache.get() );

    float s = OE_STOP_TIMER(getAerodromeNode);
    OE_INFO << LC << std::setprecision(3) << "Built \"" << icao << "\" - duration=" << s << "s" << std::endl;

    if (node.valid() && node->getBound().valid())
    {
        node->getOrCreateUserDataContainer()->addUserObject(
            new osgEarth::Util::TrackerTag(_residentTiles));
    }

    return node.release();
}

void 
AerodromeFactory::setRendererFactory(RendererFactory value)
{
    s_renderer_factory = value;
}

void
AerodromeFactory::seedAerodromes(AerodromeCatalog* catalog, const osgDB::Options* options)
{
    osg::ref_ptr<const Map> refMap;
    if (_map.lock(refMap)==false)
    {
        return;
    }
    removeChildren(0, getNumChildren());

    // set up a spatial indexing tree
    HTMGroup* tree = new HTMGroup();
    tree->setMaximumObjectsPerCell(4);
    tree->setThreshold(_lodRange);
    //tree->setMaxRange(_lodRange);

    this->addChild(tree);

    OE_INFO << LC << "Seeding aerodromes from boundaries." << std::endl;

    int aeroCount = 0;

    unsigned maxRangesLevel = ElevationRanges::getMaxLevel();

    for(BoundaryOptionsSet::iterator i = catalog->boundaryOptions().begin(); i != catalog->boundaryOptions().end(); ++i)
    {
        Status s = i->featureSource().open(options);
        if (s.isError() || !i->featureSource().getLayer())
        {
            OE_WARN << LC << "Skipping boundary source; failed to open feature source \"" 
                << i->featureSource().embeddedOptions()->name().get() << "\" (error=" << s.message() << ")" << std::endl;
            continue;
        }
        FeatureSource* featureSource = i->featureSource().getLayer();

        osg::ref_ptr<FeatureCursor> cursor = featureSource->createFeatureCursor();
        while ( cursor.valid() && cursor->hasMore() )
        {
            Feature* f = cursor->nextFeature();

            std::string icao = f->getString(i->icaoAttr().value());

            // skip the aerodrome if it is in the exclusions list..
            if (!icao.empty() && _icaoExclusions.count(toLower(icao)) == 0)
            {
                // create PagedLOD for this aerodrome
                std::string uri = s_makeURI( _uid, icao );

                if (f->getGeometry())
                {
                    GeoCircle circle = f->getExtent().computeBoundingGeoCircle();
                    GeoPoint center = circle.getCenter();

#if 1 // Attempt to get the bounding sphere more accurate
                    GeoPoint centerER = center.transform(ElevationRanges::getProfile()->getSRS());
                    TileKey rangeKey = ElevationRanges::getProfile()->createTileKey(centerER.x(), centerER.y(), maxRangesLevel);
                    short lo, hi;
                    ElevationRanges::getElevationRange(maxRangesLevel, rangeKey.getTileX(), rangeKey.getTileY(), lo, hi);
                    center.z() = (lo + hi) / 2.0;
#endif
                    center.transformInPlace(refMap->getSRS());
                    osg::Vec3d centerWorld;
                    center.toWorld(centerWorld);

                    PagedNode2* p = new PagedNode2();
                    p->setName(icao);
                    p->setSceneGraphCallbacks(getSceneGraphCallbacks());
                    p->setCenter(centerWorld);
                    
                    // Radius needs to encompass any possible location of the aerodrome within the tile's Z range:
                    p->setRadius(circle.getRadius() + 0.5 * (float)(hi - lo));

                    osg::ref_ptr<const osgDB::Options> load_options(options);
                    osg::observer_ptr<AerodromeFactory> factory_obs(this);

                    p->setLoadFunction([uri, factory_obs, icao](Cancelable* c) mutable
                        {
                            osg::ref_ptr<osg::Node> result;
                            osg::ref_ptr<AerodromeFactory> factory(factory_obs);
                            if (factory.valid())
                                result = factory->getAerodromeNode(std::string(icao));
                            return result;
                        }
                    );
                    
                    if (_lodRange.getUnits() == Units::PIXELS)
                    {
                        p->setMinPixels(_lodRange.getValue());
                    }
                    else
                    {
                        p->setMaxRange(_lodRange.as(Units::METERS));
                    }

                    tree->addChild( p );
                    aeroCount++;
                }
                else
                {
                    OE_DEBUG << LC << "Skipping boundary feature: no geometry." << std::endl;
                }
            }
        }
    }

    OE_DEBUG << LC << aeroCount << " aerodromes found and seeded." << std::endl;
}