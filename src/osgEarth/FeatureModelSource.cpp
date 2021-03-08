/* -*-c++-*- */
/* osgEarth - Geospatial SDK for OpenSceneGraph
 * Copyright 2020 Pelican Mapping
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
#include <osgEarth/FeatureModelSource>
#include <osgEarth/FeatureModelGraph>
#include <osgEarth/SpatialReference>
#include <osgEarth/ShaderFactory>
#include <osgEarth/ShaderUtils>
#include <osgEarth/Registry>
#include <osgEarth/Capabilities>
#include <osgEarth/DrapeableNode>
#include <osgEarth/ClampableNode>
#include <osgEarth/GLUtils>
#include <osgEarth/Progress>
#include <osg/Notify>

using namespace osgEarth;

#ifndef GL_CLIP_DISTANCE0
#define GL_CLIP_DISTANCE0 0x3000
#endif

#define LC "[FeatureModelSource] "

//........................................................................

FeatureModelOptions::FeatureModelOptions(const ConfigOptions& co) :
_lit               ( true ),
_maxGranularity_deg( 1.0 ),
_clusterCulling    ( false ),
_backfaceCulling   ( true ),
_alphaBlending     ( true ),
_sessionWideResourceCache( true ),
_nodeCaching(false)
{
    fromConfig(co.getConfig());
}

void
FeatureModelOptions::fromConfig(const Config& conf)
{
    styleSheet().get(conf, "styles");

    conf.get( "layout",           _layout );
    conf.get( "fading",           _fading );
    conf.get( "feature_name",     _featureNameExpr );
    conf.get( "feature_indexing", _featureIndexing );

    conf.get( "lighting",         _lit );
    conf.get( "max_granularity",  _maxGranularity_deg );
    conf.get( "cluster_culling",  _clusterCulling );
    conf.get( "backface_culling", _backfaceCulling );
    conf.get( "alpha_blending",   _alphaBlending );
    conf.get( "node_caching",     _nodeCaching );
    
    conf.get( "session_wide_resource_cache", _sessionWideResourceCache );

    const Config& filtersConf = conf.child("filters");
    for(ConfigSet::const_iterator i = filtersConf.children().begin(); i != filtersConf.children().end(); ++i)
        filters().push_back( ConfigOptions(*i) );

    // shorthand for enabling feature indexing
    if (featureIndexing().isSet() == false && conf.value("pickable", false) == true)
    {
        featureIndexing()->enabled() = true;
    }
}

Config
FeatureModelOptions::getConfig() const
{
    Config conf;

    styleSheet().set(conf, "styles");

    conf.set( "layout",           _layout );
    conf.set( "fading",           _fading );
    conf.set( "feature_name",     _featureNameExpr );
    conf.set( "feature_indexing", _featureIndexing );

    conf.set( "lighting",         _lit );
    conf.set( "max_granularity",  _maxGranularity_deg );
    conf.set( "cluster_culling",  _clusterCulling );
    conf.set( "backface_culling", _backfaceCulling );
    conf.set( "alpha_blending",   _alphaBlending );
    conf.set( "node_caching",     _nodeCaching );
    
    conf.set( "session_wide_resource_cache", _sessionWideResourceCache );

    if (filters().empty() == false)
    {
        Config temp;
        for(unsigned i=0; i<filters().size(); ++i)
            temp.add( filters()[i].getConfig() );
        conf.set( "filters", temp );
    }

    return conf;
}

//------------------------------------------------------------------------

osg::Group*
FeatureNodeFactory::getOrCreateStyleGroup(const Style& style,
                                          Session*     session)
{
    osg::Group* group = 0L;

    // If we're draping, the style group will be a DrapeableNode.
    const AltitudeSymbol* alt = style.get<AltitudeSymbol>();
    if (alt &&
        alt->clamping() == AltitudeSymbol::CLAMP_TO_TERRAIN &&
        alt->technique() == AltitudeSymbol::TECHNIQUE_DRAPE )
    {
        group = new DrapeableNode();
    }

    else if (alt &&
        alt->clamping() == alt->CLAMP_TO_TERRAIN &&
        alt->technique() == alt->TECHNIQUE_GPU)
    {
        group = new ClampableNode();
    }

    // Otherwise, a normal group.
    if ( !group )
    {
        group = new osg::Group();
    }

    // apply necessary render styles.
    const RenderSymbol* render = style.get<RenderSymbol>();
    if ( render )
    {
        if ( render->depthTest().isSet() )
        {
            group->getOrCreateStateSet()->setMode(
                GL_DEPTH_TEST, 
                (render->depthTest() == true ? osg::StateAttribute::ON : osg::StateAttribute::OFF) | osg::StateAttribute::OVERRIDE );
        }

        if ( render->lighting().isSet() )
        {
            osg::StateSet* stateset = group->getOrCreateStateSet();

            GLUtils::setLighting(
                stateset,
                (render->lighting() == true ? osg::StateAttribute::ON : osg::StateAttribute::OFF) | osg::StateAttribute::OVERRIDE );
        }

        if ( render->backfaceCulling().isSet() )
        {
            group->getOrCreateStateSet()->setMode(
                GL_CULL_FACE,
                (render->backfaceCulling() == true ? osg::StateAttribute::ON : osg::StateAttribute::OFF) | osg::StateAttribute::OVERRIDE );
        }

#if !(defined(OSG_GLES2_AVAILABLE) || defined(OSG_GLES3_AVAILABLE))
        if ( render->clipPlane().isSet() )
        {
            GLenum mode = GL_CLIP_DISTANCE0 + (render->clipPlane().value());
            group->getOrCreateStateSet()->setMode(mode, 1);
        }
#endif

        if ( render->minAlpha().isSet() )
        {
            DiscardAlphaFragments().install( group->getOrCreateStateSet(), render->minAlpha().value() );
        }
    }

    return group;
}


//------------------------------------------------------------------------


GeomFeatureNodeFactory::GeomFeatureNodeFactory( const GeometryCompilerOptions& options ) : 
_options( options ) 
{ 
    //nop
}

bool GeomFeatureNodeFactory::createOrUpdateNode(
    FeatureCursor*            features,
    const Style&              style,
    const FilterContext&      context,
    osg::ref_ptr<osg::Node>&  node,
    const Query&              query
)
{
    GeometryCompiler compiler( _options );
    node = compiler.compile( features, style, context );
    return node.valid();
}




// VRV - for backwards compatibility

FeatureModelSourceOptions::FeatureModelSourceOptions(const ConfigOptions& options) :
    ModelSourceOptions(options),
    FeatureModelOptions(options)
{
    fromConfig( _conf );
}

void
FeatureModelSourceOptions::fromConfig( const Config& conf )
{
    featureSource().get(conf, "features");
}

Config
FeatureModelSourceOptions::getConfig() const
{
    Config conf = ModelSourceOptions::getConfig();
    conf.merge(FeatureModelOptions::getConfig());
    featureSource().set(conf, "features");
    return conf;
}

//------------------------------------------------------------------------

FeatureModelSource::FeatureModelSource( const FeatureModelSourceOptions& options ) :
    ModelSource( options ),
    _options   ( options )
{
    //nop
}

void
FeatureModelSource::setFeatureSource( FeatureSource* source )
{
    if ( !_options.featureSource().getLayer())
    {
        _options.featureSource().setLayer(source);
    }
    else
    {
        OE_WARN << LC << "Illegal: cannot set a feature source after one is already set" << std::endl;
    }
}

Status
FeatureModelSource::initialize(const osgDB::Options* readOptions)
{
    if (readOptions)
        setReadOptions(readOptions);

    Status fs = _options.featureSource().open(readOptions);
    if (fs.isError())
        return fs;

    Status s =_options.styleSheet().open(readOptions);
    if (s.isError())
        return s;

    // Try to fill the DataExtent list using the FeatureProfile
    const FeatureProfile* featureProfile = getFeatureSource()->getFeatureProfile();
    if (featureProfile == NULL)
        return Status::Error("Failed to establish a feature profile");

    if (featureProfile->getTilingProfile() != NULL)
    {
        // Use specified profile's GeoExtent
        getDataExtents().push_back(DataExtent(featureProfile->getTilingProfile()->getExtent()));
    }
    else if (featureProfile->getExtent().isValid() == true)
    {
        // Use FeatureProfile's GeoExtent
        getDataExtents().push_back(DataExtent(featureProfile->getExtent()));
    }

    return Status::OK();
}

void
FeatureModelSource::setReadOptions(const osgDB::Options* readOptions)
{
    _readOptions = Registry::cloneOrCreateOptions(readOptions);

    // for texture atlas support
    _readOptions->setObjectCacheHint(osgDB::Options::CACHE_IMAGES);

    if (_options.featureSource().getLayer())
    {
        _options.featureSource().getLayer()->setReadOptions(_readOptions.get());
    }
}

osg::Node*
FeatureModelSource::createNodeImplementation(const Map* map,  ProgressCallback* progress )
{
    // trivial bailout.
    if (!getStatus().isOK())
        return 0L;

    // user must provide a valid map.
    if ( !map )
    {
        OE_WARN << LC << "NULL Map is illegal when building feature data." << std::endl;
        return 0L;
    }

    // make sure the feature source initialized properly:
    if ( !_options.featureSource().getLayer() || !_options.featureSource().getLayer()->getFeatureProfile() )
    {
        return 0L;
    }

    // create a feature node factory:
    FeatureNodeFactory* factory = createFeatureNodeFactory();
    if ( !factory )
    {
        OE_WARN << LC << "Unable to create a feature node factory!" << std::endl;
        setStatus(Status::Error(Status::ServiceUnavailable, "Failed to create a feature node factory"));
        return 0L;
    }

    // Session holds data that's shared across the life of the FMG
    Session* session = new Session( 
        map, 
        _options.styleSheet().getLayer(),
        _options.featureSource().getLayer(),
        _readOptions.get() );

    // Name the session (for debugging purposes)
    session->setName( this->getName() );

    // Graph that will render feature models. May included paged data.
    FeatureModelGraph* graph = new FeatureModelGraph(_options);
    graph->setSession(session);
    graph->setNodeFactory(factory);
    graph->setSceneGraphCallbacks(getSceneGraphCallbacks());
    graph->setStyleSheet(_options.styleSheet().getLayer());
    graph->open();

    return graph;
}