/* -*-c++-*- */
/* osgEarth - Dynamic map generation toolkit for OpenSceneGraph
 * Copyright 2008-2016 Pelican Mapping
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
#include "BuildingFactory"
#include "BuildingSymbol"
#include "BuildingVisitor"
#include "BuildContext"

#include <osgEarth/Session>
#include <osgEarth/AltitudeFilter>
#include <osgEarth/Geometry>
#include <osgEarth/ResourceLibrary>
#include <osgEarth/StyleSheet>
#include <memory>

using namespace osgEarth;
using namespace osgEarth::Buildings;

#define LC "[BuildingFactory] "

BuildingFactory::BuildingFactory():
    _parapets(true),
    _clutter(true)
{
    setSession( new Session(0L) );
}

void
BuildingFactory::setSession(Session* session)
{
    _session = session;
}

bool
BuildingFactory::cropToCentroid(const Feature* feature, const GeoExtent& extent) const
{
    if ( !extent.isValid() )
        return true;

    // make sure the centroid is in the crop-to extent: 
    GeoPoint centroid( feature->getSRS(), feature->getGeometry()->getBounds().center() );
    return extent.contains(centroid);
}

bool
BuildingFactory::create(Feature*               feature,
                        const GeoExtent&       cropTo,
                        ElevationPool::Envelope& envelope,
                        const Style*           style,
                        BuildingVector&        output,
                        const osgDB::Options*  readOptions,
                        ProgressCallback*      progress)
{
    if ( !feature || !feature->getGeometry() )
        return false;

    // TODO: for the single-feature version here, consider a token or other way
    // to compute all this common stuff up front. This was not necessary for the
    // FeatureCursor variation. -gw

    bool needToClamp =
        style &&
        style->has<AltitudeSymbol>() &&
        style->get<AltitudeSymbol>()->clamping() != AltitudeSymbol::CLAMP_NONE;

    // Find the building symbol if there is one; this will tell us how to 
    // resolve building heights, among other things.
    const BuildingSymbol* buildingSymbol =
        style ? style->get<BuildingSymbol>() :
        _session->styles() ? _session->styles()->getDefaultStyle()->get<BuildingSymbol>() :
        0L;        
  
    // Pull a resource library if one is defined.
    ResourceLibrary* reslib = 0L;
    if (buildingSymbol && buildingSymbol->library().isSet())
    {
        auto library = buildingSymbol->library()->eval(feature, FilterContext(_session.get()));
        reslib = _session->styles()->getResourceLibrary(library);
    }
    if ( !reslib )
    {
        reslib = _session->styles()->getDefaultResourceLibrary();
    }

    // Construct a context to use during the build process.
    BuildContext context;
    context.setDBOptions( readOptions );
    context.setResourceLibrary( reslib );
    context.setClutter(_clutter);
    context.setParapets(_parapets);

    // URI context for external models
    URIContext uriContext( readOptions );

    if ( progress && progress->isCanceled() )
    {
        progress->message() = "in BuildingFactory::create";
        return false;
    }

    // resolve selection values from the symbology:
    optional<URI> externalModelURI;
    float         height    = 0.0f;
    unsigned      numFloors = 0u;
    TagSet     tags;

    FilterContext cx(_session.get());

    if ( buildingSymbol )
    {
        // see if we are referencing an external model.
        if ( buildingSymbol->modelURI().isSet())
        {
            std::string modelStr = buildingSymbol->modelURI()->eval(feature, cx);

            //VRV_PATCH
            //Sometimes the ductape engine returns undefined for the modelStr.
            // lets set the modelStr to "" and ignore those cases as well
            if (modelStr == "undefined")
            {
               modelStr = "";
            }
            //END VRV_PATCH
            if (!modelStr.empty())
            {
                externalModelURI = URI(modelStr, uriContext);
            }
        }

        // calculate height from expression. We do this first because
        // a height of zero will cause us to skip the feature altogether.
        if ( !externalModelURI.isSet() && buildingSymbol->height().isSet())
        {
            height = (float)buildingSymbol->height()->eval(feature, cx);
        
            if ( height > 0.0f )
            {
                // calculate tags from expression:
                if ( buildingSymbol->tags().isSet())
                {
                    std::string tagString = trim(buildingSymbol->tags()->eval(feature, cx));
                    if (!tagString.empty())
                    {
                        StringVector tagVector = StringTokenizer()
                            .whitespaceDelims()
                            .standardQuotes()
                            .keepEmpties(false)
                            .tokenize(tagString);
                        for (auto &i : tagVector)
                        {
                            tags.insert(osgEarth::toLower(i));
                        }
                    }
                }

                else
                {
                    // If there's no tag expression, the default is "building".
                    tags.insert("building");
                }
            }
        }
    }

    if ( height > 0.0f || externalModelURI.isSet() )
    {
        // Removing co-linear points will help produce a more "true"
        // longest-edge for rotation and roof rectangle calcuations.
        feature->getGeometry()->removeColinearPoints();

        // Transform the feature into the output SRS
        if ( _outSRS.valid() )
        {
            feature->transform( _outSRS.get() );
        }

        // this ensures that the feature's centroid is in our bounding
        // extent, so that a feature doesn't end up in multiple extents
        if ( !cropToCentroid(feature, cropTo) )
        {
            return true;
        }

        // Prepare for terrain clamping by finding the minimum and 
        // maximum elevations under the feature:
        bool terrainMinMaxValid = false;
        float min = FLT_MAX, max = -FLT_MAX;
        osg::ref_ptr<const Map> map = _session->getMap();

        if (needToClamp && feature->getGeometry() && map.valid())
        {
            std::vector<osg::Vec3d> points;
            points.reserve(feature->getGeometry()->getTotalPointCount());

            ConstGeometryIterator cgi(feature->getGeometry(), false);
            while (cgi.hasMore())
            {
                auto part = cgi.next();
                for (auto& i : part->asVector())
                    points.emplace_back(i);
            }
            envelope.sampleMapCoords(points.begin(), points.end(), progress);

            for(auto& i : points)
            {
                if (i.z() != NO_DATA_VALUE)
                {
                    min = std::min(min, (float)i.z());
                    max = std::max(max, (float)i.z());
                    terrainMinMaxValid = true;
                }
            }
        }
                
        context.setTerrainMinMax(
            terrainMinMaxValid ? min : 0.0f,
            terrainMinMaxValid ? max : 0.0f );

        // If this is an external model, set up a building referencing the model
        if ( externalModelURI.isSet() )
        {
            auto building = createExternalModelBuilding( feature, externalModelURI.get(), context );
            if ( building.has_value() )
            {
                output.emplace_back( std::move(*building) );
            }
        }

        // Otherwise, we are creating a parametric building:
        else
        {
            // If using a catalog, ask it to create one or more buildings for this feature:
            if ( _catalog.valid() )
            {
                float minHeight = terrainMinMaxValid ? max-min+3.0f : 3.0f;
                height = std::max( height, minHeight );
                _catalog->createBuildings(feature, tags, height, context, output, progress);
            }

            // Otherwise, create a simple one by hand:
            else
            {
                auto building = createBuilding(feature, progress);
                if ( building.has_value() )
                {
                    output.emplace_back( std::move(*building) );
                }
            }
        }
    }

    return true;
}

std::optional<Building>
BuildingFactory::createExternalModelBuilding(Feature*      feature,
                                             const URI&    modelURI,
                                             BuildContext& context)
{
    if ( !feature || modelURI.empty() )
        return std::nullopt;

    Geometry* geometry = feature->getGeometry();
    if ( !geometry || !geometry->isValid() )
        return std::nullopt;

    Building building;
    building.setExternalModelURI( modelURI );

    // Calculate a local reference frame for this building.
    // The frame clamps the building by including the terrain elevation that was passed in.
    osg::Vec3d center2d = geometry->getBounds().center();
    GeoPoint centerPoint( feature->getSRS(), center2d.x(), center2d.y(), context.getTerrainMin(), ALTMODE_ABSOLUTE );
    osg::Matrix local2world;
    centerPoint.createLocalToWorld( local2world );
    building.setReferenceFrame( local2world );

    return building;
}

std::optional<Building>
BuildingFactory::createBuilding(Feature* feature, ProgressCallback* progress)
{
    if ( feature == nullptr )
        return std::nullopt;

    std::optional<Building> building;

    Geometry* geometry = feature->getGeometry();

    if ( geometry && geometry->getComponentType() == Geometry::TYPE_POLYGON && geometry->isValid() )
    {
        // Calculate a local reference frame for this building:
        osg::Vec3d center2d = geometry->getBounds().center();
        GeoPoint centerPoint( feature->getSRS(), center2d.x(), center2d.y(), 0.0, ALTMODE_ABSOLUTE );
        osg::Matrix local2world, world2local;
        centerPoint.createLocalToWorld( local2world );
        world2local.invert( local2world );

        // Transform feature geometry into the local frame. This way we can do all our
        // building creation in cartesian, single-precision space.
        GeometryIterator iter(geometry, true);
        while(iter.hasMore())
        {
            Geometry* part = iter.next();
            for(Geometry::iterator i = part->begin(); i != part->end(); ++i)
            {
                osg::Vec3d world;
                feature->getSRS()->transformToWorld( *i, world );
                (*i) = world * world2local;
            }
        }

        BuildContext context;
        context.setSeed( feature->getFID() );
        context.setParapets(_parapets);
        context.setClutter(_clutter);

        // Next, iterate over the polygons and set up the Building object.
        GeometryIterator iter2( geometry, false );
        while(iter2.hasMore())
        {
            Polygon* polygon = dynamic_cast<Polygon*>(iter2.next());
            if ( polygon && polygon->isValid() )
            {
                // A footprint is the minumum info required to make a building.
                building = createSampleBuilding( feature );

                if (building.has_value())
                {
                    // Install the reference frame of the footprint geometry:
                    building->setReferenceFrame( local2world );

                    // Do initial cleaning of the footprint and install is:
                    cleanPolygon( polygon );

                    // Finally, build the internal structure from the footprint.
                    building->build( polygon, context, progress );
                }
            }
            else
            {
                //OE_WARN << LC << "Feature " << feature->getFID() << " is not a polygon. Skipping..\n";
            }
        }
    }

    return building;
}

void
BuildingFactory::cleanPolygon(Polygon* polygon)
{
    polygon->open();

    polygon->removeDuplicates();

    polygon->rewind( Polygon::ORIENTATION_CCW );

    // TODO: remove colinear points? for skeleton?
}

std::optional<Building>
BuildingFactory::createSampleBuilding(const Feature* feature)
{
    Building building;
    building.setUID( feature->getFID() );

    // figure out the building's height and number of floors.
    // single-elevation building.
    float height       = 15.0f;
    unsigned numFloors = 1u;

    // Add a single elevation.
    Elevation elevation;

    Roof roof;
    roof.setType( Roof::TYPE_FLAT );

    SkinResource* wallSkin = nullptr;
    SkinResource* roofSkin = nullptr;

    if ( _session.valid() )
    {
        ResourceLibrary* reslib = _session->styles()->getDefaultResourceLibrary();
        if ( reslib )
        {
            wallSkin = reslib->getSkin( "facade.commercial.1" );
            elevation.setSkinResource( wallSkin );

            roofSkin = reslib->getSkin( "roof.commercial.1" );
            roof.setSkinResource( roofSkin );
        }
        else
        {
            //OE_WARN << LC << "No resource library\n";
        }

        const BuildingSymbol* sym = _session->styles()->getDefaultStyle()->get<BuildingSymbol>();
        if ( sym )
        {
            FilterContext cx(_session.get());
            if ( feature )
            {
                height = sym->height()->eval(const_cast<Feature*>(feature), cx);
            }

            // calculate the number of floors
            if ( wallSkin )
            {
                numFloors = (unsigned)std::max(1.0f, osg::round(height / wallSkin->imageHeight().get()));
            }
            else
            {
                numFloors = (unsigned)std::max(1.0f, osg::round(height / sym->floorHeight().get()));
            }
        }
    }

    elevation.setHeight( height );
    elevation.setNumFloors( numFloors );
    elevation.setRoof( std::move(roof) );

    // Create parapet as a child elevation
    Elevation parapet;
    parapet.setParapet(true);
    parapet.setParapetWidth( 2.0f );
    parapet.setHeight( 2.0f );
    parapet.setColor( Color::Gray.brightness(1.3f) );

    Roof parapetRoof;
    parapetRoof.setSkinResource( roofSkin );
    parapetRoof.setColor( Color::Gray.brightness(1.2f) );
    parapet.setRoof( std::move(parapetRoof) );

    elevation.getElevations().emplace_back( std::move(parapet) );
    // Fix up parent pointer after push_back
    elevation.getElevations().back().setParent( &elevation );

    building.getElevations().emplace_back( std::move(elevation) );

    return building;
}
