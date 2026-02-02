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
#include "BuildingCompiler"

#include <osg/Geometry>
#include <osg/Geode>
#include <osg/MatrixTransform>

#include "ExternalModelNode"
#include <osgEarth/OEAssert>

#define LC "[BuildingCompiler] "

using namespace osgEarth;
using namespace osgEarth::Buildings;


BuildingCompiler::BuildingCompiler(Session* session) :
  _session( session )
  , _filterUsage(FILTER_USAGE_NORMAL),
    _clutter(true),
    _parapets(true)
{
    _elevationCompiler = new ElevationCompiler( session );
    _flatRoofCompiler = new FlatRoofCompiler( session );
    _gableRoofCompiler = new GableRoofCompiler( session );
    _instancedRoofCompiler = new InstancedRoofCompiler( session );
    _instancedBuildingCompiler = new InstancedBuildingCompiler( session );
}

void BuildingCompiler::setUsage(FilterUsage usage)
{
   _filterUsage = usage;
   _flatRoofCompiler->setUsage(_filterUsage);
}

void BuildingCompiler::setClutter(bool value)
{
    _clutter = value;
}

void BuildingCompiler::setParapets(bool value)
{
    _parapets = value;
}

bool
BuildingCompiler::compile(const BuildingVector& input,
                          CompilerOutput&       output,
                          const osgDB::Options* readOptions,
                          ProgressCallback*     progress)
{
    OE_START_TIMER(total);

    for(auto& building : input)
    {
        if ( progress && progress->isCanceled() )
        {
            progress->message() = "in BuildingCompiler::compile()";
            return false;
        }

        if ( building.externalModelURI().isSet() )
        {
            addExternalModel( output, building, output.getWorldToLocal(), readOptions, progress );
        }
        else if ( building.getInstancedModelResource() )
        {
           osg::Matrixd _worldToLocal = (_filterUsage == FILTER_USAGE_NORMAL) ? output.getWorldToLocal() : osg::Matrixd::identity();
           _instancedBuildingCompiler->compile(building, output, _worldToLocal, progress);
        }
        else
        {
           addElevations( output, building, building.getElevations(), output.getWorldToLocal(), readOptions);
        }
    }

    return true;
}

bool
BuildingCompiler::addExternalModel(CompilerOutput&       output,
                                   const Building&       building,
                                   const osg::Matrix&    world2local,
                                   const osgDB::Options* readOptions,
                                   ProgressCallback*     progress) const
{
   if (_filterUsage == FILTER_USAGE_NORMAL)
   {
      // TODO: perhaps disable the image caching for an external model, 
      //       since it probably won't be shared?
      osg::ref_ptr<osg::Node> node = building.getExternalModelURI().getNode(readOptions, progress);
      if (node.valid())
      {
         osg::MatrixTransform* xform = new osg::MatrixTransform(building.getReferenceFrame() * world2local);
         xform->addChild(node.get());

         MetadataNode* metadataNode = output.getMetadata();
         Feature* feature = output.getCurrentFeature();
         if(metadataNode && feature)
         {
            auto id = metadataNode->add(feature, true);
            metadataNode->tagNode(xform, id);
         }

         output.getExternalModelsGroup()->addChild(xform);
      }
      return node.valid();
   }
   else
   {
      if (output.getExternalModelsGroup()->getNumChildren() == 0)
      {
         osg::Group* group = new osg::Group;
         group->setName("ExternalModelNodeListAttachPoint");

         ExternalModelNode* externalModelNodeList = new ExternalModelNode();
         group->addChild(externalModelNodeList);

         output.getExternalModelsGroup()->addChild(group);
      }

      ASSERT_PREDICATE(output.getExternalModelsGroup()->getNumChildren() == 1);
      ExternalModelNode* externalModelNode =
         static_cast<ExternalModelNode*>(output.getExternalModelsGroup()->getChild(0)->asGroup()->getChild(0));
      ExternalModel * externalModel = new ExternalModel();

      externalModel->modelName = building.getExternalModelURI().full();
      externalModel->xform = building.getReferenceFrame()/* * world2local*/;

      MetadataNode* metadataNode = output.getMetadata();
      Feature* feature = output.getCurrentFeature();
      if (metadataNode && feature)
      {
         externalModel->objectId = metadataNode->add(feature, true);
      }

      externalModelNode->vectorExternalModels.emplace_back(externalModel);
      return true;
   }
}

bool
BuildingCompiler::addElevations(CompilerOutput&        output,
                                const Building&        building,
                                const ElevationVector& elevations,
                                const osg::Matrix&     world2local,
                                const osgDB::Options*  readOptions) const
{
   if (!building) return false;


   // Iterator over each Elevation in this building:
   for (const auto& elevation : elevations)
   {
       // Skip parapets if they are disabled
       if (elevation.isParapet() && !_parapets)
       {
           continue;
       }

      _elevationCompiler->compile(output, building, elevation, world2local, readOptions);

      if (elevation.getRoof())
      {
         addRoof(output, building, elevation, world2local, readOptions);
      }

      if (!elevation.getElevations().empty())
      {
         addElevations(output, building, elevation.getElevations(), world2local, readOptions);
      }

      // if this is an instanced model but without a model resource, we are
      // falling back on the parametric template. In this case, if there is no
      // roof defined, add a flat one by default.
      if (building.getInstanced() &&
          building.getInstancedModelResource() == nullptr &&
          elevation.getRoof() == nullptr)
      {
          addRoof(output, building, elevation, world2local, readOptions);
      }

   } // elevations loop


   return true;
}

bool
BuildingCompiler::addRoof(
    CompilerOutput&       output, 
    const Building&       building, 
    const Elevation& elevation,
    const osg::Matrix& world2local,
    const osgDB::Options* readOptions) const
{
    auto* roof = elevation.getRoof();

    if (roof && roof->getType() == Roof::TYPE_GABLE)
    {
        return _gableRoofCompiler->compile(output, building, elevation, world2local, readOptions);
    }
    else if (roof && roof->getType() == Roof::TYPE_INSTANCED)
    {
        return _instancedRoofCompiler->compile(output, building, elevation, world2local, readOptions);
    }
    else // flat roof.
    {
        return _flatRoofCompiler->compile(output, building, elevation, _clutter, world2local, readOptions);
    }

    return false;
}
