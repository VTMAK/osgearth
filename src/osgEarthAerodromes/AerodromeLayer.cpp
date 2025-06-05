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
#include "AerodromeLayer"
#include "AerodromeFactory"
#include "AerodromeCatalog"
#include <osgEarth/NodeUtils>

using namespace osgEarth;
using namespace osgEarth::Aerodrome;

#define LC "[AerodromeLayer] "

REGISTER_OSGEARTH_LAYER(aerodrome, AerodromeLayer);
REGISTER_OSGEARTH_LAYER(aerodromes, AerodromeLayer);

void
AerodromeLayer::init()
{
    VisibleLayer::init();

    // Create the root group
    _root = new osg::Group();
    _root->setName(getName());
}

osg::Node*
AerodromeLayer::getNode() const
{
    return _root.get();
}

Status
AerodromeLayer::openImplementation()
{
    if (options().catalog().isSet())
    {
        _catalog = AerodromeCatalog::read(options().catalog().get(), getReadOptions());
    }

    if (!_catalog.valid())
    {
        _catalog = new AerodromeCatalog();
        _catalog->fromConfig(options().getConfig());
    }

    return VisibleLayer::openImplementation();
}

Status
AerodromeLayer::closeImplementation()
{
    _catalog = nullptr;
    return VisibleLayer::closeImplementation();
}

void
AerodromeLayer::addedToMap(const Map* map)
{
    // Hang on to the Map reference
    _map = map;

    // Recreate the scene graph
    createSceneGraph();
}

void
AerodromeLayer::removedFromMap(const Map* map)
{
    destroySceneGraph();
}

Layer::Stats
AerodromeLayer::reportStats() const
{
    Layer::Stats report;
    auto factory = osgEarth::findTopMostNodeOfType<AerodromeFactory>(_root.get());
    if (factory)
    {
        report.push_back({ "Resident aerodromes", std::to_string(*factory->_residentTiles) });
    }
    return report;
}

void
AerodromeLayer::createSceneGraph()
{
    // notify of removal:
    if (_root->getNumChildren() > 0)
        getSceneGraphCallbacks()->fireRemoveNode(_root->getChild(0));

    _root->removeChildren(0, _root->getNumChildren());

    if (options().renderOrder().isSet())
    {
        osg::StateSet* mss = _root->getOrCreateStateSet();
        mss->setRenderBinDetails(
            options().renderOrder().value(),
            mss->getBinName().empty() ? "DepthSortedBin" : mss->getBinName());
    }

    if (options().renderBin().isSet())
    {
        osg::StateSet* mss = _root->getOrCreateStateSet();
        mss->setRenderBinDetails(
            mss->getBinNumber(),
            options().renderBin().get());
    }

    osg::ref_ptr<const Map> map;
    if (_map.lock(map) && _catalog.valid())
    {
        osg::ref_ptr<AerodromeFactory> factory = new AerodromeFactory(
            map.get(),
            _catalog.get(),
            getSceneGraphCallbacks());

        if (options().range().isSet())
            factory->setRange(options().range().get());

        if (options().buildTerminals().isSet())
            factory->setBuildTerminals(options().buildTerminals().get());

        for (auto& icao : options().icaoExclusions())
            factory->icaoExclusions().insert(toLower(icao));

        if (options().clampToTerrain().isSet())
            factory->setClampToTerrain(options().clampToTerrain().get());

        factory->_serialLoading = options().serialLoading().get();

        factory->init(getReadOptions());

        _root->addChild(factory.get());

        // notify of addition:
        getSceneGraphCallbacks()->firePostMergeNode(factory.get());
    }
}

void
AerodromeLayer::destroySceneGraph()
{
    if (_root.valid())
    {
        _root->releaseGLObjects(nullptr);

        // notify of removal:
        if (_root->getNumChildren() > 0)
            getSceneGraphCallbacks()->fireRemoveNode(_root->getChild(0));

        _root->removeChildren(0, _root->getNumChildren());
    }
}