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
#include "Building"
#include "BuildingVisitor"
#include "BuildContext"
#include <osgEarth/Progress>
#include <sstream>

#define LC "[Building] "

using namespace osgEarth;
using namespace osgEarth::Buildings;

Building::Building() :
_zoning    ( Zoning::ZONING_UNKNOWN ),
_minHeight ( 0.0f ),
_maxHeight ( FLT_MAX ),
_minArea   ( 0.0f ),
_maxArea   ( FLT_MAX ),
_instanced ( false ),
_uid(0)
{
    //nop
}

Building::Building(const Building& rhs) :
_tags      ( rhs._tags ),
_uid       ( rhs._uid ),
_zoning    ( rhs._zoning ),
_local2world( rhs._local2world ),
_minHeight ( rhs._minHeight ),
_maxHeight ( rhs._maxHeight ),
_minArea   ( rhs._minArea ),
_maxArea   ( rhs._maxArea ),
_instanced ( rhs._instanced ),
_externalModelURI      ( rhs._externalModelURI ),
_instancedModelSymbol  ( rhs._instancedModelSymbol),
_instancedModelResource( rhs._instancedModelResource)
{
    for (const auto& e : rhs._elevations)
        _elevations.emplace_back(e);
}

Building::Building(Building&& rhs) noexcept :
_tags      ( std::move(rhs._tags) ),
_uid       ( rhs._uid ),
_zoning    ( rhs._zoning ),
_local2world( rhs._local2world ),
_minHeight ( rhs._minHeight ),
_maxHeight ( rhs._maxHeight ),
_minArea   ( rhs._minArea ),
_maxArea   ( rhs._maxArea ),
_instanced ( rhs._instanced ),
_externalModelURI      ( std::move(rhs._externalModelURI) ),
_instancedModelSymbol  ( std::move(rhs._instancedModelSymbol) ),
_instancedModelResource( std::move(rhs._instancedModelResource) ),
_elevations( std::move(rhs._elevations) )
{
}

Building& Building::operator=(const Building& rhs)
{
    if (this != &rhs)
    {
        _tags = rhs._tags;
        _uid = rhs._uid;
        _zoning = rhs._zoning;
        _local2world = rhs._local2world;
        _minHeight = rhs._minHeight;
        _maxHeight = rhs._maxHeight;
        _minArea = rhs._minArea;
        _maxArea = rhs._maxArea;
        _instanced = rhs._instanced;
        _externalModelURI = rhs._externalModelURI;
        _instancedModelSymbol = rhs._instancedModelSymbol;
        _instancedModelResource = rhs._instancedModelResource;
        _elevations.clear();
        for (const auto& e : rhs._elevations)
            _elevations.push_back(e);
    }
    return *this;
}

Building& Building::operator=(Building&& rhs) noexcept
{
    if (this != &rhs)
    {
        _tags = std::move(rhs._tags);
        _uid = rhs._uid;
        _zoning = rhs._zoning;
        _local2world = rhs._local2world;
        _minHeight = rhs._minHeight;
        _maxHeight = rhs._maxHeight;
        _minArea = rhs._minArea;
        _maxArea = rhs._maxArea;
        _instanced = rhs._instanced;
        _externalModelURI = std::move(rhs._externalModelURI);
        _instancedModelSymbol = std::move(rhs._instancedModelSymbol);
        _instancedModelResource = std::move(rhs._instancedModelResource);
        _elevations = std::move(rhs._elevations);
    }
    return *this;
}

Building Building::clone() const
{
    return Building(*this);
}

void Building::addTag(const std::string& tag)
{
    _tags.insert(normalize(tag));
}

void Building::addTags(const TagVector& tags)
{
    for (const auto& tag : tags)
        _tags.insert(normalize(tag));
}

void Building::addTags(const std::string& tagString)
{
    auto tags = StringTokenizer()
        .delim(" ")
        .standardQuotes()
        .keepEmpties(false)
        .tokenize(tagString);
    addTags(tags);
}

void Building::removeTag(const std::string& tag)
{
    _tags.erase(normalize(tag));
}

bool Building::containsTag(const std::string& tag) const
{
    return _tags.find(normalize(tag)) != _tags.end();
}

bool Building::containsTags(const TagSet& tags) const
{
    for (const auto& tag : tags)
    {
        if (_tags.find(normalize(tag)) == _tags.end())
            return false;
    }
    return true;
}

bool Building::containsTags(const TagVector& tags) const
{
    for (const auto& tag : tags)
    {
        if (_tags.find(normalize(tag)) == _tags.end())
            return false;
    }
    return true;
}

std::string Building::tagString() const
{
    std::stringstream buf;
    for (auto i = _tags.begin(); i != _tags.end(); ++i)
        buf << (i != _tags.begin() ? " " : "") << *i;
    return buf.str();
}

std::string Building::tagString(const TagSet& tags)
{
    std::stringstream buf;
    for (auto i = tags.begin(); i != tags.end(); ++i)
        buf << (i != tags.begin() ? " " : "") << *i;
    return buf.str();
}

std::string Building::tagString(const TagVector& tags)
{
    std::stringstream buf;
    for (auto i = tags.begin(); i != tags.end(); ++i)
        buf << (i != tags.begin() ? " " : "") << *i;
    return buf.str();
}

std::string Building::normalize(const std::string& input) const
{
    return osgEarth::toLower(input);
}

void
Building::setHeight(float height)
{
    for (auto& e : _elevations)
    {
        e.setHeight(height);
    }
}

bool
Building::build(const Polygon* footprint, BuildContext& bc, ProgressCallback* progress)
{
    if ( !footprint || !footprint->isValid() )
        return false;

    // Resolve an instanced building model if available.
    resolveInstancedModel(bc, progress);

    // In the absence of an instanced model, build parametric data.
    if ( getInstancedModelResource() == nullptr )
    {
        for (auto& e : _elevations)
        {
            e.build(footprint, bc);
        }
    }

    // if we are using an instanced model, we still need the rotation/AABB from the first elevation:
    else if ( !_elevations.empty() )
    {
        _elevations.front().calculateRotations(footprint);
    }

    return true;
}

void
Building::resolveInstancedModel(BuildContext& bc, ProgressCallback* progress)
{
    if ( getInstancedModelSymbol() && bc.getResourceLibrary() )
    {
        // resolve the resource.
        ModelResourceVector candidates;
        bc.getResourceLibrary()->getModels( getInstancedModelSymbol(), candidates, bc.getDBOptions() );
        if ( !candidates.empty() )
        {
            unsigned index = Random(bc.getSeed()).next( candidates.size() );
            setInstancedModelResource( candidates.at(index).get() );
        }
        else
        {
            if (progress)
            {
                progress->message() += Stringify()
                    << "Failed to instantiate a model for \""
                    << getInstancedModelSymbol()->getConfig().toJSON(true) << "\"\n";
            }

            OE_WARN << LC << "Failed to instantiate a model for \""
                << getInstancedModelSymbol()->getConfig().toJSON(true) << "\"" << std::endl;
        }
    }
}

void
Building::accept(BuildingVisitor& bv)
{
    bv.apply( this );
}

Config
Building::getConfig() const
{
    //TODO: incomplete.
    Config conf("building");
    conf.set("external_model_uri", _externalModelURI);
    if ( !getElevations().empty() )
    {
        Config evec("elevations");
        for (const auto& e : getElevations())
            evec.add("elevation", e.getConfig());
        conf.set(evec);
    }
    return conf;
}
