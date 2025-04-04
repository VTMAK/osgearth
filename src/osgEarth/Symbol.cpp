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
#include <osgEarth/Symbol>

using namespace osgEarth;

SymbolRegistry*
SymbolRegistry::instance()
{
    static SymbolRegistry* s_singleton =0L;
    static std::mutex s_singletonMutex;

    if ( !s_singleton )
    {
        std::lock_guard<std::mutex> lock(s_singletonMutex);
        if ( !s_singleton )
        {
            s_singleton = new SymbolRegistry();
        }
    }
    return s_singleton;
}

void
SymbolRegistry::add( SymbolFactory* factory )
{
    _factories.push_back( factory );
}

Symbol*
SymbolRegistry::create( const Config& conf )
{
    for (SymbolFactoryList::iterator itr = _factories.begin(); itr != _factories.end(); itr++)
    {
        Symbol* symbol = itr->get()->create( conf );
        if (symbol) return symbol;
    }
    return 0;
} 

void SymbolRegistry::parseSLD(const Config& c, class Style& style) const
{
    for (SymbolFactoryList::const_iterator itr = _factories.begin(); itr != _factories.end(); itr++)
    {
        itr->get()->parseSLD( c, style );     
    }
}

//------------------------------------------------------------------------

Symbol::Symbol(const Config& conf)
{
    _uriContext = URIContext(conf.referrer());
    mergeConfig(conf);
}

Symbol::Symbol(const Symbol& rhs, const osg::CopyOp& copyop) :
    osg::Object(rhs, copyop)
{
    _library = rhs._library;
    _script = rhs._script;
    _uriContext = rhs._uriContext;
}

void
Symbol::mergeConfig(const Config& conf)
{
    conf.get("script", script());
    conf.get("library", library());
    if (conf.hasChild("__original"))
        _ctorConfig = conf.child("__original");
    else
        _ctorConfig = conf;
}

Config
Symbol::getConfig() const
{
    Config conf;
    conf.set("script", script());
    conf.set("library", library());
    if (!_ctorConfig.empty())
        conf.add("__original", _ctorConfig);
    return conf;
}

bool
Symbol::match(const std::string& s, const char* reservedWord)
{
    if ( s.compare(reservedWord) == 0 ) return true;
    std::string temp1 = toLower(s), temp2 = toLower(reservedWord);
    replaceIn(temp1, "_", "-");
    replaceIn(temp2, "_", "-");
    return temp1.compare(temp2) == 0;
}
