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
#include "ShaderLoader"
#include "ShaderUtils"
#include "URI"
#include "VirtualProgram"

#include <osgDB/FileUtils>

#undef  LC
#define LC "[ShaderLoader] "

using namespace osgEarth;
using namespace osgEarth::Util;


namespace
{
    //TODO: unordered?
    using StringMap = std::map<std::string, std::string>;


    bool parseLocation(
        const std::string& loc,
        optional<ShaderComp::FunctionLocation>& location)
    {
        bool locationSet = true;

        if (ciEquals(loc, "vertex_model"))
            location = ShaderComp::LOCATION_VERTEX_MODEL;
        else if (ciEquals(loc, "vertex_view"))
            location = ShaderComp::LOCATION_VERTEX_VIEW;
        else if (ciEquals(loc, "vertex_clip"))
            location = ShaderComp::LOCATION_VERTEX_CLIP;
        else if (ciEquals(loc, "tess_control") || ciEquals(loc, "tessellation_control"))
            location = ShaderComp::LOCATION_TESS_CONTROL;
        else if (ciEquals(loc, "tess_eval") || ciEquals(loc, "tessellation_eval") || ciEquals(loc, "tessellation_evaluation") || ciEquals(loc, "tess_evaluation"))
            location = ShaderComp::LOCATION_TESS_EVALUATION;
        else if (ciEquals(loc, "vertex_geometry") || ciEquals(loc, "geometry"))
            location = ShaderComp::LOCATION_GEOMETRY;
        else if (ciEquals(loc, "fragment"))
            location = ShaderComp::LOCATION_FRAGMENT_COLORING;
        else if (ciEquals(loc, "fragment_coloring"))
            location = ShaderComp::LOCATION_FRAGMENT_COLORING;
        else if (ciEquals(loc, "fragment_lighting"))
            location = ShaderComp::LOCATION_FRAGMENT_LIGHTING;
        else if (ciEquals(loc, "fragment_output"))
            location = ShaderComp::LOCATION_FRAGMENT_OUTPUT;
        else
            locationSet = false;

        return locationSet;
    }


    struct VPFunction
    {
        VPFunction()
        {
            order.setDefault(1.0f);
        }

        std::string entryPoint;
        optional<ShaderComp::FunctionLocation> location;
        optional<float> order;
    };

    void getVPFunction(const std::string& source, VPFunction& f)
    {
        std::string::size_type pragmaPos = source.find("#pragma vp_function");
        if (pragmaPos != std::string::npos)
        {
            std::string line = ShaderLoader::getPragmaValue(source, "vp_function");
            if (!line.empty())
            {
                StringVector tokens;
                StringTokenizer(line, tokens, " ,\t", "", false, true);
                if (tokens.size() > 0)
                    f.entryPoint = tokens[0];
                if (tokens.size() > 1)
                    parseLocation(tokens[1], f.location);
                if (tokens.size() > 2)
                    f.order = atof(tokens[2].c_str());
            }
        }

        std::string entryPointStr = ShaderLoader::getPragmaValue(source, "vp_entryPoint");
        if (!entryPointStr.empty())
        {
            f.entryPoint = entryPointStr;
        }

        std::string locationStr = ShaderLoader::getPragmaValue(source, "vp_location");
        if (!locationStr.empty())
        {
            parseLocation(locationStr, f.location);
        }

        std::string orderStr = ShaderLoader::getPragmaValue(source, "vp_order");
        if (!orderStr.empty())
        {
            if (ciEquals(orderStr, "FLT_MAX") || ciEquals(orderStr, "last"))
                f.order = FLT_MAX;
            else if (ciEquals(orderStr, "-FLT_MAX") || ciEquals(orderStr, "first"))
                f.order = -FLT_MAX;
            else
                f.order = as<float>(orderStr, 1.0f);
        }
    }

    std::vector<std::string>
    splitAtEndOfLineStartingWith(const std::string& source, const std::string& token)
    {
        std::vector<std::string> result(2);

        std::string::size_type tokenPos = source.find(token);
        if (tokenPos == std::string::npos)
        {
            result[1] = source;
            return result;
        }

        std::string::size_type newlinePos = source.find('\n', tokenPos);
        if (newlinePos == std::string::npos)
        {
            result[0] = source;
            return result;
        }

        result[0] = source.substr(0, newlinePos+1);
        result[1] = source.substr(newlinePos+1);
        return result;
    }

    void insertStageDefine(std::string& source, osg::Shader::Type stage)
    {
        std::string define;
        if (stage == osg::Shader::VERTEX) define = "VP_STAGE_VERTEX";
        else if (stage == osg::Shader::TESSCONTROL) define = "VP_STAGE_TESSCONTROL";
        else if (stage == osg::Shader::TESSEVALUATION) define = "VP_STAGE_TESSEVALUATION";
        else if (stage == osg::Shader::GEOMETRY) define = "VP_STAGE_GEOMERTY";
        else if (stage == osg::Shader::COMPUTE) define = "VP_STAGE_COMPUTE";
        else if (stage == osg::Shader::FRAGMENT) define = "VP_STAGE_FRAGMENT";
        else define = "UNDEFINED";

        std::vector<std::string> parts = splitAtEndOfLineStartingWith(source, "#version");
        if (parts.size() == 2)
        {
            source =
                parts[0] +
                "#define " + define + "\n" +
                parts[1];
        }
    }

    osg::Shader::Type getShaderTypeFromLocation(ShaderComp::FunctionLocation loc)
    {
        // If a location is set, install in that location only
        if (loc == ShaderComp::LOCATION_VERTEX_MODEL ||
            loc == ShaderComp::LOCATION_VERTEX_VIEW ||
            loc == ShaderComp::LOCATION_VERTEX_CLIP)
        {
            return osg::Shader::VERTEX;
        }
        else if (
            loc == ShaderComp::LOCATION_FRAGMENT_COLORING ||
            loc == ShaderComp::LOCATION_FRAGMENT_LIGHTING ||
            loc == ShaderComp::LOCATION_FRAGMENT_OUTPUT)
        {
            return osg::Shader::FRAGMENT;
        }
        else if (
            loc == ShaderComp::LOCATION_GEOMETRY)
        {
            return osg::Shader::GEOMETRY;
        }
        else if (
            loc == ShaderComp::LOCATION_TESS_CONTROL)
        {
            return osg::Shader::TESSCONTROL;
        }
        else if (
            loc == ShaderComp::LOCATION_TESS_EVALUATION)
        {
            return osg::Shader::TESSEVALUATION;
        }
        else
        {
            return osg::Shader::UNDEFINED;
        }
    }
}


// find the value of a pragma, e.g.:
//   #pragma oe_key value
// returns the string "value" (without the quotes).
std::string
ShaderLoader::getPragmaValue(const std::string& source, const std::string& key)
{
    std::string token("#pragma " + key);
    std::string::size_type statementPos = source.find(token);
    if ( statementPos == std::string::npos )
        return "";

    // no quotes; parse to newline.
    std::string::size_type startPos = source.find_first_not_of(" \t", statementPos+token.length());
    if ( startPos == std::string::npos )
        return ""; // no whitespace after the pragma key

    std::string::size_type newlinePos = source.find('\n', startPos);
    if ( newlinePos == std::string::npos )
        return ""; // new newline found after pragma

    return trim(source.substr(startPos, newlinePos-startPos));
}

bool
ShaderLoader::getPragmaValueAsTokens(
    const std::string& input,
    const std::string& key,
    std::string& line_out,
    std::vector<std::string>& tokens_out)
{
    std::string::size_type statementPos = input.find(key);
    if (statementPos == std::string::npos)
        return 0;

    std::string::size_type startPos = input.find_first_not_of(" \t(", statementPos + key.length());
    if (startPos == std::string::npos)
        return 0;

    std::string::size_type endPos = input.find_first_of(")\n", startPos);
    if (endPos == std::string::npos)
        return 0;

    std::string::size_type nlPos = input.find('\n', startPos);
    if (nlPos == std::string::npos)
        return 0;

    line_out = input.substr(statementPos, nlPos - statementPos);
    std::string statement(input.substr(statementPos, endPos - statementPos));
    std::string value(trim(input.substr(startPos, endPos - startPos)));

    StringTokenizer(value, tokens_out, ", \t", "", false, true);
    return tokens_out.size();
}

void
ShaderLoader::getAllPragmaValues(const std::string&     source,
                                 const std::string&     key,
                                 std::set<std::string>& output)
{
    std::string token("#pragma " + key);
    std::string::size_type pragmaPos = 0;
    while( pragmaPos != std::string::npos )
    {
        pragmaPos = source.find(token, pragmaPos);
        if ( pragmaPos != std::string::npos )
        {
            std::string::size_type startPos = source.find_first_not_of(" \t", pragmaPos+token.length());
            if ( startPos != std::string::npos )
            {
                std::string::size_type newlinePos = source.find('\n', startPos);
                if ( newlinePos != std::string::npos )
                {
                    const size_t len = newlinePos - startPos;
                    if ( len > 0 )
                    {
                        output.insert( trim(source.substr(startPos, len)) );
                    }
                    pragmaPos = newlinePos;
                }
                else
                {
                    break;
                }
            }
            else
            {
                break;
            }
        }
    }
}

bool
ShaderLoader::load(
    VirtualProgram* vp,
    const std::string& source)
{
    ShaderPackage pkg;
    pkg.add("", source);
    return load(vp, "", pkg, nullptr);


}

std::string
ShaderLoader::load(const std::string&    filename,
                   const ShaderPackage&  package,
                   const osgDB::Options* dbOptions)
{
    std::string output;
    bool useInlineSource = false;

    URIContext context( dbOptions );
    URI uri(filename, context);

    std::string inlineSource;
    ShaderPackage::SourceMap::const_iterator source = package._sources.find(filename);
    if ( source != package._sources.end() )
        inlineSource = source->second;

    if (!filename.empty())
    {
        std::string path = osgDB::findDataFile(uri.full(), dbOptions);
        if (path.empty())
        {
            output = inlineSource;
            useInlineSource = true;
            if (inlineSource.empty())
            {
                OE_WARN << LC << "Inline source for \"" << filename << "\" is empty, and no external file could be found.\n";
            }
        }
        else
        {
            std::string externalSource = URI(path, context).getString(dbOptions);
            if (!externalSource.empty())
            {
                OE_DEBUG << LC << "Loaded external shader " << filename << " from " << path << "\n";
                output = externalSource;
            }
            else
            {
                output = inlineSource;
                useInlineSource = true;
            }
        }
    }
    else
    {
        output = inlineSource;
        useInlineSource = true;
    }

    // replace common tokens:
    osgEarth::replaceIn(output, "$GLSL_VERSION_STR", GLSL_VERSION_STR);
    osgEarth::replaceIn(output, "$GLSL_DEFAULT_PRECISION_FLOAT", GLSL_DEFAULT_PRECISION_FLOAT);

    // If we're using inline source, we have to post-process the string.
    if ( useInlineSource )
    {
        // Replace tokens inserted in the CMakeModules/ConfigureShaders.cmake.in script.
        osgEarth::replaceIn(output, "%EOL%",   "\n");
        osgEarth::replaceIn(output, "%QUOTE%", "\"");
    }

    // Process any user-defined replacements.
    for (ShaderPackage::ReplaceMap::const_iterator i = package._replaces.begin();
        i != package._replaces.end();
        ++i)
    {
        osgEarth::replaceIn(output, i->first, i->second);
    }

    // Run the "pre" callbacks before includes
    ShaderPreProcessor::runPre(output);

    // Process any "#pragma include" statements
    while(true)
    {
        const std::string token("#pragma include");
        std::string::size_type statementPos = output.find(token);
        if ( statementPos == std::string::npos )
            break;

        std::string::size_type startPos = output.find_first_not_of(" \t", statementPos+token.length());
        if ( startPos == std::string::npos )
            break;

        std::string::size_type endPos = output.find('\n', startPos);
        if ( endPos == std::string::npos )
            break;

        std::string statement( output.substr(statementPos, endPos-statementPos) );
        std::string fileToInclude( trim(output.substr(startPos, endPos-startPos)) );

        // load the source of the included file, and append a newline so we
        // don't break the MULTILINE macro if the last line of the include
        // file is a comment.
        std::string fileSource = Stringify()
            << load(fileToInclude, package, dbOptions)
            << "\n";

        Strings::replaceIn(output, statement, fileSource);
    }

    // Process any "#pragma define" statements
    while (true)
    {
        const std::string token("#pragma vp_define");
        std::string::size_type statementPos = output.find(token);
        if (statementPos == std::string::npos)
            break;

        std::string::size_type startPos = output.find_first_not_of(" \t", statementPos + token.length());
        if (startPos == std::string::npos)
            break;

        std::string::size_type endPos = output.find('\n', startPos);
        if (endPos == std::string::npos)
            break;

        std::string statement(output.substr(statementPos, endPos - statementPos));
        std::string varName(trim(output.substr(startPos, endPos - startPos)));

        ShaderPackage::DefineMap::const_iterator d = package._defines.find(varName);

        bool defineIt =
            d != package._defines.end() &&
            d->second == true;

        std::string newStatement = Stringify()
            << (defineIt ? "#define " : "#undef ")
            << varName;

        Strings::replaceIn(output, statement, newStatement);
    }

    // Lastly, remove any CRs
    osgEarth::replaceIn(output, "\r", "");

    return output;
}

std::string
ShaderLoader::load(const std::string&    filename,
                   const std::string&    inlineSource,
                   const osgDB::Options* dbOptions)
{
    std::string output;
    bool useInlineSource = false;

    URIContext context(dbOptions);
    URI uri(filename, context);

    std::string path = osgDB::findDataFile(filename, dbOptions);
    if (path.empty())
    {
        output = inlineSource;
        useInlineSource = true;
    }
    else
    {
        std::string externalSource = URI(path, context).getString(dbOptions);
        if (!externalSource.empty())
        {
            OE_DEBUG << LC << "Loaded external shader " << filename << " from " << path << "\n";
            output = externalSource;
        }
        else
        {
            output = inlineSource;
            useInlineSource = true;
        }
    }

    // replace common tokens:
    osgEarth::replaceIn(output, "$GLSL_VERSION_STR", GLSL_VERSION_STR);
    osgEarth::replaceIn(output, "$GLSL_DEFAULT_PRECISION_FLOAT", GLSL_DEFAULT_PRECISION_FLOAT);

    // If we're using inline source, we have to post-process the string.
    if (useInlineSource)
    {
        // Replace tokens inserted in the CMakeModules/ConfigureShaders.cmake.in script.
        osgEarth::replaceIn(output, "%EOL%", "\n");
        osgEarth::replaceIn(output, "%QUOTE%", "\"");
    }

    // Lastly, remove any CRs
    osgEarth::replaceIn(output, "\r", "");

    return output;
}

void
ShaderLoader::split(const std::string& multisource,
                    std::vector<std::string>& output)
{
#define SPLIT_DELIM "[break]"
#define SPLIT_DELIM_LEN 7
    std::string::size_type offset = 0, pos = 0;
    while ((pos = multisource.find(SPLIT_DELIM, offset)) != std::string::npos)
    {
        std::string source = multisource.substr(offset, pos-offset);
        output.push_back(source);
        offset = pos + SPLIT_DELIM_LEN;
    }
    output.push_back(multisource.substr(offset));
}

bool
ShaderLoader::load(VirtualProgram*       vp,
                   const std::string&    filename,
                   const ShaderPackage&  package,
                   const osgDB::Options* dbOptions)
{
    if ( !vp )
    {
        OE_WARN << LC << "Illegal: null VirtualProgram\n";
        return false;
    }

    // load the source string:
    std::string multisource = load(filename, package, dbOptions);
    if ( multisource.empty() )
    {
        OE_WARN << LC << "Failed to load shader source from \"" << filename << "\"\n";
        return false;
    }

    // split the multisource string into one or more shader sources:
    std::vector<std::string> sources;
    split(multisource, sources);

    for (unsigned i = 0; i < sources.size(); ++i)
    {
        std::string source = sources[i];

        // Remove the quotation marks from the source since they are illegal in GLSL
        replaceIn( source, "\"", " ");

        // Named?
        if (vp->getName().empty())
        {
            std::string name = getPragmaValue(source, "vp_name");
            vp->setName(name);
        }

        VPFunction f;
        getVPFunction(source, f);

        if (!f.entryPoint.empty())
        {
            if (f.location.isSet() == false)
                f.location = ShaderComp::LOCATION_FRAGMENT_COLORING;

            insertStageDefine(source, getShaderTypeFromLocation(f.location.get()));

            // set the function!
            vp->setFunction(
                f.entryPoint, 
                source, 
                f.location.get(),
                f.order.get());
        }

        else // no entry point - library shader.
        {
            // install as a simple shader.
            if (f.location.isSet())
            {
                // If a location is set, install in that location only
                osg::Shader::Type type = getShaderTypeFromLocation(f.location.get());

                insertStageDefine(source, type);

                osg::Shader* shader = new osg::Shader(type, source);
                shader->setName( filename );
                vp->setShader( filename, shader );
            }

            else
            {
                // If no location was set, install in all stages.
                const osg::Shader::Type types[5] = { 
                    osg::Shader::VERTEX, 
                    osg::Shader::FRAGMENT, 
                    osg::Shader::GEOMETRY, 
                    osg::Shader::TESSCONTROL, 
                    osg::Shader::TESSEVALUATION
                };

                for(int i=0; i<5; ++i)
                {
                    std::string new_source = source;
                    insertStageDefine(new_source, types[i]);
                    osg::Shader* shader = new osg::Shader(types[i], new_source);
                    std::string name = Stringify() << filename + "_" + shader->getTypename();
                    shader->setName( name );
                    vp->setShader( name, shader );
                }
            }
        }
    }

    return true;
}

bool
ShaderLoader::unload(
    VirtualProgram* vp,
    const std::string& source)
{
    ShaderPackage pkg;
    pkg.add("", source);
    return unload(vp, "", pkg, nullptr);
}

bool
ShaderLoader::unload(VirtualProgram*       vp,
                     const std::string&    filename,
                     const ShaderPackage&  package,
                     const osgDB::Options* dbOptions)
{
    if ( !vp )
    {
        // fail quietly
        return false;
    }
    
    // load the source string:
    std::string multisource = load(filename, package, dbOptions);
    if ( multisource.empty() )
    {
        OE_WARN << LC << "Failed to load shader source from \"" << filename << "\"\n";
        return false;
    }

    // split the multisource string into one or more shader sources:
    std::vector<std::string> sources;
    split(multisource, sources);

    for (unsigned i = 0; i < sources.size(); ++i)
    {
        const std::string& source = sources[i];

        VPFunction f;
        getVPFunction(source, f);

        if ( !f.entryPoint.empty() )
        {
            vp->removeShader( f.entryPoint );
        }
        else
        {
            vp->removeShader( filename );
        }
    }

    return true;
}

//...................................................................

void
ShaderPackage::define(const std::string& name,
                      bool               defOrUndef)
{
    _defines[name] = defOrUndef;
}

void
ShaderPackage::replace(const std::string& pattern,
                       const std::string& value)
{
    _replaces[pattern] = value;
}

bool
ShaderPackage::load(VirtualProgram*       vp,
                    const std::string&    filename,
                    const osgDB::Options* dbOptions) const
{
    return ShaderLoader::load(vp, filename, *this, dbOptions);
}

bool
ShaderPackage::unload(VirtualProgram*       vp,
                      const std::string&    filename,
                      const osgDB::Options* dbOptions) const
{
    return ShaderLoader::unload(vp, filename, *this, dbOptions);
}

bool
ShaderPackage::loadAll(VirtualProgram*       vp,
                       const osgDB::Options* dbOptions) const
{
    int oks = 0;
    for(SourceMap::const_iterator i = _sources.begin(); i != _sources.end(); ++i)
    {
        oks += load( vp, i->first ) ? 1 : 0;
    }
    return oks == _sources.size();
}

bool
ShaderPackage::unloadAll(VirtualProgram*       vp,
                          const osgDB::Options* dbOptions) const
{
    int oks = 0;
    for(SourceMap::const_iterator i = _sources.begin(); i != _sources.end(); ++i)
    {
        oks += unload( vp, i->first ) ? 1 : 0;
    }
    return oks == _sources.size();
}
