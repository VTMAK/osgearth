/* osgEarth
 * Copyright 2026 Pelican Mapping
 * MIT License
 */
#include <osgEarthPrestige/SelectFilter>
#include <osgEarth/StringUtils>
#include <algorithm>

using namespace osgEarth;
using namespace osgEarthPrestige;

OSGEARTH_REGISTER_SIMPLE_FEATUREFILTER(prestige_select, SelectFilter);

SelectFilter::SelectFilter(const Config& conf)
{
    _attribute = conf.value("attribute");
    for (auto& value : StringTokenizer().delim(",").tokenize(conf.value("values")))
    {
        value = trim(value);
        if (!value.empty()) _values.push_back(value);
    }
}

Status SelectFilter::initialize(const osgDB::Options*)
{
    if (_attribute.empty() || _values.empty())
        return Status(Status::ConfigurationError, "prestige_select requires an attribute and at least one value");
    return Status::NoError;
}

Config SelectFilter::getConfig() const
{
    Config conf("prestige_select");
    conf.set("attribute", _attribute);
    conf.set("values", joinStrings(_values, ','));
    return conf;
}

FilterContext SelectFilter::push(FeatureList& input, FilterContext& context)
{
    input.erase(std::remove_if(input.begin(), input.end(), [&](const osg::ref_ptr<Feature>& feature)
    {
        if (!feature.valid() || !feature->hasAttr(_attribute)) return true;
        const std::string value = feature->getString(_attribute);
        return std::find(_values.begin(), _values.end(), value) == _values.end();
    }), input.end());
    return context;
}
