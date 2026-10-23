/* osgEarth
 * Copyright 2026 Pelican Mapping
 * MIT License
 */
#include <osgEarth/Registry>
#include <osgEarth/URI>
#include <osgDB/FileNameUtils>
#include <osgDB/ReaderWriter>
#include <osgDB/Registry>

namespace
{
    constexpr const char* FLIP_VERTICAL_EXTENSION = "flipvertical";

    class FlipVerticalPseudoLoader : public osgDB::ReaderWriter
    {
    public:
        FlipVerticalPseudoLoader()
        {
            supportsExtension(
                FLIP_VERTICAL_EXTENSION,
                "Vertical image flip pseudoloader");
        }

        const char* className() const override
        {
            return "osgEarth vertical image flip pseudoloader";
        }

        ReadResult readImage(
            const std::string& filename,
            const osgDB::Options* options) const override
        {
            if (!osgDB::equalCaseInsensitive(
                    osgDB::getFileExtension(filename),
                    FLIP_VERTICAL_EXTENSION))
            {
                return ReadResult::FILE_NOT_HANDLED;
            }

            const std::string source = osgDB::getNameLessExtension(filename);
            osg::ref_ptr<osgDB::Options> sourceOptions =
                osgEarth::Registry::cloneOrCreateOptions(options);

            const unsigned cacheHints =
                static_cast<unsigned>(sourceOptions->getObjectCacheHint()) &
                ~static_cast<unsigned>(osgDB::Options::CACHE_IMAGES);
            sourceOptions->setObjectCacheHint(
                static_cast<osgDB::Options::CacheHintOptions>(cacheHints));
            sourceOptions->removePluginData("osgEarth::URIResultCache");

            osgEarth::ReadResult result =
                osgEarth::URI(source).readImage(sourceOptions.get());

            if (result.succeeded() && result.getImage())
            {
                osg::ref_ptr<osg::Image> image = result.releaseImage();
                image->flipVertical();
                return image.release();
            }

            return result.code() == osgEarth::ReadResult::RESULT_NOT_FOUND ?
                ReadResult::FILE_NOT_FOUND :
                ReadResult::ERROR_IN_READING_FILE;
        }
    };

    REGISTER_OSGPLUGIN(flipvertical, FlipVerticalPseudoLoader)
}
