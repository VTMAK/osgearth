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
#include "BiomeManager"

#include <osgEarth/TextureArena>
#include <osgEarth/Elevation>
#include <osgEarth/GLUtils>
#include <osgEarth/Metrics>
#include <osgEarth/MaterialLoader>
#include <osgEarth/Utils>
#include <osg/ComputeBoundsVisitor>
#include <osg/MatrixTransform>

using namespace osgEarth;
using namespace osgEarth::Procedural;

#undef LC
#define LC "[BiomeManager] "

#define NORMAL_MAP_TEX_UNIT 1

//...................................................................

ResidentModelAsset::Ptr
ResidentModelAsset::create()
{
    return Ptr(new ResidentModelAsset());
}

ResidentModelAsset::ResidentModelAsset() :
    _assetDef(nullptr)
{
    //nop
}

//...................................................................

namespace
{
    osg::Image* convertNormalMapFromRGBToRG(const osg::Image* in)
    {
        osg::Image* out = new osg::Image();
        out->allocateImage(in->s(), in->t(), 1, GL_RG, GL_UNSIGNED_BYTE);
        out->setInternalTextureFormat(GL_RG8);
        osg::Vec4 v;
        osg::Vec4 packed;
        ImageUtils::PixelReader read(in);
        ImageUtils::PixelWriter write(out);

        ImageUtils::ImageIterator i(read);
        i.forEachPixel([&]()
            {
                read(v, i.s(), i.t());
                osg::Vec3 normal(v.r()*2.0f - 1.0f, v.g()*2.0f - 1.0f, v.b()*2.0f - 1.0f);
                NormalMapGenerator::pack(normal, packed);
                write(packed, i.s(), i.t());
            }
        );

        return out;
    }
}

//...................................................................

BiomeManager::BiomeManager() :
    _revision(0),
    _refsAndRevision_mutex("BiomeManager.refsAndRevision(OE)"),
    _residentData_mutex("BiomeManager.residentData(OE)"),
    _lodTransitionPixelScale(16.0f),
    _locked(false)
{
    // this arena will hold all the textures for loaded assets.
    _textures = new TextureArena();
    _textures->setName("Biomes");
    _textures->setBindingPoint(1);
}

void
BiomeManager::ref(const Biome* biome)
{
    ScopedMutexLock lock(_refsAndRevision_mutex);

    auto item = _refs.emplace(biome, 0);
    ++item.first->second;
    if (item.first->second == 1) // ref count of 1 means it's new
    {
        ++_revision;
        OE_INFO << LC << "Hello, " << biome->name().get() << " (" << biome->index() << ")" << std::endl;
    }
}

void
BiomeManager::unref(const Biome* biome)
{
    ScopedMutexLock lock(_refsAndRevision_mutex);

    auto iter = _refs.find(biome);
    
    // silent assertion
    if (iter == _refs.end() || iter->second == 0)
    {
        OE_SOFT_ASSERT(iter != _refs.end() && iter->second != 0);
        return;
    }

    if (!_locked)
    {
        --iter->second;
        if (iter->second == 0)
        {
            // Technically yes, this creates a new revision, but it
            // would also result in re-loading assets that are already
            // resident... think on this -gw
            //++_revision;
            OE_INFO << LC << "Goodbye, " << biome->name().get() << "(" << biome->index() << ")" << std::endl;
        }
    }

    //OE_INFO << LC << "(" << iter->second << ")" << biome->name().get() << std::endl;
}

int
BiomeManager::getRevision() const
{
    return _revision;
}

void
BiomeManager::setLODTransitionPixelScale(float value)
{
    _lodTransitionPixelScale = value;

    // We need to rebuild the assets, so clear everything out,
    // recalculate the biome set, and bump the revision.
    _residentData_mutex.lock();
    _residentModelAssets.clear();
    _residentData_mutex.unlock();

    _refsAndRevision_mutex.lock();
    ++_revision;
    _refsAndRevision_mutex.unlock();
}

float
BiomeManager::getLODTransitionPixelScale() const
{
    return _lodTransitionPixelScale;
}

void
BiomeManager::reset()
{
    // reset the reference counts, and bump the revision so the
    // next call to update will remove any resident data
    {
        ScopedMutexLock lock(_refsAndRevision_mutex);

        for (auto& iter : _refs)
        {
            const Biome* biome = iter.first;
            OE_INFO << LC << "Goodbye, " << biome->name().get() << std::endl;
            iter.second = 0;
        }

        ++_revision;
    }

    // Resolve the references and unload any resident assets from memory.
    recalculateResidentBiomes();
}

void
BiomeManager::recalculateResidentBiomes()
{
    std::vector<const Biome*> biomes_to_add;
    std::vector<const Biome*> biomes_to_remove;

    // Figure out which biomes we need to load and which we can discard.
    {
        ScopedMutexLock lock(_refsAndRevision_mutex);

        for (auto& ref : _refs)
        {
            const Biome* biome = ref.first;
            int refcount = ref.second;

            if (refcount > 0)
            {
                biomes_to_add.push_back(biome);
            }
            else if (!_locked)
            {
                biomes_to_remove.push_back(biome);
            }
        }
    }

    // Update the resident biome data structure:
    {
        ScopedMutexLock lock(_residentData_mutex);

        // add biomes that might need adding
        for (auto biome : biomes_to_add)
        {
            auto& dummy = _residentBiomes[biome];
        }

        // get rid of biomes we no longer need.
        for (auto biome : biomes_to_remove)
        {
            _residentBiomes.erase(biome);
        }


        // finally, update the collection of resident assets to
        // reflect the reference counts.
        std::list<const ModelAsset*> to_delete;

        for (auto& entry : _residentModelAssets)
        {
            if (entry.second.use_count() == 1)
            {
                to_delete.push_back(entry.first);
            }
        }

        for (auto& asset : to_delete)
        {
            _residentModelAssets.erase(asset);
            OE_DEBUG << LC << "Unloaded asset " << asset->name() << std::endl;
        }
    }
}

std::vector<const Biome*>
BiomeManager::getActiveBiomes() const
{
    ScopedMutexLock lock(_refsAndRevision_mutex);

    std::vector<const Biome*> result;

    for (auto& ref : _refs)
    {
        if (ref.second > 0)
            result.push_back(ref.first);
    }
    return std::move(result);
}

std::vector<const ModelAsset*>
BiomeManager::getResidentAssets() const
{
    std::vector<const ModelAsset*> result;

    ScopedMutexLock lock(_residentData_mutex);

    result.reserve(_residentModelAssets.size());
    for (auto& entry : _residentModelAssets)
        result.push_back(entry.first);

    return std::move(result);
}

std::vector<const ModelAsset*>
BiomeManager::getResidentAssetsIfNotLocked() const
{
    std::vector<const ModelAsset*> result;
    if (_residentData_mutex.try_lock())
    {
        result.reserve(_residentModelAssets.size());
        for (auto& entry : _residentModelAssets)
            result.push_back(entry.first);

        _residentData_mutex.unlock();
    }
    return std::move(result);
}

namespace
{
    struct ModelCacheEntry
    {
        osg::ref_ptr<osg::Node> _node;
        osg::BoundingBox _modelAABB;
    };
}

void
BiomeManager::materializeNewAssets(
    const osgDB::Options* readOptions)
{
    OE_PROFILING_ZONE;

    // VRV_PATCH
    // The VRV DtOsgFileCache will fail to cache the textures from our vegetation
    // models correctly, and they will show up without textures. It prints out
    // this message:
    //    "Warning: makVrv::DtOsgFileCache::ts_readImageFile: Attempting to save
    //    accelerated file for [imagename] which also has myInstanceNodesFlag set
    //    to false. These two options are incompatible and the accelerated file
    //    will not be saved."
    // We work around this by disabling the use of the acceleration cache.
    osg::ref_ptr<osgDB::Options> local =
        readOptions ? osg::clone(readOptions, osg::CopyOp::DEEP_COPY_ALL) :
        new osgDB::Options();
    std::string os = local->getOptionString();
    Strings::ciReplaceIn(os, "vrvUseAcceleratedCache", "");
    local->setOptionString(os);
    readOptions = local.get();
    // END VRV_PATCH

    // exclusive access to the resident dataset
    ScopedMutexLock lock(_residentData_mutex);

    std::set<std::string> asset_groups;

    // Any billboard that doesn't have its own normal map will use this one.
    osg::ref_ptr<osg::Texture> defaultNormalMap = osgEarth::createEmptyNormalMapTexture();

    // Some caches to avoid duplicating data
    using TextureShareCache = std::map<URI, ResidentModelAsset::Ptr>;
    TextureShareCache texcache;
    using ModelCache = std::map<URI, ModelCacheEntry>;
    ModelCache modelcache;

    // Factory for loading chonk data. It will use our texture arena.
    ChonkFactory factory(_textures.get());

    // Clear out each biome's instances so we can start fresh.
    // This is a low-cost operation since anything we can re-use
    // will already by in the _residentModelAssetData collection.
    OE_DEBUG << LC << "Found " << _residentBiomes.size() << " resident biomes..." << std::endl;
    for (auto& iter : _residentBiomes)
    {
        iter.second.clear();
    }

    // This loader will find material textures and install them on
    // secondary texture image units.. in this case, normal maps.
    // We can expand this later to include other types of material maps.
    Util::MaterialLoader materialLoader;

    materialLoader.setMangler(NORMAL_MAP_TEX_UNIT,
        [](const std::string& filename)
        {
            return osgDB::getNameLessExtension(filename) + "_NML.png";
        }
    );

    materialLoader.setTextureFactory(NORMAL_MAP_TEX_UNIT,
        [](osg::Image* image)
        {
            osg::Texture2D* tex = nullptr;

            // compresses the incoming texture if necessary
            if (image->getPixelFormat() != GL_RG)
                tex = new osg::Texture2D(convertNormalMapFromRGBToRG(image));
            else
                tex = new osg::Texture2D(image);

            return tex;
        }
    );

    // Go through the residency list and materialize any model assets
    // that are not already loaded (and present in _residentModelAssets);
    // Along the way, build the instances for each biome.
    for (auto& e : _residentBiomes)
    {
        const Biome* biome = e.first;

        auto& assetsToUse = biome->getModelAssets();

        // this is the collection of instances we're going to populate
        ResidentModelAssetInstances& assetInstances = e.second;

        // The group points to multiple assets, which we will analyze and load.
        for (auto& asset_ptr : assetsToUse)
        {
            const ModelAsset* assetDef = asset_ptr->asset();

            OE_SOFT_ASSERT(assetDef != nullptr);
            if (assetDef == nullptr)
                continue;

            // Look up this model asset. If it's already in the resident set,
            // do nothing; otherwise make it resident by loading all the data.
            ResidentModelAsset::Ptr& residentAsset = _residentModelAssets[assetDef];

            // First reference to this instance? Populate it:
            if (residentAsset == nullptr)
            {
                OE_INFO << LC << "  Loading asset " << assetDef->name() << std::endl;

                residentAsset = ResidentModelAsset::create();

                residentAsset->assetDef() = assetDef;

                osg::BoundingBox bbox;

                if (assetDef->modelURI().isSet())
                {
                    const URI& uri = assetDef->modelURI().get();

                    ModelCache::iterator ic = modelcache.find(uri);
                    if (ic != modelcache.end())
                    {
                        residentAsset->model() = ic->second._node.get();
                        residentAsset->boundingBox() = ic->second._modelAABB;
                    }
                    else
                    {
                        residentAsset->model() = uri.getNode(readOptions);
                        if (residentAsset->model().valid())
                        {
                            // apply a static scale:
                            if (assetDef->scale().isSet())
                            {
                                osg::MatrixTransform* mt = new osg::MatrixTransform();
                                mt->setMatrix(osg::Matrix::scale(assetDef->scale().get(), assetDef->scale().get(), assetDef->scale().get()));
                                mt->addChild(residentAsset->model().get());
                                residentAsset->model() = mt;
                            }

                            // find materials:
                            materialLoader.setReferrer(uri.full());
                            residentAsset->model()->accept(materialLoader);

                            // add flexors:
                            bool isUndergrowth = assetDef->group() == "undergrowth";
                            addFlexors(residentAsset->model(), assetDef->stiffness().get(), isUndergrowth);

                            OE_DEBUG << LC << "Loaded model: " << uri.base() << std::endl;
                            modelcache[uri]._node = residentAsset->model().get();

                            osg::ComputeBoundsVisitor cbv;
                            residentAsset->model()->accept(cbv);
                            residentAsset->boundingBox() = cbv.getBoundingBox();
                            modelcache[uri]._modelAABB = residentAsset->boundingBox();
                        }
                        else
                        {
                            OE_WARN << LC << "Failed to load model " << uri.full() << std::endl;
                        }
                    }

                    bbox = residentAsset->boundingBox();
                }

                URI sideBB;
                if (assetDef->sideBillboardURI().isSet())
                    sideBB = assetDef->sideBillboardURI().get();
                else if (assetDef->modelURI().isSet())
                    sideBB = URI(assetDef->modelURI()->full() + ".side.png", assetDef->modelURI()->context());
                
                float impostorFarPixelScale = 1.0f;

                if (!sideBB.empty())
                {
                    const URI& uri = sideBB;

                    auto ic = texcache.find(uri);
                    if (ic != texcache.end())
                    {
                        residentAsset->sideBillboardTex() = ic->second->sideBillboardTex();
                        residentAsset->sideBillboardNormalMap() = ic->second->sideBillboardNormalMap();
                    }
                    else
                    {
                        osg::ref_ptr<osg::Image> image = uri.getImage(readOptions);
                        if (image.valid())
                        {
                            residentAsset->sideBillboardTex() = new osg::Texture2D(image.get());

                            OE_DEBUG << LC << "Loaded BB: " << uri.base() << std::endl;
                            texcache[uri] = residentAsset;

                            osg::ref_ptr<osg::Image> normalMap;
                            for (int i = 0; i < 2 && !normalMap.valid(); ++i)
                            {
                                // normal map is the same file name but with _NML inserted before the extension
                                URI normalMapURI(
                                    osgDB::getNameLessExtension(uri.full()) +
                                    (i == 0 ? "_NML." : ".normal.") +
                                    osgDB::getFileExtension(uri.full()));

                                // silenty fail if no normal map found.
                                normalMap = normalMapURI.getImage(readOptions);
                                if (normalMap.valid())
                                {
                                    residentAsset->sideBillboardNormalMap() = new osg::Texture2D(normalMap.get());
                                }
                            }
                        }
                        else
                        {
                            OE_WARN << LC << "Failed to load side billboard " << uri.full() << std::endl;
                        }
                    }

                    URI topBB;
                    if (assetDef->topBillboardURI().isSet())
                        topBB = assetDef->topBillboardURI().get();
                    else if (assetDef->modelURI().isSet())
                        topBB = URI(assetDef->modelURI()->full() + ".top.png", assetDef->modelURI()->context());

                    if (!topBB.empty())
                    {
                        const URI& uri = topBB;

                        auto ic = texcache.find(uri);
                        if (ic != texcache.end())
                        {
                            residentAsset->topBillboardTex() = ic->second->topBillboardTex();
                            residentAsset->topBillboardNormalMap() = ic->second->topBillboardNormalMap();
                        }
                        else
                        {
                            osg::ref_ptr<osg::Image> image = uri.getImage(readOptions);
                            if (image.valid())
                            {
                                residentAsset->topBillboardTex() = new osg::Texture2D(image.get());

                                OE_DEBUG << LC << "Loaded BB: " << uri.base() << std::endl;
                                texcache[uri] = residentAsset;

                                osg::ref_ptr<osg::Image> normalMap;
                                for (int i = 0; i < 2 && !normalMap.valid(); ++i)
                                {
                                    // normal map is the same file name but with _NML inserted before the extension
                                    URI normalMapURI(
                                        osgDB::getNameLessExtension(uri.full()) +
                                        (i == 0 ? "_NML." : ".normal.") +
                                        osgDB::getFileExtension(uri.full()));

                                    // silenty fail if no normal map found.
                                    normalMap = normalMapURI.getImage(readOptions);
                                    if (normalMap.valid())
                                    {
                                        residentAsset->topBillboardNormalMap() = new osg::Texture2D(normalMap.get());
                                    }
                                }
                            }
                            else
                            {
                                OE_WARN << LC << "Failed to load top billboard " << uri.full() << std::endl;
                            }
                        }
                    }

                    std::vector<osg::Texture*> textures(4);
                    textures[0] = residentAsset->sideBillboardTex().get();
                    textures[1] = residentAsset->sideBillboardNormalMap().get();
                    textures[2] = residentAsset->topBillboardTex().get();
                    textures[3] = residentAsset->topBillboardNormalMap().get();

                    // if this group has an impostor creation function, call it
                    auto iter = _createImpostorFunctions.find(assetDef->group());
                    if (iter != _createImpostorFunctions.end())
                    {
                        auto& createImpostor = iter->second;

                        if (createImpostor != nullptr)
                        {
                            if (!bbox.valid())
                            {
                                bbox.set(
                                    -residentAsset->assetDef()->width().get(),
                                    -residentAsset->assetDef()->width().get(),
                                    0,
                                    residentAsset->assetDef()->width().get(),
                                    residentAsset->assetDef()->width().get(),
                                    residentAsset->assetDef()->height().get());
                            }

                            Impostor imp = createImpostor(
                                bbox,
                                textures);

                            residentAsset->impostor() = imp._node;
                            impostorFarPixelScale = imp._farPixelScale;

                            // if no MAIN model exists, just re-use the impostor as the 
                            // main model!
                            if (!residentAsset->model().valid())
                            {
                                residentAsset->model() = residentAsset->impostor().get();
                            }
                        }
                    }
                }

                // Finally, chonkify.
                if (residentAsset->model().valid())
                {
                    // Replace impostor with 3D model "N" times closer than the SSE:
                    float far_pixel_scale = getLODTransitionPixelScale();

                    // Real model can get as close as it wants:
                    float near_pixel_scale = FLT_MAX;

                    if (residentAsset->chonk() == nullptr)
                        residentAsset->chonk() = Chonk::create();

                    residentAsset->chonk()->add(
                        residentAsset->model().get(),
                        far_pixel_scale,
                        near_pixel_scale,
                        factory);
                }

                if (residentAsset->impostor().valid())
                {
                    // Fade out the impostor at this pixel scale (i.e., multiple of SSE)
                    float far_pixel_scale = impostorFarPixelScale;

                    // Fade the impostor to the 3D model at this pixel scale:
                    float near_pixel_scale = residentAsset->model().valid() ?
                        getLODTransitionPixelScale() : 
                        FLT_MAX;

                    if (residentAsset->chonk() == nullptr)
                        residentAsset->chonk() = Chonk::create();

                    residentAsset->chonk()->add(
                        residentAsset->impostor().get(),
                        far_pixel_scale,
                        near_pixel_scale,
                        factory);
                }
            }

            // If this data successfully materialized, add it to the
            // biome's instance collection.
            if (residentAsset->sideBillboardTex().valid() || residentAsset->model().valid())
            {
                ResidentModelAssetInstance instance;
                instance.residentAsset() = residentAsset;
                instance.weight() = asset_ptr->weight();
                instance.coverage() = asset_ptr->coverage();
                assetInstances.push_back(std::move(instance));
            }
        }   
    }
}

void
BiomeManager::setCreateImpostorFunction(
    const std::string& group,
    BiomeManager::CreateImpostorFunction func)
{
    ScopedMutexLock lock(_residentData_mutex);
    _createImpostorFunctions[group] = func;
}

BiomeManager::ResidentBiomes
BiomeManager::getResidentBiomes(
    const osgDB::Options* readOptions)
{
    OE_PROFILING_ZONE;

    // First refresh the resident biome collection based on current refcounts
    recalculateResidentBiomes();

    // Next go through and load any assets that are not yet loaded
    materializeNewAssets(readOptions);

    // Make a copy:
    ResidentBiomes result = _residentBiomes;

    return std::move(result);
}

namespace
{
}

void
BiomeManager::addFlexors(
    osg::ref_ptr<osg::Node>& node,
    float stiffness,
    bool isUndergrowth)
{
    // note: stiffness is [0..1].

    osg::ComputeBoundsVisitor cb;
    node->accept(cb);
    auto& bbox = cb.getBoundingBox();

    if (isUndergrowth)
    {
        // for undergrowth, modify the normals so they all point in the skyward
        // direction of the blade card
        auto fixNormals = [&bbox](
            osg::Geometry& geom,
            unsigned i0, unsigned i1, unsigned i2,
            const osg::Matrix& l2w)
        {
            auto verts = dynamic_cast<osg::Vec3Array*>(geom.getVertexArray());
            if (!verts) return;

            // make the proper normal array on demand:
            osg::Vec3Array* normals = dynamic_cast<osg::Vec3Array*>(geom.getNormalArray());
            if (!normals || normals->getBinding() != osg::Array::BIND_PER_VERTEX || normals->size() != verts->size())
            {
                normals = new osg::Vec3Array(osg::Array::BIND_PER_VERTEX, verts->size());
                geom.setNormalArray(normals);
            }

            // make a flexor array on demand
            osg::Vec3Array* flexors = dynamic_cast<osg::Vec3Array*>(geom.getTexCoordArray(3));
            if (!flexors)
            {
                flexors = new osg::Vec3Array(osg::Array::BIND_PER_VERTEX, verts->size());
                geom.setTexCoordArray(3, flexors);
            }

            osg::Vec3 v0 = (*verts)[i0] * l2w;
            osg::Vec3 v1 = (*verts)[i1] * l2w;
            osg::Vec3 v2 = (*verts)[i2] * l2w;

            osg::Vec3 tri_normal = (v1 - v0) ^ (v2 - v0);
            tri_normal.normalize();

            const osg::Vec3 zaxis(0, 0, 1);
            const osg::Vec3 xaxis(1, 0, 0);
            const osg::Vec3 yaxis(0, 1, 0);

            osg::Vec3 tangent = tri_normal ^ zaxis;
            osg::Vec3 normal = tri_normal ^ tangent;

            if (normal.z() < 0.0f)
                normal = -normal;

            normal.normalize();

            // back to local space
            if (!l2w.isIdentity())
            {
                osg::Matrix w2l;
                w2l.invert(l2w);
                normal = osg::Matrix::transform3x3(w2l, normal);
            }

            (*normals)[i0] = normal;
            (*normals)[i1] = normal;
            (*normals)[i2] = normal;

            (*flexors)[i0] = normal * accel(v0.z() / bbox.zMax());
            (*flexors)[i1] = normal * accel(v1.z() / bbox.zMax());
            (*flexors)[i2] = normal * accel(v2.z() / bbox.zMax());

        };

        TriangleVisitor visitor(fixNormals);
        node->accept(visitor);
    }

    else // not undergrorwth (normal tree, bush models)
    {
        auto addFlexors = [&bbox, stiffness](osg::Geometry& geom, const osg::Matrix& l2w)
        {
            osg::Vec3Array* verts = dynamic_cast<osg::Vec3Array*>(geom.getVertexArray());
            if (!verts) return;

            osg::Vec3Array* flexors = new osg::Vec3Array(osg::Array::BIND_PER_VERTEX);
            flexors->reserve(verts->size());
            geom.setTexCoordArray(3, flexors);

            osg::Vec3f base(bbox.center().x(), bbox.center().y(), bbox.zMin());
            float xy_radius = std::max(bbox.xMax() - bbox.xMin(), bbox.yMax() - bbox.yMin());

            for (auto& local_vert : *verts)
            {
                auto vert = local_vert * l2w;
                osg::Vec3f anchor(0.0f, 0.0f, vert.z());
                float height_ratio = harden((vert.z() - bbox.zMin()) / (bbox.zMax() - bbox.zMin()));
                auto lateral_vec = (vert - osg::Vec3f(0, 0, vert.z()));
                float lateral_ratio = harden(lateral_vec.length() / xy_radius);
                auto flex =
                    normalize(lateral_vec) *
                    (lateral_ratio * 3.0f) *
                    (height_ratio * 1.5f) *
                    (1.0 - stiffness);

                flexors->push_back(flex);
            }
        };

        TypedNodeVisitor<osg::Geometry> visitor(addFlexors);
        node->accept(visitor);
    }
}
