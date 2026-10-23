/* osgEarth
 * Copyright 2026 Pelican Mapping
 * MIT License
 */
#include "AssetLighting.h"
#include <osgEarth/VirtualProgram>

void osgEarthPrestige::installAssetLighting(osg::StateSet* state)
{
    state->setDefine("OE_CHONK_BAKED_CROWN");
    state->setDefine("OE_PRESTIGE_SKY_SPECULAR_FACTOR");
    osgEarth::VirtualProgram::getOrCreate(state)->setFunction("oe_p2_leaf_normal", R"glsl(
        in vec3 vp_Normal;
        flat in uint oe_normal_technique;
        float oe_prestige_sky_specularFactor;
        void oe_p2_leaf_normal(inout vec4 color)
        {
            // Crown-volume normals stand in for many differently oriented leaves. A full surface
            // reflection makes their shared lobe wash out whole branches under a grazing sun.
            // Match authored crowns (3) and their baked proxies (5); ordinary surfaces and VRV
            // source-normal impostors retain their material response.
            oe_prestige_sky_specularFactor = (oe_normal_technique == 3u || oe_normal_technique == 5u) ? 0.1 : 1.0;

            // Only real leaf surfaces turn their normal on the back face. An impostor's
            // volume normals describe the crown, so turning them would point the crown downward.
            if (!gl_FrontFacing && oe_normal_technique == 0u) vp_Normal = -vp_Normal;
        }
    )glsl", osgEarth::VirtualProgram::LOCATION_FRAGMENT_LIGHTING, -1.0f);
}
