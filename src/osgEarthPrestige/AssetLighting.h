/* osgEarth
 * Copyright 2026 Pelican Mapping
 * MIT License
 */
#pragma once
#include <osgEarthPrestige/Export>
#include <osg/StateSet>

namespace osgEarthPrestige
{
    //! Installs two-sided leaf lighting and view-corrected baked crown normals before scene publication.
    //! Volume normals retain their orientation with reduced Prestige specular reflection; baked cards fade edge-on.
    OSGEARTHPRESTIGE_EXPORT void installAssetLighting(osg::StateSet* state);
}
