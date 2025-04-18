/* osgEarth
 * Copyright 2025 Pelican Mapping
 * MIT License
 */
#include <SilverLining.h>
#include "SilverLiningSkyDrawable"
#include "SilverLiningContext"
#include "SilverLiningContextNode"
#include <osgEarth/SpatialReference>
#include <osgEarth/GLUtils>

#define LC "[SilverLining:SkyDrawable] "

using namespace osgEarth::SilverLining;


SkyDrawable::SkyDrawable(SilverLiningContextNode* contexNode) :
_SL(contexNode->getSLContext()),
_contextNode(contexNode)

{
    // call this to ensure draw() gets called every frame.
    setSupportsDisplayList( false );

    // not MT-safe (camera updates, etc)
    this->setDataVariance( osg::Object::DYNAMIC );

    setName("SilverLining::SkyDrawable");;
}

void
SkyDrawable::drawImplementation(osg::RenderInfo& renderInfo) const
{
    OE_GL_ZONE;

    osg::Camera* camera = renderInfo.getCurrentCamera();
#ifndef SL_USE_CULL_MASK
	//Check if this is the target camera
	if (_contextNode->getTargetCamera() == camera) 
#endif 
	{
	if ( camera)
    {
        renderInfo.getState()->disableAllVertexArrays();

        _SL->initialize( renderInfo );

        // convey the sky box size (far plane) to SL:
        double fovy, ar, znear, zfar;
        _SL->setCamera(camera);
        camera->getProjectionMatrixAsPerspective(fovy, ar, znear, zfar);
        _SL->setSkyBoxSize( zfar < 100000.0 ? zfar : 100000.0 );

        // invoke the user callback if it exists
        if (_SL->getCallback())
            _SL->getCallback()->onDrawSky(_SL->getAtmosphereWrapper());

        osg::Matrix projMat = renderInfo.getState()->getProjectionMatrix();
        _SL->getAtmosphere()->SetProjectionMatrix(projMat.ptr());
        osg::Matrix viewMat = renderInfo.getCurrentCamera()->getViewMatrix();
        _SL->getAtmosphere()->SetCameraMatrix(viewMat.ptr());

        // draw the sky.
        _SL->getAtmosphere()->DrawSky(
            true,
            _SL->getSRS()->isGeographic(),
            _SL->getSkyBoxSize(),
            true,
            false );

        // Dirty the state and the program tracking to prevent GL state conflicts.
        renderInfo.getState()->dirtyAllVertexArrays();
        renderInfo.getState()->dirtyAllAttributes();

        // Reset the saved program.  SilverLining exits its functionality with a glUseProgram(0). Without this line,
        // GL Core 3.3 rendering will attempt to load uniforms without an active program, which is an error.  This
        // tells the state that there is currently no installed program, so if it needs one, to load one.
        renderInfo.getState()->setLastAppliedProgramObject(0L);
        renderInfo.getState()->apply();
    }
	}
}

osg::BoundingBox
SkyDrawable::computeBoundingBox() const
{
    osg::BoundingBox skyBoundBox;
    if ( !_SL->ready() )
        return skyBoundBox;

    ::SilverLining::Atmosphere* atmosphere = _SL->getAtmosphere();
    double skyboxSize = _SL->getSkyBoxSize();
    if ( skyboxSize == 0.0 )
        skyboxSize = 1000.0;

    osg::Vec3d radiusVec = osg::Vec3d(skyboxSize, skyboxSize, skyboxSize) * 0.5;
    osg::Vec3d camPos = _SL->getCameraPosition();
    if (_SL->getCamera())
    {
        osg::Vec3f eye, center, up;
        _SL->getCamera()->getViewMatrixAsLookAt(eye, center, up);
        camPos = osg::Vec3d(eye.x(), eye.y(), eye.z());
    }

    skyBoundBox.set( camPos-radiusVec, camPos+radiusVec );

    // this enables the "blue ring" around the earth when viewing from space.
    bool hasLimb = atmosphere->GetConfigOptionBoolean("enable-atmosphere-from-space");
    if ( hasLimb )
    {
        // Compute bounds of atmospheric limb centered at (0,0,0)
        double earthRadius = atmosphere->GetConfigOptionDouble("earth-radius-meters");
        double atmosphereHeight = earthRadius + atmosphere->GetConfigOptionDouble("atmosphere-height");
        double atmosphereThickness = atmosphere->GetConfigOptionDouble("atmosphere-scale-height-meters") + earthRadius;

        osg::BoundingBox atmosphereBox;
        osg::Vec3d atmMin(-atmosphereThickness, -atmosphereThickness, -atmosphereThickness);
        osg::Vec3d atmMax(atmosphereThickness, atmosphereThickness, atmosphereThickness);
        atmosphereBox.set( atmMin, atmMax );
        skyBoundBox.expandBy( atmosphereBox );
    }
    return skyBoundBox;
}
