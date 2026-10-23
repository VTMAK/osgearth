#include <osgEarth/MapNode>
#include <osgEarth/TMS>
#include <osgEarth/EarthManipulator>
#include <osgEarth/GLUtils>
#include <osgViewer/Viewer>

int main(int argc, char** argv)
{
    osg::ArgumentParser args(&argc, argv);

    osgEarth::initialize();
    osgEarth::setNotifyLevel(osg::INFO);

    osgViewer::Viewer viewer(args);
    viewer.setRealizeOperation(new osgEarth::GL3RealizeOperation());
    viewer.setCameraManipulator(new osgEarth::Util::EarthManipulator(args));    
    
    auto layer = new osgEarth::TMSImageLayer();
    layer->setURL("https://readymap.org/readymap/tiles/1.0.0/7/");
    
    auto mapNode = new osgEarth::MapNode();
    mapNode->getMap()->addLayer(layer);    
    
    viewer.setSceneData(mapNode);
    return viewer.run();
}