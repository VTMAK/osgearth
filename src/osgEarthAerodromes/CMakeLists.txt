
IF(OSGEARTH_BUILD_SHARED_LIBS)
    ADD_DEFINITIONS(-DOSGEARTHAERODROME_LIBRARY)
ELSE (OSGEARTH_BUILD_SHARED_LIBS)
    ADD_DEFINITIONS(-DOSGEARTHAERODROME_LIBRARY_STATIC)
ENDIF(OSGEARTH_BUILD_SHARED_LIBS)

set(HEADERS
    AerodromeCatalog
    AerodromeFactory
    AerodromeFeatureNode
    AerodromeFeatureOptions
    AerodromeLayer
    AerodromeNode
    AerodromeRenderer
    BoundaryNode
    Common
    Export
    LightBeaconNode
    LightIndicatorNode
    LinearFeatureNode
    PavementNode
    RunwayNode
    RunwayThresholdNode
    StartupLocationNode
    StopwayNode
    TaxiwayNode
    TaxiwaySignNode
    TerminalNode
    WindsockNode )

set(SOURCES
    AerodromeCatalog.cpp
    AerodromeFactory.cpp
    AerodromeLayer.cpp
    AerodromeNode.cpp
    AerodromeRenderer.cpp )

add_osgearth_library(
    TARGET osgEarthAerodromes
    PUBLIC_HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    FOLDER "NodeKits" )
