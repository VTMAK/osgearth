
IF(OSGEARTH_BUILD_SHARED_LIBS)
    ADD_DEFINITIONS(-DOSGEARTHBUILDINGS_LIBRARY)
ELSE ()
    ADD_DEFINITIONS(-DOSGEARTHBUILDINGS_LIBRARY_STATIC)
ENDIF()

set(HEADERS
    Analyzer
    BuildContext
    Building
    BuildingCatalog
    BuildingCompiler
    BuildingFactory
    BuildingLayer
    BuildingLayerOptions
    BuildingPager
    BuildingSymbol
    BuildingVisitor
    Common
    Compiler
    CompilerOutput
    CompilerSettings
    InstancedBuildingCompiler
    InstancedRoofCompiler
    Elevation
    ElevationCompiler
    Export
    FlatRoofCompiler
    GableRoofCompiler
    Parapet
    Roof
    Zoning
	ExternalModelNode
	InstancedModelNode
	ElevationsLodNode
	FilterUsage
)

set(SOURCES
    Analyzer.cpp
    Building.cpp
    BuildingCatalog.cpp
    BuildingCompiler.cpp
    BuildingFactory.cpp
    BuildingLayer.cpp
    BuildingPager.cpp
    BuildingSymbol.cpp
    BuildingVisitor.cpp
    Compiler.cpp
    CompilerOutput.cpp
    CompilerSettings.cpp
    InstancedBuildingCompiler.cpp
    InstancedRoofCompiler.cpp
    Elevation.cpp
    ElevationCompiler.cpp
    FeaturePlugin.cpp
    FlatRoofCompiler.cpp
    GableRoofCompiler.cpp
    Parapet.cpp
    Roof.cpp
    ExternalModelNode.cpp
    InstancedModelNode.cpp
    ElevationsLodNode.cpp
)

add_osgearth_library(
    TARGET osgEarthBuildings
    PUBLIC_HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    FOLDER "NodeKits" )
