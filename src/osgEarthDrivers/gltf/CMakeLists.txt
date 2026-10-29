SET(TARGET_H
    GLTFReader.h
    GLTFWriter.h
    B3DMReader.h
    B3DMWriter.h
)

SET(TARGET_SRC
    ReaderWriterGLTF.cpp
)

add_osgearth_plugin(
    TARGET osgdb_gltf
    HEADERS ${TARGET_H}
    SOURCES ${TARGET_SRC} )
    
include_directories(
    ${OSGEARTH_EMBEDDED_THIRD_PARTY_DIR}/tinygltf 
    ${OSGEARTH_EMBEDDED_THIRD_PARTY_DIR}/rapidjson/include/rapidjson )

find_package(draco QUIET)
if(TARGET draco::draco)
    target_link_libraries(osgdb_gltf PRIVATE draco::draco)
    target_compile_definitions(osgdb_gltf PRIVATE OSGEARTH_HAVE_DRACO)
endif()
