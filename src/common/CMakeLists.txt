SET(PUBGEO_HEADER_FILES
        util.h
        Image.h
        orthoimage.h
        PointCloud.h
        disjoint_set.h
        geo_polygon.h)

SET(PUBGEO_SOURCE_FILES
        PointCloud.cpp
        )

FIND_PACKAGE(Threads REQUIRED)

ADD_LIBRARY(PUBGEO_LIB STATIC ${PUBGEO_HEADER_FILES} ${PUBGEO_SOURCE_FILES})
TARGET_INCLUDE_DIRECTORIES(PUBGEO_LIB
        PUBLIC
        ${DAL_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR})
TARGET_LINK_LIBRARIES(PUBGEO_LIB
        PUBLIC
        ${DAL_LIBS}
        ${CMAKE_THREAD_LIBS_INIT})