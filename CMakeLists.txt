CMAKE_MINIMUM_REQUIRED(VERSION 3.0.2)
PROJECT(pubgeo)

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(PUBGEO_INSTALL_ALIGN3D "Install align3d executable." OFF)
option(PUBGEO_INSTALL_SHR3D "Install shr3d executable." OFF)

IF(WIN32)
    # Default install directory
    MESSAGE(WARNING "Advance INSTALL OSGeo4W64. Include: (Command Line Tools) pdal and its dependencies. Potentially requires change to PDALConfig.cmake.")
    SET(CMAKE_PREFIX_PATH "C:\\OSGeo4W64\\")

    SET(BOOST_INCLUDEDIR "C:\\OSGeo4W64\\include\\boost-1_55")
    IF(NOT EXISTS ${BOOST_INCLUDEDIR})
        message(SEND_ERROR "Boost directory ${BOOST_INCLUDEDIR} does not exist!")
    ENDIF()

    SET(CMAKE_BUILD_TYPE RelWithDebInfo)

    # This part is for auto-defines in windows libraries that cause macro errors in our code
    add_definitions(-DWIN32_LEAN_AND_MEAN -DNOMINMAX)
ELSE()
    SET(CMAKE_INSTALL_PREFIX "/usr")
ENDIF()

FIND_PACKAGE(GDAL 1.9.0 REQUIRED)# Using Open Source Geo for Windows installer
if (GDAL_FOUND)
    include_directories(SYSTEM "${GDAL_INCLUDE_DIR}")
    mark_as_advanced(CLEAR GDAL_INCLUDE_DIR)
    mark_as_advanced(CLEAR GDAL_LIBRARY)
else()
    message(FATAL_ERROR "GDAL support is required")
endif()
FIND_PACKAGE(PDAL 1.4.0 REQUIRED)
LINK_DIRECTORIES(${PDAL_LIBRARY_DIRS})
IF(WIN32)
    # PDAL also has a requirement for winsock dlls
    set(PDAL_LIBRARIES ${PDAL_LIBRARIES} ws2_32)
ENDIF()

if (CMAKE_COMPILER_IS_GNUCXX)
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -O3" )
ENDIF()

# This makes linking/including a little cleaner looking later
SET(DAL_LIBS ${GDAL_LIBRARY} ${PDAL_LIBRARIES})
SET(DAL_INCLUDE_DIRS ${GDAL_INCLUDE_DIR} ${PDAL_INCLUDE_DIRS})
MESSAGE(STATUS "${DAL_INCLUDE_DIRS}")

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

ADD_SUBDIRECTORY(src/common)
ADD_SUBDIRECTORY(src/align3d)
ADD_SUBDIRECTORY(src/shr3d)
