CMAKE_MINIMUM_REQUIRED(VERSION 3.0.2)
PROJECT(pubgeo)

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_POSITION_INDEPENDENT_CODE ON)

IF(WIN32)
       MESSAGE(ERROR "Windows builds are currently not working due to erroneous linking errors.")
       # Default install directory
       MESSAGE(WARNING "Advance INSTALL OSGeo4W64. Include: (Command Line Tools) pdal and its dependencies. Potentially requires change to PDALConfig.cmake.")
       SET(CMAKE_PREFIX_PATH "C:\\OSGeo4W64\\")

       # This part is for auto-defines in windows libraries that cause macro errors in our code
       add_definitions(-DWIN32_LEAN_AND_MEAN -DNOMINMAX)
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
