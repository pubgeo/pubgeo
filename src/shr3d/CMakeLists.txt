# Copyright (c) 2017-2019, The Johns Hopkins University /
# Applied Physics Laboratory (JHU/APL)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Please reference the following when reporting any results using this software:
#
# M. Bosch, A. Leichtman, D. Chilcott, H. Goldberg, M. Brown, “Metric
# Evaluation Pipeline for 3D Modeling of Urban Scenes,” ISPRS Archives, 2017.
#
# S. Almes, S. Hagstrom, D. Chilcott, H. Goldberg, M. Brown, “Open Source
# Geospatial Tools to Enable Large Scale 3D Scene Modeling,” FOSS4G, 2017.
#
# For more information, please see: http://www.jhuapl.edu/pubgeo.html

SET(SHR3D_HEADER_FILES
        shr3d.h)

SET(SHR3D_SOURCE_FILES
        shr3d.cpp)

IF(WIN32)
    include_directories(${BOOST_INCLUDEDIR})
ENDIF()

ADD_LIBRARY(SHR3D_LIB STATIC ${SHR3D_HEADER_FILES} ${SHR3D_SOURCE_FILES})
TARGET_LINK_LIBRARIES(SHR3D_LIB
        PUBLIC
        PUBGEO_LIB)

### Add some executables to play with.
### Main exe will showcase all features
ADD_EXECUTABLE(shr3d main.cpp)
TARGET_LINK_LIBRARIES(shr3d
        PUBLIC
        SHR3D_LIB)

if(PUBGEO_INSTALL_SHR3D)
    INSTALL(TARGETS shr3d
        RUNTIME DESTINATION bin)
endif()

IF(TESTING)
    ## This executable will be used to test IO
    ADD_EXECUTABLE(shr3dTest test-main.cpp)
    TARGET_LINK_LIBRARIES(shr3dTest
            PUBLIC
            SHR3D_LIB)
ENDIF()

IF(WIN32)
    SET(SHR3D_WRITER_NAME libpdal_plugin_writer_shr3d)
ELSE(WIN32)
    SET(SHR3D_WRITER_NAME pdal_plugin_writer_shr3d)
ENDIF(WIN32)

ADD_LIBRARY(${SHR3D_WRITER_NAME} SHARED plugin.cpp)
TARGET_LINK_LIBRARIES(${SHR3D_WRITER_NAME} ${PDAL_LIBRARIES} SHR3D_LIB)
SET_TARGET_PROPERTIES(${SHR3D_WRITER_NAME} PROPERTIES SOVERSION "0.1.0")
INSTALL(TARGETS ${SHR3D_WRITER_NAME}
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
