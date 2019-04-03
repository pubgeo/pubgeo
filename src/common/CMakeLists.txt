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