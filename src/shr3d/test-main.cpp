// Copyright (c) 2017, The Johns Hopkins University /
// Applied Physics Laboratory (JHU/APL)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// Please reference the following when reporting any results using this software:
//
// M. Bosch, A. Leichtman, D. Chilcott, H. Goldberg, M. Brown, “Metric
// Evaluation Pipeline for 3D Modeling of Urban Scenes,” ISPRS Archives, 2017.
//
// S. Almes, S. Hagstrom, D. Chilcott, H. Goldberg, M. Brown, “Open Source
// Geospatial Tools to Enable Large Scale 3D Scene Modeling,” FOSS4G, 2017.
//
// For more information, please see: http://www.jhuapl.edu/pubgeo.html

#include <iostream>
#include <random>
#include "PointCloud.h"

int main(int argc, char *argv[]) {
    const char *fileName = "test.las";
    switch (argc) {
        case 2:
            fileName = argv[1];
        case 1: {
            pubgeo::PointCloud pointCloud;
            std::cerr << "Processing file: " << fileName << std::endl;
            if (pointCloud.Read(fileName)) {
                // Grab a truly random c++ 11 number
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_int_distribution<> dis(0, pointCloud.numPoints);
                long int index = dis(gen);
                std::cerr << pointCloud.numPoints << " points in file." << std::endl
                          << "Here is a sample point from index " << index << ":"
                          << pointCloud.x(index) << ", " << pointCloud.y(index) << ", "
                          << pointCloud.z(index) << std::endl;
            }
            break;
        }
        default:
            std::cerr << "Invalid number of arguments" << std::endl;
    }

#ifdef WIN32
    system("pause");
#endif
}