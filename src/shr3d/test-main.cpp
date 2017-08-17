// Copyright 2017 The Johns Hopkins University Applied Physics Laboratory.
// Licensed under the MIT License. See LICENSE.txt in the project root for full license information.

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