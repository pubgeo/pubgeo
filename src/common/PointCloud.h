// Copyright (c) 2017-2019, The Johns Hopkins University /
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

//
// Created by almess1 on 4/20/17.
// This file was created in order to adapt the SHR3D and ALIGN3D to open source IO.
//

#ifndef PUBGEO_NOT_POINT_SETS_H
#define PUBGEO_NOT_POINT_SETS_H

#include <pdal/PointView.hpp>
#include <pdal/PipelineExecutor.hpp>

namespace pubgeo {
    struct MinMaxXYZ {
        double xMin;
        double xMax;
        double yMin;
        double yMax;
        double zMin;
        double zMax;
    };


    class PointCloud {
    public:
        PointCloud();

        ~PointCloud();

        static bool TransformPointCloud(std::string inputFileName, std::string outputFileName,
                                        float translateX, float translateY, float translateZ);

        PointCloud CropToClass(int keep_class);

        bool Read(const char *fileName);
        bool Read(pdal::PointViewPtr view);

        MinMaxXYZ bounds;
        int zone;
        unsigned long numPoints;
        int xOff;
        int yOff;
        int zOff;

        inline float x(int i) const {
            if (pv != nullptr)
                return (float) (pv->getFieldAs<double>(pdal::Dimension::Id::X, i) - xOff);
            throw "Point set has not been initialized.";
        }

        inline float y(int i) const {
            if (pv != nullptr)
                return (float) (pv->getFieldAs<double>(pdal::Dimension::Id::Y, i) - yOff);
            throw "Point set has not been initialized.";
        }

        inline float z(int i) const {
            if (pv != nullptr)
                return (float) (pv->getFieldAs<double>(pdal::Dimension::Id::Z, i) - zOff);
            throw "Point set has not been initialized.";
        }

        inline float i(int i) const {
            if (pv != nullptr)
                return pv->getFieldAs<float>(pdal::Dimension::Id::Intensity, i);
            throw "Point set has not been initialized.";
        }

        inline int c(int i) const {
            if (pv != nullptr)
                return pv->getFieldAs<int>(pdal::Dimension::Id::Classification, i);
            throw "Point set has not been initialized.";
        }

    private:
        pdal::PipelineExecutor *executor;
        pdal::PointViewPtr pv;

        void CleanupPdalPointers();
    };
}
#endif //PUBGEO_NOT_POINT_SETS_H
