// Copyright 2017 The Johns Hopkins University Applied Physics Laboratory.
// Licensed under the MIT License. See LICENSE.txt in the project root for full license information.

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
