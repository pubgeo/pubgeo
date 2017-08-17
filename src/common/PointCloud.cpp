// Copyright 2017 The Johns Hopkins University Applied Physics Laboratory.
// Licensed under the MIT License. See LICENSE.txt in the project root for full license information.

#include "PointCloud.h"

// Pipeline needs to read in point cloud file of any type, and read it in. ideally in meters
static std::string PDAL_PIPELINE_OPEN_ENGINE = R"({ "pipeline": [ ")";
static std::string PDAL_PIPELINE_OPEN_CABOOSE = R"("] } )";

std::string buildPipelineStr(const char *fileName) {
    return PDAL_PIPELINE_OPEN_ENGINE + fileName + PDAL_PIPELINE_OPEN_CABOOSE;
}

namespace pubgeo {
    PointCloud::PointCloud() : executor(nullptr), pv(nullptr), zone(0), xOff(0), yOff(0), zOff(0), numPoints(0) {
        bounds = MinMaxXYZ{0, 0, 0, 0, 0, 0};
    }

    PointCloud::~PointCloud() {
        CleanupPdalPointers();
    }

    bool PointCloud::Read(const char *fileName) {
        try {
            CleanupPdalPointers();
            executor = new pdal::PipelineExecutor(buildPipelineStr(fileName));
            executor->execute();
            const pdal::PointViewSet &pvs = executor->getManagerConst().views();
            pv = pvs.begin()->get();
            if (pvs.size() > 1) {
                std::cerr << "[PUBGEO::PointCloud::READ] File contains additional unread sets." << std::endl;
            }

            numPoints = pv->size();
            if (numPoints < 1) {
                std::cerr << "[PUBGEO::PointCloud::READ] No points found in file." << std::endl;
                return false;
            }

            pdal::BOX3D box;
            pv->calculateBounds(box);
            pdal::SpatialReference sr = pv->spatialReference();
            zone = sr.computeUTMZone(box);

            // used later to return points
            xOff = (int) floor(box.minx);
            yOff = (int) floor(box.miny);
            zOff = (int) floor(box.minz);

            bounds = {box.minx, box.maxx, box.miny, box.maxy, box.minz, box.maxz};
            return true;
        }
        catch (pdal::pdal_error &pe) {
            std::cerr << pe.what() << std::endl;
            return false;
        }
    }

    bool PointCloud::TransformPointCloud(const char *inputFileName, const char *outputFileName,
                                         float translateX = 0, float translateY = 0, float translateZ = 0) {
        std::ostringstream pipeline;
        pipeline << "{\n\t\"pipeline\":[\n\t\t\"" << inputFileName
                 << "\",\n\t\t{\n\t\t\t\"type\":\"filters.transformation\",\n"
                 << "\t\t\t\"matrix\":\""
                 << "1 0 0 " << translateX << " "
                 << "0 1 0 " << translateY << " "
                 << "0 0 1 " << translateZ << " "
                 << "0 0 0 1\"\n\t\t},\n\t\t{"
                 << "\n\t\t\t\"filename\":\"" << outputFileName << "\"\n\t\t}\n\t]\n}";
        try {
            const std::string pipe = pipeline.str();
            pdal::PipelineExecutor executor(pipe);
            executor.execute();
        } catch (pdal::pdal_error &pe) {
            std::cerr << pe.what() << std::endl;
            return false;
        }
        return true;
    }

    void PointCloud::CleanupPdalPointers() {
        if (pv != nullptr) {
            pv = nullptr;
        }
        if (executor != nullptr) {
            delete executor;
            executor = nullptr;
        }
    }

};
