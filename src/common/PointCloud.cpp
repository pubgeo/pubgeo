// Copyright 2017 The Johns Hopkins University Applied Physics Laboratory.
// Licensed under the MIT License. See LICENSE.txt in the project root for full license information.

#include "PointCloud.h"
#ifdef WIN32
#include <regex>
#endif

// Pipeline needs to read in point cloud file of any type, and read it in. ideally in meters
static std::string PDAL_PIPELINE_OPEN_ENGINE = R"({ "pipeline": [ ")";
static std::string PDAL_PIPELINE_OPEN_CABOOSE = R"("] } )";

std::string buildPipelineStr(std::string fileName) {
#ifdef WIN32
	fileName = std::regex_replace(fileName, std::regex("\\\\"), "/"); // Replace single backslash with double
#endif
    return PDAL_PIPELINE_OPEN_ENGINE + fileName + PDAL_PIPELINE_OPEN_CABOOSE;
}

namespace pubgeo {
    PointCloud::PointCloud() : zone(0), numPoints(0), xOff(0), yOff(0), zOff(0), executor(nullptr), pv(nullptr) {
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
            if (pvs.size() > 1) {
                std::cerr << "[PUBGEO::PointCloud::READ] File contains additional unread sets." << std::endl;
            }

            return Read(*pvs.begin());
        }
        catch (pdal::pdal_error &pe) {
            std::cerr << pe.what() << std::endl;
            return false;
        }
    }

    bool PointCloud::Read(pdal::PointViewPtr view) {
        pv = view;
        numPoints = pv->size();
        if (numPoints < 1) {
            std::cerr << "[PUBGEO::PointCloud::READ] No points found in file." << std::endl;
            return false;
        }

        pdal::BOX3D box;
        pv->calculateBounds(box);
        zone = pv->spatialReference().getUTMZone();

        // used later to return points
        xOff = (int) floor(box.minx);
        yOff = (int) floor(box.miny);
        zOff = (int) floor(box.minz);

        bounds = {box.minx, box.maxx, box.miny, box.maxy, box.minz, box.maxz};
        return true;
    }


    bool PointCloud::TransformPointCloud(std::string inputFileName, std::string outputFileName,
                                         float translateX = 0, float translateY = 0, float translateZ = 0) {
#ifdef WIN32
		inputFileName = std::regex_replace(inputFileName, std::regex("\\\\"), "/"); // Replace single backslash with double
		outputFileName = std::regex_replace(outputFileName, std::regex("\\\\"), "/"); // Replace single backslash with double
#endif
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

    PointCloud PointCloud::CropToClass(int keep_class) {
        pdal::PointViewPtr outView = pv->makeNew();
        for (pdal::PointId idx = 0; idx < pv->size(); ++idx)
            if (c(idx)==keep_class)
                outView->appendPoint(*pv,idx);
        PointCloud out;
        out.Read(outView);
        return out;
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
