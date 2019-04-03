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

#include <cstdio>
#include <chrono>
#include <regex>
#include "orthoimage.h"
#include "shr3d.h"

// Print command line arguments.
void printArguments() {
    printf("Command line arguments: <Input File (LAS|TIF)> <Options>\n");
    printf("Required Options:\n");
    printf("  DH=    horizontal uncertainty (meters)\n");
    printf("  DZ=    vertical uncertainty (meters)\n");
    printf("  AGL=   minimum building height above ground level (meters)\n");
    printf("Options:\n");
    printf("  AREA=  minimum building area (meters)\n");
    printf("  EGM96  set this flag to write vertical datum = EGM96\n");
    printf("  BOUNDS=MINX,MAXX,MINY,MAXY set to define image bounds\n");
    printf("  DTM=   Path to optional input DTM geotiff file\n");
    printf("  GND_LABEL= Set to the value of the point cloud ground classification label (typically 2)\n");
    printf("             if the point cloud has already been partially classified.\n");
    printf("Examples:\n");
    printf("  For EO DSM:    shr3d dsm.tif DH=5.0 DZ=1.0 AGL=2 AREA=50.0 EGM96\n");
    printf("  For lidar DSM: shr3d dsm.tif DH=1.0 DZ=1.0 AGL=2 AREA=50.0\n");
    printf("  For lidar LAS: shr3d pts.las DH=1.0 DZ=1.0 AGL=2 AREA=50.0\n");
}


// Main program for bare earth classification.
int main(int argc, char **argv) {
    // If no parameters, then print command line arguments.
    if (argc < 4) {
        printf("Number of arguments = %d\n", argc);
        for (int i = 0; i < argc; i++) {
            printf("ARG[%d] = %s\n", i, argv[i]);
        }
        printArguments();
        return -1;
    }

    // Get command line arguments and confirm they are valid.
    shr3d::Shr3dder shr3dder;
    bool convert = false;
    char inputFileName[1024] = {0};
    strcpy(inputFileName, argv[1]);
    std::string dtmFileName;
    for (int i = 2; i < argc; i++) {
        if (strstr(argv[i], "DH=")) { shr3dder.dh_meters = atof(&(argv[i][3])); }
        if (strstr(argv[i], "DZ=")) { shr3dder.dz_meters = atof(&(argv[i][3])); }
        if (strstr(argv[i], "AGL=")) { shr3dder.agl_meters = atof(&(argv[i][4])); }
        if (strstr(argv[i], "AREA=")) { shr3dder.min_area_meters = atof(&(argv[i][5])); }
        if (strstr(argv[i], "EGM96")) { shr3dder.egm96 = true; }
        if (strstr(argv[i], "BOUNDS=")) {
            std::string values(&(argv[i][7]));
            std::regex fmt("^(.+),(.+),(.+),(.+)$");
            if (!std::regex_match(values,fmt)) {
                printf("Error: Bounds are in an incorrect format.\n");
                printArguments();
                return -1;
            }
            std::istringstream str_bnds(std::regex_replace(values,fmt,"([$1,$2],[$3,$4])"));
            str_bnds >> shr3dder.bounds;
        }
        if (strstr(argv[i], "CONVERT")) { convert = true; }
        if (strstr(argv[i], "GND_LABEL=")) { shr3dder.gnd_label = atoi(&(argv[i][10])); }
        if (strstr(argv[i], "DTM=")) { dtmFileName = std::string(&(argv[i][4])); }
    }
    if ((shr3dder.dh_meters == 0.0) || (shr3dder.dz_meters == 0.0) || (shr3dder.agl_meters == 0.0)) {
        printf("DH_METERS = %f\n", shr3dder.dh_meters);
        printf("DZ_METERS = %f\n", shr3dder.dz_meters);
        printf("AGL_METERS = %f\n", shr3dder.agl_meters);
        printf("Error: All three values must be nonzero.\n");
        printArguments();
        return -1;
    }

    // Initialize the timer.
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();

    // If specified, then convert to GDAL TIFF.
    char readFileName[1024];
    strcpy(readFileName, inputFileName);
    if (convert) {
        char cmd[4096];
        sprintf(readFileName, "temp.tif");
        sprintf(cmd, ".\\gdal\\gdal_translate %s temp.tif\n", inputFileName);
        int rval = system(cmd);
        if (rval != 0) return rval;
    }

    std::map<shr3d::ImageType,std::string> outputFilenames;
    std::string basename(inputFileName);

    // Read DSM as SHORT.
    // Input can be GeoTIFF or LAS/BPF.
    char *ext = strrchr(inputFileName,'.');
    printf("File Type = %s\n", ext);
    if (strcmp(ext, ".tif") == 0) {
        if (!shr3dder.setDSM(readFileName))
            return -1;
    } else {
        if (!shr3dder.setPSET(inputFileName))
            return -1;

        outputFilenames[shr3d::DSM] = basename + "_DSM.tif";
        outputFilenames[shr3d::INTENSITY] = basename + "_INT.tif";
    }
    
    if (dtmFileName.empty()) {
        outputFilenames[shr3d::DTM] = basename + "_DTM.tif";
    } else {
        if (!shr3dder.setDTM0(dtmFileName))
            return -1;
    }

    // Convert horizontal and vertical uncertainty values to bin units.
    int dh_bins = MAX(1, (int) floor(shr3dder.dh_meters / shr3dder.getDSM().gsd));
    printf("DZ_METERS = %f\n", shr3dder.dz_meters);
    printf("DH_METERS = %f\n", shr3dder.dh_meters);
    printf("DH_BINS = %d\n", dh_bins);
    unsigned int dz_short = (unsigned int) (shr3dder.dz_meters / shr3dder.getDSM().scale);
    printf("DZ_SHORT = %d\n", dz_short);
    printf("AGL_METERS = %f\n", shr3dder.agl_meters);
    unsigned int agl_short = (unsigned int) (shr3dder.agl_meters / shr3dder.getDSM().scale);
    printf("AGL_SHORT = %d\n", agl_short);
    printf("AREA_METERS = %f\n", shr3dder.min_area_meters);

    // Set outputs
#ifdef DEBUG
    outputFilenames[shr3d::DSM] = basename + "_DSM.tif";
    outputFilenames[shr3d::MIN] = basename + "_MIN.tif";
    outputFilenames[shr3d::DSM2] = basename + "_DSM2.tif";
    outputFilenames[shr3d::MINAGL] = basename + "_MINAGL.tif";
    outputFilenames[shr3d::DTM0] = basename + "_DTM0.tif";
    outputFilenames[shr3d::DTM] = basename + "_DTM.tif";
    outputFilenames[shr3d::LABEL0] = basename + "_label0.tif";
    outputFilenames[shr3d::LABEL] = basename + "_label.tif";
    outputFilenames[shr3d::LABELED_BUILDINGS] = basename + "_building_labels.tif";
    outputFilenames[shr3d::LABELED_BUILDINGS_3] = basename + "_building_labels_3.tif";
#endif
    outputFilenames[shr3d::CLASS] = basename + "_class.tif";
    outputFilenames[shr3d::BUILDING] = basename + "_buildings.tif";
    outputFilenames[shr3d::BUILDING_OUTLINES] = basename + "_buildings.shp";

    // Process
    try {
        shr3dder.process(outputFilenames);
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    // Report total elapsed time.
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    printf("Total time elapsed = %f seconds\n", std::chrono::duration<double>(t1-t0).count());
}
