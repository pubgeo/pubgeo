// Copyright 2017 The Johns Hopkins University Applied Physics Laboratory.
// Licensed under the MIT License. See LICENSE.txt in the project root for full license information.

#include <cstdio>
#include <ctime>
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
    printf("  AREA=	 minimum building area (meters)\n");
    printf("  EGM96  set this flag to write vertical datum = EGM96\n");
    printf("Examples:\n");
    printf("  For EO DSM:    shr3d dsm.tif DH=5.0 DZ=1.0 AGL=2 AREA=50.0 EGM96\n");
    printf("  For lidar DSM: shr3d dsm.tif DH=1.0 DZ=1.0 AGL=2.0 AREA=50.0\n");
    printf("  For lidar LAS: shr3d pts.las DH=1.0 DZ=1.0 AGL=2.0 AREA=50.0\n");
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
    double dh_meters = 0.0;
    double dz_meters = 0.0;
    double agl_meters = 0.0;
    double min_area_meters = 50.0;
    bool egm96 = false;
    bool convert = false;
    char inputFileName[1024];
    strcpy(inputFileName, argv[1]);
    for (int i = 2; i < argc; i++) {
        if (strstr(argv[i], "DH=")) { dh_meters = atof(&(argv[i][3])); }
        if (strstr(argv[i], "DZ=")) { dz_meters = atof(&(argv[i][3])); }
        if (strstr(argv[i], "AGL=")) { agl_meters = atof(&(argv[i][4])); }
        if (strstr(argv[i], "AREA=")) { min_area_meters = atof(&(argv[i][5])); }
        if (strstr(argv[i], "EGM96")) { egm96 = true; }
        if (strstr(argv[i], "CONVERT")) { convert = true; }
    }
    if ((dh_meters == 0.0) || (dz_meters == 0.0) || (agl_meters == 0.0)) {
        printf("DH_METERS = %f\n", dh_meters);
        printf("DZ_METERS = %f\n", dz_meters);
        printf("AGL_METERS = %f\n", agl_meters);
        printf("Error: All three values must be nonzero.\n");
        printArguments();
        return -1;
    }

    // Initialize the timer.
    time_t t0;
    time(&t0);

    // If specified, then convert to GDAL TIFF.
    char readFileName[1024];
    strcpy(readFileName, inputFileName);
    if (convert) {
        char cmd[4096];
        sprintf(readFileName, "temp.tif");
        sprintf(cmd, ".\\gdal\\gdal_translate %s temp.tif\n", inputFileName);
        system(cmd);
    }

    // Read DSM as SHORT.
    // Input can be GeoTIFF or LAS.
    shr3d::OrthoImage<unsigned short> dsmImage;
    shr3d::OrthoImage<unsigned short> minImage;
    printf("Reading DSM as SHORT.\n");
    int len = (int) strlen(inputFileName);
    char *ext = &inputFileName[len - 3];
    printf("File Type = .%s\n", ext);
    if (strcmp(ext, "tif") == 0) {
        bool ok = dsmImage.read(readFileName);
        if (!ok) return -1;

        // Copy DSM
        minImage = dsmImage;

        // Min filter, replacing only points differing by more than the AGL threshold.
        minImage.minFilter(4, (unsigned int) (agl_meters / minImage.scale));

        // Fill small voids in the DSM.
        minImage.fillVoidsPyramid(true, 2);
#ifdef DEBUG
        // Write the MIN image as FLOAT.
        char minOutFileName[1024];
        sprintf(minOutFileName, "%s_MIN.tif", inputFileName);
        minImage.write(minOutFileName, true);
#endif
    } else if ((strcmp(ext, "las") == 0) || (strcmp(ext, "bpf") == 0)) {
        // First get the max Z values for the DSM.
        // Read a PSET file (e.g., BPF or LAS).
        shr3d::PointCloud pset;
        bool ok = pset.Read(inputFileName);
        if (!ok) return -1;

        ok = dsmImage.readFromPointCloud(pset, (float) dh_meters, shr3d::MAX_VALUE);
        if (!ok) return -1;

        // Median filter, replacing only points differing by more than the AGL threshold.
        dsmImage.medianFilter(1, (unsigned int) (agl_meters / dsmImage.scale));

        // Fill small voids in the DSM.
        dsmImage.fillVoidsPyramid(true, 2);

        // Write the DSM image as FLOAT.
        char dsmOutFileName[1024];
        sprintf(dsmOutFileName, "%s_DSM.tif", inputFileName);
        dsmImage.write(dsmOutFileName, true);

        // Now get the minimum Z values for the DTM.
        ok = minImage.readFromPointCloud(pset, (float) dh_meters, shr3d::MIN_VALUE);
        if (!ok) return -1;

        // Median filter, replacing only points differing by more than the AGL threshold.
        minImage.medianFilter(1, (unsigned int) (agl_meters / minImage.scale));

        // Fill small voids in the DSM.
        minImage.fillVoidsPyramid(true, 2);
#ifdef DEBUG
        // Write the MIN image as FLOAT.
        char minOutFileName[1024];
        sprintf(minOutFileName, "%s_MIN.tif", inputFileName);
        minImage.write(minOutFileName, true);
#endif
        // Find many of the trees by comparing MIN and MAX. Set their values to void.
        shr3d::OrthoImage<unsigned short> varImage = dsmImage - minImage;
        unsigned short maxTreeHeightScaled = (40.0 / minImage.scale);
        unsigned short threshold = (dz_meters / dsmImage.scale);

        // Apply tree filter
        shr3d::Image<unsigned short>::filter(&dsmImage, &varImage,
                [&](unsigned short* val, const unsigned short& ref, std::vector<unsigned short> &ngbrs) {
            // CAUTION: This is a hack to address an observed lidar sensor issue and may not generalize well.
            if (ref <= maxTreeHeightScaled) {
                // Set dsm to void if none of the neighbors are solid (var is < threshold)
                if (std::none_of(ngbrs.begin(), ngbrs.end(), [&](unsigned short v){ return v<=threshold; }))
                    *val = 0;
            }
        }, 1, 0, false);

        // Write the DSM2 image as FLOAT.
#ifdef DEBUG
        char dsm2OutFileName[1024];
        sprintf(dsm2OutFileName, "%s_DSM2.tif", inputFileName);
        dsmImage.write(dsm2OutFileName, true);
#endif
    } else {
        printf("Error: Unrecognized file type.");
        return -1;
    }

    // Convert horizontal and vertical uncertainty values to bin units.
    int dh_bins = MAX(1, (int) floor(dh_meters / dsmImage.gsd));
    printf("DZ_METERS = %f\n", dz_meters);
    printf("DH_METERS = %f\n", dh_meters);
    printf("DH_BINS = %d\n", dh_bins);
    unsigned int dz_short = (unsigned int) (dz_meters / dsmImage.scale);
    printf("DZ_SHORT = %d\n", dz_short);
    printf("AGL_METERS = %f\n", agl_meters);
    unsigned int agl_short = (unsigned int) (agl_meters / dsmImage.scale);
    printf("AGL_SHORT = %d\n", agl_short);
    printf("AREA_METERS = %f\n", min_area_meters);

    // Generate label image.
    shr3d::OrthoImage<unsigned long> labelImage;
    labelImage.Allocate(dsmImage.width, dsmImage.height);
    labelImage.easting = dsmImage.easting;
    labelImage.northing = dsmImage.northing;
    labelImage.zone = dsmImage.zone;
    labelImage.gsd = dsmImage.gsd;

    // Allocate a DTM image as SHORT and copy in the Min DSM values.
    shr3d::OrthoImage<unsigned short> dtmImage(minImage);

    // Classify ground points.
    shr3d::Shr3dder::classifyGround(labelImage, dsmImage, dtmImage, dh_bins, dz_short);

    // For DSM voids, also set DTM value to void.
    // Note: because we've changed the DSM by this point (setting voids where all the trees are),
    //  use the minImage which will have the same voids as the original DSM
    printf("Setting DTM values to VOID where DSM is VOID...\n");
    for (unsigned int j = 0; j < minImage.height; j++) {
        for (unsigned int i = 0; i < minImage.width; i++) {
            if (minImage.data[j][i] == 0) dtmImage.data[j][i] = 0;
        }
    }

    // Median filter, replacing only points differing by more than the DZ threshold.
    dtmImage.medianFilter(1, (unsigned int) (dz_meters / dsmImage.scale));

    // Refine the object label image and export building outlines.
    shr3d::Shr3dder::classifyNonGround(dsmImage, dtmImage, labelImage, dz_short, agl_short, (float) min_area_meters);

    // Fill small voids in the DTM after all processing is complete.
    dtmImage.fillVoidsPyramid(true, 2);

    // Write the DTM image as FLOAT.
    char dtmOutFileName[1024];
    sprintf(dtmOutFileName, "%s_DTM.tif", inputFileName);
    dtmImage.write(dtmOutFileName, true, egm96);

    // Produce a classification raster image with LAS standard point classes.
    shr3d::OrthoImage<unsigned char> classImage;
    classImage.Allocate(labelImage.width, labelImage.height);
    classImage.easting = dsmImage.easting;
    classImage.northing = dsmImage.northing;
    classImage.zone = dsmImage.zone;
    classImage.gsd = dsmImage.gsd;
    for (unsigned int j = 0; j < classImage.height; j++) {
        for (unsigned int i = 0; i < classImage.width; i++) {
            // Set default as unlabeled.
            classImage.data[j][i] = LAS_UNCLASSIFIED;

            // Label trees.
            if ((dsmImage.data[j][i] == 0) ||
                (fabs((float) dsmImage.data[j][i] - (float) dtmImage.data[j][i]) > agl_short))
                classImage.data[j][i] = LAS_TREE;

            // Label buildings.
            if (labelImage.data[j][i] == 1) classImage.data[j][i] = LAS_BUILDING;

            // Label ground.
            if (fabs((float) dsmImage.data[j][i] - (float) dtmImage.data[j][i]) < dz_short)
                classImage.data[j][i] = LAS_GROUND;
        }
    }

    // Fill missing labels inside building regions.
    shr3d::Shr3dder::fillInsideBuildings(classImage);

    // Write the classification image.
    char classOutFileName[1024];
    sprintf(classOutFileName, "%s_class.tif", inputFileName);
    classImage.write(classOutFileName, false, egm96);
    for (unsigned int j = 0; j < classImage.height; j++) {
        for (unsigned int i = 0; i < classImage.width; i++) {
            if (classImage.data[j][i] != LAS_BUILDING) classImage.data[j][i] = 0;
        }
    }
    sprintf(classOutFileName, "%s_buildings.tif", inputFileName);
    classImage.write(classOutFileName, false, egm96);

    // Report total elapsed time.
    time_t t1;
    time(&t1);
    printf("Total time elapsed = %f seconds\n", double(t1 - t0));
}
