// Copyright 2017 The Johns Hopkins University Applied Physics Laboratory.
// Licensed under the MIT License. See LICENSE.txt in the project root for full license information.

#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <random>
#include "align3d.h"

namespace align3d {
    bool computeRMS(float dx, float dy, long numSamples, long maxSamples, std::vector<double> &xlist,
                    std::vector<double> &ylist, OrthoImage<unsigned short> &referenceDSM,
                    OrthoImage<unsigned short> &targetDSM, float &medianDZ, float &rms, long &ndx,
                    float &completeness) {
        // Loop on number of samples and compute RMS.
        // Just in case there aren't many valid points, don't let this loop forever.
        long count = 0;
        ndx = 0;
        std::vector<float> differences;
        while ((count < numSamples) && (ndx < maxSamples)) {
            // Get the next point for matching.
            double x = xlist[ndx];
            double y = ylist[ndx];
            ndx++;

            // Map the point into the target.
            // Skip if this isn't a valid point.
            long col = long((x - targetDSM.easting + 0.5*dx) / targetDSM.gsd);
            long row = targetDSM.height - 1 - long((y - targetDSM.northing + 0.5*dy) / targetDSM.gsd);
            if (col < 0) continue;
            if (row < 0) continue;
            if (col >= targetDSM.width - 1) continue;
            if (row >= targetDSM.height - 1) continue;
            if (targetDSM.data[row][col] == 0) continue;
            float targetZ = float(targetDSM.data[row][col]) * targetDSM.scale + targetDSM.offset;

            // Map the point into the reference.
            // Skip if this isn't a valid point.
            col = long((x - referenceDSM.easting - 0.5*dx) / referenceDSM.gsd);
            row = referenceDSM.height - 1 - long((y - referenceDSM.northing - 0.5*dy) / referenceDSM.gsd);
            if (col < 0) continue;
            if (row < 0) continue;
            if (col >= referenceDSM.width - 1) continue;
            if (row >= referenceDSM.height - 1) continue;
            if (referenceDSM.data[row][col] == 0) continue;
            float referenceZ = float(referenceDSM.data[row][col]) * referenceDSM.scale + referenceDSM.offset;

            // Keep going until we have enough points.
            float difference = referenceZ - targetZ;
            differences.push_back(difference);
            count++;
        }

        // Skip if not enough sampled points.
        if (count < numSamples) return false;

        // Compute median Z offset and a robust estimate of the RMS difference.
        rms = 0.0;
        sort(differences.begin(), differences.end());
        medianDZ = differences[count / 2];
        for (long k = 0; k < count; k++) {
            differences[k] = fabs(differences[k] - medianDZ);
        }
        sort(differences.begin(), differences.end());
        rms = differences[(long) (count * 0.67)];

        // Compute the completeness.
        long good = 0;
        for (long k = 0; k < count; ++k) {
            if (differences[k] < 1.0) ++good;
        }
        completeness = good / (float) numSamples;

        return true;
    }

// Estimate 3D rigid body transform parameters to align target points with reference.
    void EstimateRigidBody(OrthoImage<unsigned short> &referenceDSM, OrthoImage<unsigned short> &targetDSM, float maxt,
                           AlignBounds &bounds, AlignResult &result) {
        float step = MIN(referenceDSM.gsd, targetDSM.gsd);
        long numSamples = 10000;
        long maxSamples = numSamples * 10;

        // Allocate RMS array.
        maxt = step * ceil(maxt / step);
        long bins = long(maxt / step * 2) + 1;
        float **rmsArray = new float *[bins];
        for (long i = 0; i < bins; i++) rmsArray[i] = new float[bins];

        // Get random samples.
        std::mt19937_64 rgen(0); // Create random number generator with seed of 0
        std::uniform_real_distribution<double> randx(bounds.xmin, bounds.xmax);
        std::uniform_real_distribution<double> randy(bounds.ymin, bounds.ymax);

        std::vector<double> xlist;
        std::vector<double> ylist;
        for (long i = 0; i < maxSamples; i++) {
            // Get a random point in the overlap area.
            xlist.push_back(randx(rgen));
            ylist.push_back(randy(rgen));
        }

        // Start with brute force, but sample points to reduce timeline.
        // Loop on X and Y translation within search distance.
        float threshold = MAX_FLOAT;
        float bestDX = 0.0;
        float bestDY = 0.0;
        float bestDZ = 0.0;
        long besti = 0;
        long bestj = 0;
        float bestRMS = MAX_FLOAT;
        float medianDZ = 0.0;
        bestRMS = threshold;
        float bestCompleteness = 0.0;
        long numSampled = 0;
        for (long i = 0; i < bins; i++) {
            float dx = -maxt + i * step;
            for (long j = 0; j < bins; j++) {
                float dy = -maxt + j * step;
                rmsArray[i][j] = 0.0;
                float rms = 0.0;
                float completeness = 0.0;
                bool ok = computeRMS(dx, dy, numSamples, maxSamples, xlist, ylist, referenceDSM, targetDSM, medianDZ,
                                     rms, numSampled, completeness);
                if (!ok) continue;
                rmsArray[i][j] = rms;
                if (rms < bestRMS || (rms == bestRMS && (dx*dx+dy*dy < bestDX*bestDX+bestDY*bestDY))) {
                    bestCompleteness = completeness;

                    bestRMS = rms;
                    bestDX = dx;
                    bestDY = dy;
                    bestDZ = medianDZ;
                    besti = i;
                    bestj = j;
                }
            }
        }

        // Apply quadratic interpolation to localize the peak.
        if ((besti > 0) && (besti < bins - 1) && (bestj > 0) && (bestj < bins - 1)) {
            float dx = (rmsArray[besti + 1][bestj] - rmsArray[besti - 1][bestj]) / 2.f;
            float dy = (rmsArray[besti][bestj + 1] - rmsArray[besti][bestj - 1]) / 2.f;
            float dxx = (rmsArray[besti + 1][bestj] + rmsArray[besti - 1][bestj] - 2 * rmsArray[besti][bestj]);
            float dyy = (rmsArray[besti][bestj + 1] + rmsArray[besti][bestj - 1] - 2 * rmsArray[besti][bestj]);
            float dxy =
                    (rmsArray[besti + 1][bestj + 1] - rmsArray[besti + 1][bestj - 1] - rmsArray[besti - 1][bestj + 1] +
                     rmsArray[besti - 1][bestj - 1]) / 4.f;
            float det = dxx * dyy - dxy * dxy;
            if (det != 0.0) {
                float ix = besti - (dyy * dx - dxy * dy) / det;
                float iy = bestj - (dxx * dy - dxy * dx) / det;
                bestDX = -maxt + ix * step;
                bestDY = -maxt + iy * step;
            }
        }

        // Deallocate RMS array.
        for (long i = 0; i < bins; i++) delete[]rmsArray[i];
        delete[]rmsArray;

        // Update the result and return.
        result.rms = bestRMS;
        result.tx = -bestDX;
        result.ty = -bestDY;
        result.tz = bestDZ;

        printf("Percent less than 1m Z difference = %6.2f%%\n", bestCompleteness * 100.0);
        printf("X offset = %f m\n", result.tx);
        printf("Y offset = %f m\n", result.ty);
        printf("Z offset = %f m\n", result.tz);
        printf("Z RMS    = %f m\n", result.rms);
    }

    // Loads file into raster, and applies common filterings to it
    bool load_file(OrthoImage<unsigned short> &dsm, const char *inputFileName, const AlignParameters &params) {
        int len = (int) strlen(inputFileName);
        const char *ext = &inputFileName[len - 3];
        if (strcmp(ext, "tif") == 0) {
            printf("File Type = .%s; loading as raster\n", ext);
            if (!dsm.read(inputFileName)) {
                printf("Failed to read %s\n", inputFileName);
                return false;
            }
        } else {
            // Convert the point cloud file to a DSM.
            printf("File Type = .%s; loading as point cloud\n", ext);
            if (!dsm.readFromPointCloud(inputFileName, params.gsd, MAX_VALUE)) {
                printf("Failed to read %s\n", inputFileName);
                return false;
            }
        }
        // Fill small voids.
        dsm.fillVoidsPyramid(true, 2);
        printf("Filtering data\n");
        // Remove points along edges which are difficult to match.
        dsm.edgeFilter((long) (params.maxdz / dsm.scale));
        return true;
    }
    
// Align target file to match reference file.
    bool AlignTarget2Reference(char *referenceFileName, char *targetFileName, AlignParameters params) {
        printf("Reading reference file: %s\n", referenceFileName);
        OrthoImage<unsigned short> referenceDSM;
        if (!load_file(referenceDSM, referenceFileName, params))
            return false;
        if (params.gsd != referenceDSM.gsd) {
            printf("Changing gsd to %f to match reference DSM\n",referenceDSM.gsd);
            params.gsd = referenceDSM.gsd;
        }
        
        printf("Reading target file: %s\n", targetFileName);
        OrthoImage<unsigned short> targetDSM;
        if (!load_file(targetDSM, targetFileName, params))
            return false;

        // Get overlapping bounds.
        AlignBounds bounds;
        bounds.xmin = MAX(referenceDSM.easting, targetDSM.easting);
        bounds.ymin = MAX(referenceDSM.northing, targetDSM.northing);
        bounds.xmax = MIN(referenceDSM.easting + (referenceDSM.width * referenceDSM.gsd),
                          targetDSM.easting + (targetDSM.width * targetDSM.gsd));
        bounds.ymax = MIN(referenceDSM.northing + (referenceDSM.height * referenceDSM.gsd),
                          targetDSM.northing + (targetDSM.height * targetDSM.gsd));
        bounds.width = bounds.xmax - bounds.xmin;
        bounds.height = bounds.ymax - bounds.ymin;
        double overlap_km = bounds.width / 1000.0 * bounds.height / 1000.0;
        printf("Overlap = %ld m x %ld m = %f km\n", (long) bounds.width, (long) bounds.height, overlap_km);
        if (overlap_km == 0.0) return false;

        // Estimate rigid body transform to align target points to reference.
        AlignResult result;
        printf("Estimating rigid body transformation.\n");
        EstimateRigidBody(referenceDSM, targetDSM, params.maxt, bounds, result);

        // Write offsets text file.
        printf("Writing offsets text file.\n");
        int len = strlen(targetFileName);
        char outFileName[1024];
        sprintf(outFileName, "%s", targetFileName);
        sprintf(&outFileName[len - 4], "_offsets.txt");
        FILE *fptr = fopen(outFileName, "w");
        if (fptr) {
            fprintf(fptr, "X Offset  Y Offset  Z Offset  Z RMS\n");
            fprintf(fptr, "%08.3f  %08.3f  %08.3f  %08.3f\n", result.tx, result.ty, result.tz, result.rms);
            fclose(fptr);
        } else {
            printf("Failed to write %s\n", outFileName);
            return false;
        }

        // Write aligned TIF file.
        printf("Writing aligned TIF file.\n");
        sprintf(&outFileName[len - 4], "_aligned.tif");
        targetDSM.offset += result.tz;
        targetDSM.easting += result.tx;
        targetDSM.northing += result.ty;
        if (!targetDSM.write(outFileName, true)) {
            printf("Failed to write %s\n", outFileName);
            return false;
        }

        // Write aligned BPF file.
        // For now, this requires an extra read of the point cloud file.
        printf("Writing aligned LAS file.\n");
        sprintf(&outFileName[len - 4], "_aligned.las");
        if (!PointCloud::TransformPointCloud(targetFileName, outFileName, result.tx, result.ty, result.tz)) {
            printf("Failed to write %s\n", outFileName);
            return false;
        }

        return true;
    }
}
