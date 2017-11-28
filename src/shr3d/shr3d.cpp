// Copyright 2017 The Johns Hopkins University Applied Physics Laboratory.
// Licensed under the MIT License. See LICENSE.txt in the project root for full license information.

// shr3d.cpp
//

#include "shr3d.h"

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <climits>
#include <vector>

#include "geo_polygon.h"

using namespace shr3d;

void Shr3dder::process(const OrthoImage<unsigned short> &dsmImage, const OrthoImage<unsigned short> &minImage,
        std::map<ImageType,std::string> outputFilenames) {
    if (outputFilenames.empty())
        return;
    if (!outputFilenames[DSM].empty())
        dsmImage.write(outputFilenames[DSM].c_str(), true, egm96);
    if (!outputFilenames[MIN].empty())
        minImage.write(outputFilenames[MIN].c_str(), true, egm96);

    ImageType last_output = DSM;
    for (const std::pair<ImageType,std::string>& o : outputFilenames)
        if (o.first > last_output && o.first <= BUILDING_OUTLINES)
            last_output = o.first;

    if (last_output <= MIN)
        return;

    shr3d::OrthoImage<unsigned short> dsm2Image, dtmImage;
    shr3d::OrthoImage<unsigned long> labelImage;
    createDTM(dsmImage,minImage,dtmImage,dsm2Image,labelImage);

    if (!outputFilenames[DSM2].empty())
        dsm2Image.write(outputFilenames[DSM2].c_str(), true, egm96);
    if (!outputFilenames[DTM].empty())
        dtmImage.write(outputFilenames[DTM].c_str(), true, egm96);
    if (!outputFilenames[LABEL].empty())
        labelImage.write(outputFilenames[LABEL].c_str(), false);

    if (last_output <= DTM)
        return;

    shr3d::OrthoImage<unsigned char> classImage = labelClasses(dsmImage,dtmImage,dsm2Image,labelImage);

    if (!outputFilenames[CLASS].empty())
        classImage.write(outputFilenames[CLASS].c_str(), false);

    if (last_output <= CLASS)
        return;

    shr3d::OrthoImage<unsigned char> bldgImage = labelBuildings(classImage);

    if (!outputFilenames[BUILDING].empty())
        bldgImage.write(outputFilenames[BUILDING].c_str(), false);

    if (last_output <= BUILDING)
        return;

    shr3d::OrthoImage<int> bldgLabels(&bldgImage);
    bldgLabels.labelConnectedComponentsFrom(&bldgImage);

    if (!outputFilenames[LABELED_BUILDINGS].empty())
        bldgLabels.write(outputFilenames[LABELED_BUILDINGS].c_str(), false);

    shr3d::OrthoImage<int> bldgLabels3;
    bldgLabels3.nn_upsample(&bldgLabels,3);

    if (!outputFilenames[LABELED_BUILDINGS_3].empty())
        bldgLabels3.write(outputFilenames[LABELED_BUILDINGS_3].c_str(), false);

    std::map<int,GeoPolygon<double>> bounds = GeoPolygon<double>::traceBoundaries(bldgLabels3);
    printf("Traced %lu building outlines.\n",bounds.size());

    std::map<int,GeoPolygon<double>> new_bounds;
    for (std::pair<int,GeoPolygon<double>> pr : bounds) {
        GeoPolygon<double> poly = pr.second.buildingSimplify();
        if (poly.ring.empty())
            continue;
        new_bounds[pr.first] = poly;
    }
    printf("After simplifying, left with %lu building outlines.\n",new_bounds.size());

    if (!outputFilenames[BUILDING_OUTLINES].empty())
        GeoPolygon<double>::write(outputFilenames[BUILDING_OUTLINES],new_bounds);
}

bool Shr3dder::createDSM(const PointCloud& pset, OrthoImage<unsigned short> &dsmImage) {
    if (!dsmImage.readFromPointCloud(pset, (float) dh_meters, shr3d::MAX_VALUE))
        return false;

    // Median filter, replacing only points differing by more than the AGL threshold.
    dsmImage.quantileFilter(1, (unsigned int) (agl_meters / dsmImage.scale), 0.4);

    // Filter wells
    dsmImage.filter([&](unsigned short* val, const unsigned short& ref, std::vector<unsigned short> &ngbrs) {
        if (count_if(ngbrs.begin(),ngbrs.end(),[&](unsigned short ngbr) {
            return ngbr > ref + (sqrt(max_tree_height_meters*agl_meters) / dsmImage.scale); }) >= 4)
            *val = 0;
    });

    // Fill small voids in the DSM.
    dsmImage.fillVoidsPyramid(true, 2);

    return true;
}

bool Shr3dder::createMIN(const PointCloud& pset, OrthoImage<unsigned short> &minImage) {
    // Now get the minimum Z values for the DTM.
    if (!minImage.readFromPointCloud(pset, (float) dh_meters, shr3d::MIN_VALUE))
        return false;

    // Median filter, replacing only points differing by more than the AGL threshold.
    minImage.quantileFilter(2, (unsigned int) (agl_meters / minImage.scale), 0.33);

    // Fill small voids in the DSM.
    minImage.fillVoidsPyramid(true, 2);

    return true;
}

void Shr3dder::createDTM(const OrthoImage<unsigned short> &dsmImage, const OrthoImage<unsigned short> &minImage,
        OrthoImage<unsigned short> &dtmImage, OrthoImage<unsigned short> &dsm2Image, OrthoImage<unsigned long> &labelImage) {

    // Convert horizontal and vertical uncertainty values to bin units.
    int dh_bins = MAX(1, (int) floor(dh_meters / dsmImage.gsd));
    unsigned int dz_short = (unsigned int) (dz_meters / dsmImage.scale);
    unsigned int agl_short = (unsigned int) (agl_meters / dsmImage.scale);
    unsigned int maxTreeHeightScaled = (max_tree_height_meters / minImage.scale);
    unsigned short threshold = (dz_meters / dsmImage.scale);

    // Find many of the trees by comparing MIN and MAX. Set their values to void.
    shr3d::OrthoImage<unsigned short> varImage = dsmImage - minImage;

    // Copy DSM, and apply tree filter so any location where the last return image differs from the DSM
    // by more than THRESHOLD is set to void.
    dsm2Image = dsmImage;
    shr3d::Image<unsigned short>::filter(&dsm2Image, &varImage,
            [&](unsigned short* val, const unsigned short& ref, std::vector<unsigned short> &ngbrs) {
        // CAUTION: This is a hack to address an observed lidar sensor issue and may not generalize well.
        if (ref <= maxTreeHeightScaled) {
            // Set dsm to void if none of the neighbors are solid (var is < threshold)
            if (std::none_of(ngbrs.begin(), ngbrs.end(), [&](unsigned short v){ return v<=threshold; }))
                *val = 0;
        }
    }, 1, 0, false);


    // Generate label image.
    labelImage = OrthoImage<unsigned long>(&dsmImage);

    // Allocate a DTM image as SHORT and copy in the Min DSM values.
    dtmImage = minImage;

    // Classify ground points.
    classifyGround(labelImage, dsm2Image, dtmImage, dh_bins, dz_short);

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
    dtmImage.medianFilter(1, dz_short);

    // Refine the object label image and export building outlines.
    classifyNonGround(dsm2Image, dtmImage, labelImage, dz_short, agl_short, (float) min_area_meters);

    // Fill small voids in the DTM after all processing is complete.
    dtmImage.fillVoidsPyramid(true, 2);
}

OrthoImage<unsigned char> Shr3dder::labelClasses(
        const OrthoImage<unsigned short> &dsmImage, const OrthoImage<unsigned short> &dtmImage,
        const OrthoImage<unsigned short> &dsm2Image, const OrthoImage<unsigned long> &labelImage) {
    // Produce a classification raster image with LAS standard point classes.

    // Convert horizontal and vertical uncertainty values to bin units.
    int dh_bins = MAX(1, (int) floor(dh_meters / dsmImage.gsd));
    unsigned int maxTreeHeightScaled = (max_tree_height_meters / dtmImage.scale);
    unsigned int dz_short = (unsigned int) (dz_meters / dsmImage.scale);
    unsigned int agl_short = (unsigned int) (agl_meters / dsmImage.scale);

    OrthoImage<unsigned char> classImage(&dsmImage);
    for (unsigned int j = 0; j < classImage.height; j++) {
        for (unsigned int i = 0; i < classImage.width; i++) {
            // Set default as unlabeled.
            classImage.data[j][i] = LAS_UNCLASSIFIED;

            // Label trees.
            if ((dsm2Image.data[j][i] == 0) ||
                (fabs((float) dsm2Image.data[j][i] - (float) dtmImage.data[j][i]) > agl_short))
                classImage.data[j][i] = LAS_TREE;

            // Label buildings.
            if (labelImage.data[j][i] == 1 ||
                    (dsmImage.data[j][i] > dtmImage.data[j][i] + maxTreeHeightScaled))
                classImage.data[j][i] = LAS_BUILDING;

            // Label ground.
            if ((dsmImage.data[j][i] == 0) ||
                    fabs((float) dsm2Image.data[j][i] - (float) dtmImage.data[j][i]) < dz_short)
                classImage.data[j][i] = LAS_GROUND;
        }
    }

    // Fill in building labels on edge
    for (unsigned int i = 0; i < 5; ++i) {
        typedef std::pair<unsigned char, unsigned short> FType;
        shr3d::OrthoImage<unsigned char> classImageRef(classImage);
        shr3d::Image<unsigned char>::filter2(&classImage, &classImageRef, &dsmImage,
                [&](const FType& ref) { return ref.first == LAS_TREE; },
                [&](unsigned char* val, const FType& ref, std::vector<FType> &ngbrs) {
                    if (any_of(ngbrs.begin(),ngbrs.end(), [&](FType ngbr) {
                        return ngbr.first == LAS_BUILDING && (unsigned short) abs(ngbr.second-ref.second) < dz_short;
                    }))
                        *val = LAS_BUILDING;
                },dh_bins);
    }

    // Filter building class
    classImage.filter([&](unsigned char* val, const unsigned char& ref, std::vector<unsigned char> &ngbrs) {
        if ((ref != LAS_BUILDING) && ((size_t) std::count(ngbrs.begin(),ngbrs.end(),LAS_BUILDING) >= ngbrs.size()/2))
            *val = LAS_BUILDING;
    });

    // Fill missing labels inside building regions.
    fillInsideBuildings(classImage);

    return classImage;
}

OrthoImage<unsigned char> Shr3dder::labelBuildings(const OrthoImage<unsigned char> &classImage) {
    OrthoImage<unsigned char> bldgImage(classImage);
    for (unsigned int j = 0; j < classImage.height; j++) {
        for (unsigned int i = 0; i < classImage.width; i++) {
            if (classImage.data[j][i] != LAS_BUILDING) bldgImage.data[j][i] = 0;
        }
    }
    return bldgImage;
}


// Extend object boundaries to capture points missed around the edges.
void extendObjectBoundaries(OrthoImage<unsigned short> &dsmImage, OrthoImage<unsigned long> &labelImage,
                            int edgeResolution, unsigned int minDistanceShortValue) {
    // Loop enough to capture the edge resolution.
    for (int k = 0; k < edgeResolution; k++) {
        // First, label any close neighbor LABEL_TEMP.
        for (unsigned int j = 1; j < labelImage.height - 1; j++) {
            for (unsigned int i = 1; i < labelImage.width - 1; i++) {
                // For any labeled point, check all neighbors.
                if (labelImage.data[j][i] == 1) {
                    for (unsigned int jj = j - 1; jj <= j + 1; jj++) {
                        for (unsigned int ii = i - 1; ii <= i + 1; ii++) {
                            if (labelImage.data[jj][ii] == 1) continue;
                            if (((float) dsmImage.data[j][i] - (float) dsmImage.data[jj][ii]) <
                                minDistanceShortValue / 2.0) {
                                labelImage.data[jj][ii] = LABEL_TEMP;
                            }
                        }
                    }
                }
            }
        }

        // Then label any high points labeled LABEL_TEMP as an object of interest.
        for (unsigned int j = 0; j < labelImage.height; j++) {
            for (unsigned int i = 0; i < labelImage.width; i++) {
                if (labelImage.data[j][i] == LABEL_TEMP) {
                    // Check to make sure this point is also higher than one of its neighbors.
                    unsigned int j1 = j<1 ? 0 : j-1;
                    unsigned int j2 = MIN(j + 1, labelImage.height - 1);
                    unsigned int i1 = i<1 ? 0 : i-1;
                    unsigned int i2 = MIN(i + 1, labelImage.width - 1);
                    for (unsigned int jj = j1; jj <= j2; jj++) {
                        for (unsigned int ii = i1; ii <= i2; ii++) {
                            if (((float) dsmImage.data[j][i] - (float) dsmImage.data[jj][ii]) >
                                minDistanceShortValue / 2.0) {
                                labelImage.data[j][i] = 1;
                            }
                        }
                    }
                }
            }
        }
    }

    // Reset any temporary values.
    for (unsigned int j = 0; j < labelImage.height; j++) {
        for (unsigned int i = 0; i < labelImage.width; i++) {
            if (labelImage.data[j][i] == LABEL_TEMP) {
                labelImage.data[j][i] = LABEL_GROUND;
            }
        }
    }

}

// Label boundaries of objects above ground level.
void labelObjectBoundaries(OrthoImage<unsigned short> &dsmImage, OrthoImage<unsigned long> &labelImage,
                           int edgeResolution, unsigned int minDistanceShortValue) {
    // Initialize the labels to LABEL_GROUND.
    for (unsigned int j = 0; j < labelImage.height; j++) {
        for (unsigned int i = 0; i < labelImage.width; i++) {
            labelImage.data[j][i] = LABEL_GROUND;
        }
    }

    // Mark the label image with object boundaries.
    float threshold = (float) minDistanceShortValue;
    for (int j = 0; j < (int) labelImage.height; j++) {
        for (int i = 0; i < (int) labelImage.width; i++) {
            // Look for Z steps greater than a threshold.
            // Interestingly, this works about as well as checking every step.
            for (int dj = -edgeResolution; dj <= edgeResolution; dj += edgeResolution) {
                for (int di = -edgeResolution; di <= edgeResolution; di += edgeResolution) {
                    int j2 = MIN(MAX(0, j + dj), dsmImage.height - 1);
                    int i2 = MIN(MAX(0, i + di), dsmImage.width - 1);
                    if (dsmImage.data[j2][i2] != 0) {
                        // Remove local slope to avoid tagging rough terrain.
                        int j3 = MIN(MAX(0, j + dj * 2), dsmImage.height - 1);
                        int i3 = MIN(MAX(0, i + di * 2), dsmImage.width - 1);
                        float myGradient = (float) dsmImage.data[j][i] - (float) dsmImage.data[j2][i2];
                        float neighborGradient = (float) dsmImage.data[j2][i2] - (float) dsmImage.data[j3][i3];
                        float distance = (myGradient - neighborGradient);
                        if (distance > threshold) labelImage.data[j][i] = 1;
                    }
                }
            }
        }
    }
}

bool findObjectBoundsInColumn(OrthoImage<unsigned long> &labelImage, ObjectType &obj, unsigned int column, int &min_row, int &max_row) {
    min_row = -1;
    max_row = -1;

    int j1 = MAX(0, obj.ymin - 1);
    int j2 = MIN(obj.ymax + 1, labelImage.height - 1);

    // Find 1st row with matching label
    for (int j = j1; j <= j2; j++) {
        if (labelImage.data[j][column] == obj.label) {
            min_row = j;
            break;
        }
    }

    // If no labels in this row, then continue.
    if (min_row == -1)
        return false;

    // Find last row with matching label
    for (int j = j2; j >= min_row; j--) {
        if (labelImage.data[j][column] == obj.label) {
            max_row = j;
            break;
        }
    }
    return true;
}

// Fill inside the object countour labels if points are above the nearby ground level.
void fillObjectBounds(OrthoImage<unsigned long> &newLabelImage, OrthoImage<unsigned long> &labelImage, OrthoImage<unsigned short> &dsmImage, ObjectType &obj,
                      int edgeResolution) {
    unsigned int label = obj.label;

    // For the case when the object is on the left or right edge of the image, find which rows on the edge are within the image
    int lmin,lmax,rmin,rmax;
    bool ledge = (obj.xmin-1 <= 0) && findObjectBoundsInColumn(labelImage, obj, 0, lmin, lmax);
    bool redge = (obj.xmax+1 >= labelImage.width-1) && findObjectBoundsInColumn(labelImage, obj, labelImage.width-1, rmin, rmax);

    // Loop on rows, filling in labels.
    for (unsigned int j = MAX(0, obj.ymin - 1); j <= MIN(obj.ymax + 1, labelImage.height - 1); j++) {
        // Get start index.
        int startIndex = -1;
        if (ledge && (int) j >= lmin && (int) j <= lmax) {
            startIndex = 0; // If the object is on the left edge and this row is within the bounds, start the index at the edge
        } else {
            for (unsigned int i = MAX(0, obj.xmin - 1); i <= MIN(obj.xmax + 1, labelImage.width - 1); i++) {
                if (labelImage.data[j][i] == label) {
                    startIndex = i;
                    break;
                }
            }
        }

        // If no labels in this row, then continue.
        if (startIndex == -1) continue;

        // Get stop index.
        int stopIndex = -1;
        if (redge && (int) j >= rmin && (int) j <= rmax) {
            stopIndex = labelImage.width - 1; // If the object is on the right edge and this row is within the bounds, have the index stop at the edge
        } else {
            for (int i = MIN(obj.xmax + 1, labelImage.width - 1); i >= MAX(0, obj.xmin - 1); i--) {
                if (labelImage.data[j][i] == label) {
                    stopIndex = i;
                    break;
                }
            }
        }

		// If entire column is labeled, then continue.
		if ((startIndex == 0) && (stopIndex == (int) labelImage.width-1)) continue;

        // Get max ground level height for this row.
        // If the DSM image value is void, then the ground level value is zero, so that's ok.
        unsigned short groundLevel;
        if (startIndex == 0)
            groundLevel = dsmImage.data[j][stopIndex + 1];
        else if (stopIndex == (int) labelImage.width - 1)
            groundLevel = dsmImage.data[j][startIndex - 1];
        else
            groundLevel = MAX(dsmImage.data[j][startIndex - 1], dsmImage.data[j][stopIndex + 1]);

        // Fill in the label for any point in between that has height above ground level.
        for (unsigned int i = startIndex; (int) i <= stopIndex; i++) {
            if (dsmImage.data[j][i] > groundLevel) {
                if (labelImage.data[j][i] != label && newLabelImage.data[j][i] != LABEL_OBJECT) newLabelImage.data[j][i] = LABEL_IN_ONE;
            } else {
                if (labelImage.data[j][i] == label && newLabelImage.data[j][i] != LABEL_OBJECT) newLabelImage.data[j][i] = LABEL_GROUND;
            }
        }
    }

    // Loop on columns, filling in labels.
    for (unsigned int i = MAX(0, obj.xmin - 1); i <= MIN(obj.xmax + 1, labelImage.width - 1); i++) {
        // Get start & stop indices
        int startIndex, stopIndex;
        if (!findObjectBoundsInColumn(labelImage, obj, i, startIndex, stopIndex))
            continue;

		// If entire row is labeled, then continue.
		if ((startIndex == 0) && (stopIndex == (int) labelImage.height-1)) continue;

        // Get max ground level height for this row.
        unsigned short groundLevel;
        if (startIndex == 0)
            groundLevel = dsmImage.data[stopIndex + 1][i];
        else if (stopIndex == (int) labelImage.height - 1)
            groundLevel = dsmImage.data[startIndex - 1][i];
        else
            groundLevel = MAX(dsmImage.data[startIndex - 1][i], dsmImage.data[stopIndex + 1][i]);

        // Fill in the label for any point in between that has height above ground level.
        // This time make sure both the horizontal and vertical check pass and set to LABEL_ACCEPTED.
        for (unsigned int j = startIndex; (int) j <= stopIndex; j++) {
            if (dsmImage.data[j][i] > groundLevel) {
                if (((labelImage.data[j][i] == label) || (newLabelImage.data[j][i] == LABEL_IN_ONE)) && newLabelImage.data[j][i] != LABEL_OBJECT) {
                    newLabelImage.data[j][i] = LABEL_ACCEPTED;
                }
            }
        }
    }

    // Erode the labels with a kernel size based on edge resolution.
    int rad = edgeResolution;
    for (int j = MAX(0, obj.ymin - 1); j <= MIN(obj.ymax + 1, (int) labelImage.height - 1); j++) {
        for (int i = MAX(0, obj.xmin - 1); i <= MIN(obj.xmax + 1, (int) labelImage.width - 1); i++) {
            if (newLabelImage.data[j][i] == LABEL_ACCEPTED) {
                int i1 = MAX(i - rad, 0);
                int i2 = MIN(i + rad, labelImage.width - 1);
                int j1 = MAX(j - rad, 0);
                int j2 = MIN(j + rad, labelImage.height - 1);
                for (int jj = j1; jj <= j2; jj++) {
                    for (int ii = i1; ii <= i2; ii++) {
                        if (newLabelImage.data[jj][ii] != LABEL_ACCEPTED && newLabelImage.data[j][i] != LABEL_OBJECT) newLabelImage.data[jj][ii] = LABEL_TEMP;
                    }
                }
            }
        }
    }

    // Update object bounds to include erosion.
    obj.xmin = MAX(0, obj.xmin - edgeResolution - 1);
    obj.ymin = MAX(0, obj.ymin - edgeResolution - 1);
    obj.xmax = MIN(obj.xmax + edgeResolution + 1, labelImage.width - 1);
    obj.ymax = MIN(obj.ymax + edgeResolution + 1, labelImage.height - 1);

    for (unsigned int j = obj.ymin; j <= obj.ymax; j++) {
        for (unsigned int i = obj.xmin; i <= obj.xmax; i++) {
            if (newLabelImage.data[j][i] == LABEL_TEMP) {
                newLabelImage.data[j][i] = LABEL_ACCEPTED;
            }
        }
    }

    // Finish up the labels.
    for (unsigned int j = MAX(0, obj.ymin - 1); j <= MIN(obj.ymax + 1, labelImage.height - 1); j++) {
        for (unsigned int i = MAX(0, obj.xmin - 1); i <= MIN(obj.xmax + 1, labelImage.width - 1); i++) {
            if (newLabelImage.data[j][i] == label) newLabelImage.data[j][i] = LABEL_GROUND;
            if (newLabelImage.data[j][i] == LABEL_ACCEPTED) newLabelImage.data[j][i] = LABEL_OBJECT;
        }
    }
}

// Add neighboring pixels to an object.
bool addNeighbors(std::vector<PixelType> &neighbors, OrthoImage<unsigned long> &labelImage,
                  OrthoImage<unsigned short> &dsmImage, ObjectType &obj, unsigned int dzShort) {
    // Get neighbors for all pixels in the list.
    std::vector<PixelType> newNeighbors;
    for (size_t k = 0; k < neighbors.size(); k++) {
        unsigned int i = neighbors[k].i;
        unsigned int j = neighbors[k].j;
        for (unsigned int jj = j<1 ? 0 : j-1; jj <= MIN(j + 1, labelImage.height - 1); jj++) {
            for (unsigned int ii = i<1 ? 0 : i-1; ii <= MIN(i + 1, labelImage.width - 1); ii++) {
                // Skip if pixel is already labeled or if it is LABEL_GROUND.
                // Note that non-ground labels are initialized with 1.
                if (labelImage.data[jj][ii] > 1) continue;

                // Skip if height is too different.
                if (fabs((float) dsmImage.data[jj][ii] - (float) dsmImage.data[j][i]) > (float) dzShort) continue;

                // Update the label.
                labelImage.data[jj][ii] = labelImage.data[j][i];

                // Add to the new list.
                PixelType pixel = {ii, jj};
                newNeighbors.push_back(pixel);

                // Update object bounds.
                obj.xmin = MIN(obj.xmin, ii);
                obj.xmax = MAX(obj.xmax, ii);
                obj.ymin = MIN(obj.ymin, jj);
                obj.ymax = MAX(obj.ymax, jj);
                obj.count++;
            }
        }
    }

    // Replace the neighbors list with the new one. No need to merge them.
    if (newNeighbors.size() > 0) {
        neighbors = newNeighbors;
        return true;
    } else return false;
}

// Group connected labeled pixels into objects.
void groupObjects(OrthoImage<unsigned long> &labelImage, OrthoImage<unsigned short> &dsmImage,
                  std::vector<ObjectType> &objects, long maxCount, unsigned int dzShort) {
    // Sweep from top left to bottom right, assigning object labels.
    long maxGroupSize = 0;
    unsigned int label = 1;
    int bigCount = 0;
    for (unsigned int j = 0; j < labelImage.height; j++) {
        for (unsigned int i = 0; i < labelImage.width; i++) {
            // Skip unlabeled pixels.
            if (labelImage.data[j][i] == LABEL_GROUND) continue;

            // Skip already labeled pixels.
            if (labelImage.data[j][i] > 1) continue;

            // Create a new label.
            label++;

            // Create a new object structure.
            ObjectType obj;

            // Initialize the object and neighbors structures.
            obj.label = label;
            obj.xmin = i;
            obj.ymin = j;
            obj.xmax = i;
            obj.ymax = j;
            obj.count = 1;
            std::vector<PixelType> neighbors;
            PixelType pixel = {i, j};
            neighbors.push_back(pixel);
            labelImage.data[j][i] = label;

            // Get points in new group.
            bool keepSearching = true;
            while (keepSearching) {
                keepSearching = addNeighbors(neighbors, labelImage, dsmImage, obj, dzShort);

                // This is very quick but not especially smart.
                // Try to find a reasonably quick way to split the regions more sensibly.
                // That said, this also seems to work well enough.
                if (obj.count > maxCount) {
                    bigCount++;
                    break;
                }
            }

            // Add this object to the list.
            objects.push_back(obj);
            maxGroupSize = MAX(maxGroupSize, obj.count);
        }
    }
    printf("Max group size = %ld\n", maxGroupSize);
    printf("Number of cropped groups = %d\n", bigCount);
}

// Finish label image for display as image overlay in QT Reader.
// Set all labeled values to 1. Leave all unlabeled values LABEL_GROUND.
void finishLabelImage(OrthoImage<unsigned long> &labelImage) {
    for (unsigned int j = 0; j < labelImage.height; j++) {
        for (unsigned int i = 0; i < labelImage.width; i++) {
            if ((labelImage.data[j][i] != LABEL_GROUND)) {
                labelImage.data[j][i] = 1;
            }
        }
    }
}

// Classify ground points, fill the voids, and generate a bare earth terrain model. 
void Shr3dder::classifyGround(OrthoImage<unsigned long> &labelImage, OrthoImage<unsigned short> &dsmImage,
                              OrthoImage<unsigned short> &dtmImage, int dhBins, unsigned int dzShort) {
    // Fill voids.
    printf("Filling voids...\n");
    dtmImage.fillVoidsPyramid(true);

    // Allocate a binary label image to indicate voids to be filled.
    // The long integer label image has unique labels for objects detected in each iteration.
    OrthoImage<unsigned char> voidImage(&labelImage);

    // Iteratively label and remove objects from the DEM.
    // Each new iteration removes debris not identified by the previous iteration.
    unsigned int numIterations = 5;
    long maxCount = 10000 / (dsmImage.gsd * dsmImage.gsd);    // max count is in meters
    for (unsigned int k = 0; k < numIterations; k++) {
        printf("Iteration #%d\n", k + 1);

        // Label the object boundaries.
        printf("Labeling object boundaries...\n");
        labelObjectBoundaries(dtmImage, labelImage, dhBins, dzShort);

        // Extend labels for object boundaries..
        printf("Extending object boundaries...\n");
        extendObjectBoundaries(dtmImage, labelImage, dhBins, dzShort);

        // Group the objects.
        printf("Grouping objects...\n");
        std::vector<ObjectType> objects;
        groupObjects(labelImage, dtmImage, objects, maxCount, dzShort);
        printf("Number of objects = %ld\n", objects.size());

        // Generate object groups and void fill them in the DEM image.
        printf("Labeling and removing objects...\n");
        OrthoImage<unsigned long> newlabelImage(labelImage);
        for (size_t i = 0; i < objects.size(); i++) {
            fillObjectBounds(newlabelImage,labelImage, dtmImage, objects[i], dhBins);
        }
        swap(newlabelImage,labelImage);

        // Update the label image values for easy viewing.
        printf("Finishing label image for display...\n");
        finishLabelImage(labelImage);

        // Update the void image.
        printf("Updating void image...\n");
        for (unsigned int j = 0; j < voidImage.height; j++) {
            for (unsigned int i = 0; i < voidImage.width; i++) {
                if (labelImage.data[j][i] == 1) voidImage.data[j][i] = 1;
            }
        }

        // Fill voids.
        for (unsigned int j = 0; j < labelImage.height; j++) {
            for (unsigned int i = 0; i < labelImage.width; i++) {
                if (voidImage.data[j][i] == 1) dtmImage.data[j][i] = 0;
            }
        }
        bool noSmoothing = true;
        if (k == numIterations - 1) noSmoothing = false;
        printf("Filling voids with noSmoothing = %d\n", noSmoothing);
        dtmImage.fillVoidsPyramid(noSmoothing);
    }

    // If any DTM points are above the DSM, then restore the DSM values.
    for (unsigned int j = 0; j < dtmImage.height; j++) {
        for (unsigned int i = 0; i < dtmImage.width; i++) {
            if (dtmImage.data[j][i] >= dsmImage.data[j][i] && dsmImage.data[j][i] != 0) {
                dtmImage.data[j][i] = dsmImage.data[j][i];
                labelImage.data[j][i] = LABEL_GROUND;
                voidImage.data[j][i] = 0;
            }
        }
    }

    // Remove any leftover single point spikes.
    printf("Removing spikes...\n");
    for (int j = 0; j < (int) dtmImage.height; j++) {
        for (int i = 0; i < (int) dtmImage.width; i++) {
            float minDiff = FLT_MAX;
            for (int jj = -1; jj <= 1; jj++) {
                int j2 = MAX(0, MIN((int) dsmImage.height - 1, j + jj));
                for (int ii = -1; ii <= 1; ii++) {
                    if ((ii == 0) && (jj == 0)) continue;
                    int i2 = MAX(0, MIN((int) dsmImage.width - 1, i + ii));
                    float diff = MAX(0, (float) dtmImage.data[j][i] - (float) dtmImage.data[j2][i2]);
                    minDiff = MIN(minDiff, diff);
                }
            }
            if (minDiff > dzShort / 2.0) {
                labelImage.data[j][i] = 1;
                voidImage.data[j][i] = 1;
                dtmImage.data[j][i] = 0.0;
            }
        }
    }

    // Fill voids.
    printf("Filling voids...\n");
    for (unsigned int j = 0; j < labelImage.height; j++) {
        for (unsigned int i = 0; i < labelImage.width; i++) {
            if (voidImage.data[j][i] == 1) dtmImage.data[j][i] = 0;
        }
    }
    dtmImage.fillVoidsPyramid(false);

    // Mark all voids.
    printf("Marking voids in label image after all iterations are complete...\n");
    for (unsigned int j = 0; j < voidImage.height; j++) {
        for (unsigned int i = 0; i < voidImage.width; i++) {
            if (voidImage.data[j][i] == 1) labelImage.data[j][i] = 1;
            else labelImage.data[j][i] = LABEL_GROUND;
        }
    }
}

// Classify non-ground points.
void Shr3dder::classifyNonGround(OrthoImage<unsigned short> &dsmImage, OrthoImage<unsigned short> &dtmImage,
                                 OrthoImage<unsigned long> &labelImage, unsigned int dzShort, unsigned int aglShort,
                                 float minAreaMeters) {
    // Compute minimum number of points based on threshold given for area.
    // Note that ISPRS challenges indicate that performance is dramatically better for structures larger than 50m area.
    int minPointCount = int(minAreaMeters / (dsmImage.gsd * dsmImage.gsd));
    printf("Min points for removing small objects = %d\n", minPointCount);

    // Apply AGL threshold to individual points to reduce clutter.
    // Note that ground level clutter tends to be less than 2m AGL.
    for (unsigned int j = 0; j < labelImage.height; j++) {
        for (unsigned int i = 0; i < labelImage.width; i++) {
            if (labelImage.data[j][i] != LABEL_GROUND) {
                if (dsmImage.data[j][i] == 0) labelImage.data[j][i] = LABEL_GROUND;
                else if (((float) dsmImage.data[j][i] - (float) dtmImage.data[j][i]) < (float) aglShort)
                    labelImage.data[j][i] = LABEL_GROUND;
            }
        }
    }

    // Group the labeled objects based on height similarity.
    {
        std::vector<ObjectType> objects;
        groupObjects(labelImage, dsmImage, objects, LONG_MAX, dzShort / 2);

        // Accept or reject each object independently based on boundary values and flatness.
        int numRejected = 0;
        for (unsigned int k = 0; k < objects.size(); k++) {
            bool reject = false;

            // Reject if boundary gradients with neighbors labeled ground are too small.
            float meanGradient = 0.0;
            int count = 0;
            for (int j = objects[k].ymin; j <= objects[k].ymax; j++) {
                for (int i = objects[k].xmin; i <= objects[k].xmax; i++) {
                    if (labelImage.data[j][i] == objects[k].label) {
                        for (int jj = -1; jj <= 1; jj++) {
                            unsigned int j2 = MAX(0, MIN((int) dsmImage.height - 1, j + jj));
                            for (int ii = -1; ii <= 1; ii++) {
                                unsigned int i2 = MAX(0, MIN((int) dsmImage.width - 1, i + ii));
                                //if (labelImage.data[j2][i2] != objects[k].label)
                                if (labelImage.data[j2][i2] == LABEL_GROUND) {
                                    unsigned int j3 = MAX(0, MIN((int) dsmImage.height - 1, j + jj * 2));
                                    unsigned int i3 = MAX(0, MIN((int) dsmImage.width - 1, i + ii * 2));

                                    // These assume I'm higher than my neighbors.
                                    float myGradient = MAX(0, ((float) dsmImage.data[j][i] -
                                                               (float) dsmImage.data[j2][i2]));
                                    float neighborGradient = MAX(0, ((float) dsmImage.data[j2][i2] -
                                                                     (float) dsmImage.data[j3][i3]));
                                    meanGradient += MAX(0, (myGradient - neighborGradient));
                                    count++;
                                }
                            }
                        }
                    }
                }
            }
            meanGradient /= count;
            if ((meanGradient < dzShort / 2.0) && (meanGradient != 0.0)) reject = true;

            // If rejected, then label each point.
            if (reject) {
                numRejected++;
                for (unsigned int j = objects[k].ymin; j <= objects[k].ymax; j++) {
                    for (unsigned int i = objects[k].xmin; i <= objects[k].xmax; i++) {
                        if (labelImage.data[j][i] == objects[k].label) {
                            labelImage.data[j][i] = LABEL_GROUND;
                        }
                    }
                }
            }
        }
    }

    // Erode and then dilate labels to remove narrow objects.
    {
        OrthoImage<unsigned long> tempImage(&labelImage);
        for (unsigned int j = 0; j < labelImage.height; j++) {
            for (unsigned int i = 0; i < labelImage.width; i++) {
                tempImage.data[j][i] = labelImage.data[j][i];
            }
        }
        for (int j = 0; j < (int) labelImage.height; j++) {
            for (int i = 0; i < (int) labelImage.width; i++) {
                if (labelImage.data[j][i] != LABEL_GROUND) {
                    // Define neighbor bounds;
                    int i1 = MAX(0, i - 1);
                    int i2 = MIN(i + 1, labelImage.width - 1);
                    int j1 = MAX(0, j - 1);
                    int j2 = MIN(j + 1, labelImage.height - 1);

                    // Unlabel any point with an unlabeled neighbor.
                    for (int jj = j1; jj <= j2; jj++) {
                        for (int ii = i1; ii <= i2; ii++) {
                            if (labelImage.data[jj][ii] == LABEL_GROUND) {
                                tempImage.data[j][i] = LABEL_GROUND;
                            }
                        }
                    }
                }
            }
        }
        for (int j = 0; j < (int) labelImage.height; j++) {
            for (int i = 0; i < (int) labelImage.width; i++) {
                if (labelImage.data[j][i] != LABEL_GROUND) {
                    // Define neighbor bounds;
                    int i1 = MAX(0, i - 1);
                    int i2 = MIN(i + 1, labelImage.width - 1);
                    int j1 = MAX(0, j - 1);
                    int j2 = MIN(j + 1, labelImage.height - 1);

                    // Unlabel any point with no labeled neighbors after erosion.
                    bool found = false;
                    for (int jj = j1; jj <= j2; jj++) {
                        for (int ii = i1; ii <= i2; ii++) {
                            if (tempImage.data[jj][ii] != LABEL_GROUND) found = true;
                        }
                    }
                    if (!found) labelImage.data[j][i] = LABEL_GROUND;
                }
            }
        }
    }

    // Reset all non-ground labels to one.
    for (unsigned int j = 0; j < labelImage.height; j++) {
        for (unsigned int i = 0; i < labelImage.width; i++) {
            if (labelImage.data[j][i] != LABEL_GROUND) labelImage.data[j][i] = 1;
        }
    }

    // Group labeled points to remove small objects.
    // Do not split groups based on height or point count.
    {
        printf("Grouping to remove small objects...\n");
        std::vector<ObjectType> objects;
        groupObjects(labelImage, dsmImage, objects, LONG_MAX, INT_MAX);
        int numRejected = 0;
        for (unsigned int k = 0; k < objects.size(); k++) {
            bool reject = false;

            // Reject if too small.
            if (objects[k].count < minPointCount) reject = true;

            // If rejected, then label each point.
            if (reject) {
                numRejected++;
                for (unsigned int j = objects[k].ymin; j <= objects[k].ymax; j++) {
                    for (unsigned int i = objects[k].xmin; i <= objects[k].xmax; i++) {
                        if (labelImage.data[j][i] == objects[k].label) {
                            labelImage.data[j][i] = LABEL_GROUND;
                        }
                    }
                }
            }
        }
        printf("Number of small objects rejected = %d\n", numRejected);
    }

    // Reset all non-ground labels to one.
    for (unsigned int j = 0; j < labelImage.height; j++) {
        for (unsigned int i = 0; i < labelImage.width; i++) {
            if (labelImage.data[j][i] != LABEL_GROUND) labelImage.data[j][i] = 1;
        }
    }
}

// Add neighboring pixels to a void group.
bool
addClassNeighbors(std::vector<PixelType> &neighbors, OrthoImage<unsigned char> &classImage,
                  OrthoImage<unsigned char> &labeled, unsigned int label) {
    // Get neighbors for all pixels in the list.
    std::vector<PixelType> newNeighbors;
    for (size_t k = 0; k < neighbors.size(); k++) {
        unsigned int i = neighbors[k].i;
        unsigned int j = neighbors[k].j;
        for (unsigned int jj = MAX(0, (int) j - 1); jj <= MIN(j + 1, classImage.height - 1); jj++) {
            for (unsigned int ii = MAX(0, (int) i - 1); ii <= MIN(i + 1, classImage.width - 1); ii++) {
                // If already labeled, then skip.
                if (labeled.data[jj][ii] == 1) continue;

                // If not the same label, then skip.
                if (classImage.data[jj][ii] != label) continue;

                // Update the label.
                labeled.data[jj][ii] = 1;

                // Add to the new list.
                PixelType pixel = {ii, jj};
                newNeighbors.push_back(pixel);
            }
        }
    }

    // Update the neighbors list.
    if (newNeighbors.size() > 0) {
        for (size_t k = 0; k < newNeighbors.size(); k++) {
            neighbors.push_back(newNeighbors[k]);
        }
//		neighbors = newNeighbors;
        return true;
    } else return false;
}

// Fill in any pixels labeled tree that fall entirely within a labeled building group.
void Shr3dder::fillInsideBuildings(OrthoImage<unsigned char> &classImage) {
    int numFilled = 0;
    OrthoImage<unsigned char> labeled(&classImage);
    for (unsigned int j = 0; j < classImage.height; j++) {
        for (unsigned int i = 0; i < classImage.width; i++) {
            bool consider = ((labeled.data[j][i] == 0) && (classImage.data[j][i] == LAS_TREE));
            if (!consider) continue;

            // Get all pixels in this contiguous group.
            std::vector<PixelType> neighbors;
            PixelType pixel = {i, j};
            neighbors.push_back(pixel);
            bool keepSearching = true;
            unsigned int label = (unsigned int) classImage.data[j][i];
            labeled.data[j][i] = 1;
            while (keepSearching) {
                keepSearching = addClassNeighbors(neighbors, classImage, labeled, label);
            }

            // Check all pixels for neighboring labels.
            bool inside = true;
            for (size_t k = 0; k < neighbors.size(); k++) {
                long i1 = neighbors[k].i;
                long j1 = neighbors[k].j;
                for (long jj = MAX(0, j1 - 1); jj <= MIN(j1 + 1, classImage.height - 1); jj++) {
                    for (long ii = MAX(0, i1 - 1); ii <= MIN(i1 + 1, classImage.width - 1); ii++) {
                        if ((labeled.data[jj][ii] == 0) && (classImage.data[jj][ii] != LAS_BUILDING)) {
                            inside = false;
                        }
                    }
                }
            }

            // If the region is completely inside a building region, then fill it.
            if (inside) {
                for (size_t k = 0; k < neighbors.size(); k++) {
                    numFilled++;
                    classImage.data[neighbors[k].j][neighbors[k].i] = LAS_BUILDING;
                }
            }

        }
    }
    printf("Removed %d tree pixels inside building label groups.\n", numFilled);
}
