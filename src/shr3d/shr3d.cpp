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

#include "shr3d.h"

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <climits>
#include <vector>

using namespace shr3d;

void Shr3dder::process(std::map<ImageType,std::string> outputFilenames) {
    for (std::pair<ImageType,std::string> elem : outputFilenames) {
        if (elem.second.empty())
            continue;
        auto fname = elem.second.c_str();
        bool b = true;
        switch (elem.first) {
            case DSM:
                b = getDSM().write(fname, true, egm96); break;
            case MIN:
                b = getMIN().write(fname, true, egm96); break;
            case DTM:
                b = getDTM().write(fname, true, egm96); break;
            case INTENSITY:
                b = getINT().write(fname, false);   break;
            case CLASS:
                b = getCLS().write(fname, false);   break;
            case BUILDING:
                b = getBLDG().write(fname, false);  break;
            case BUILDING_OUTLINES:
                b = GeoPolygon<double>::write(fname, getBLDGPOLY()); break;
            case DSM2:
                b = getDSM2().write(fname, true, egm96);    break;
            case MINAGL:
                b = getMINAGL().write(fname, true, egm96);  break;
            case DTM0:
                b = getDTM0().write(fname, true, egm96);  break;
            case LABEL0:
                b = getLBL0().write(fname, false);  break;
            case LABEL:
                b = getLBL().write(fname, false);   break;
            case LABELED_BUILDINGS:
                b = getBLDGLBL().write(fname, false);   break;
            case LABELED_BUILDINGS_3:
                b = getBLDGLBL3().write(fname, false);  break;
            default:
                throw pdal::pdal_error("Unknown output image type\n");
        }
        if (!b)
            std::cerr << "Failed to write " << elem.second << std::endl;
    }
}

void sizeImageFromBox(OrthoImage<unsigned short> &im, pdal::BOX2D box, float dh_meters) {
    if (!box.empty()) {
        im.easting = box.minx;
        im.northing = box.miny;
        im.width = (unsigned int) ((box.maxx - box.minx) / dh_meters);
        im.height = (unsigned int) ((box.maxy - box.miny) / dh_meters);
        im.Allocate(im.width, im.height);
    }
}

void Shr3dder::createDSM() {
    printf("Creating DSM\n");
    if (is_pset_empty())
        throw pdal::pdal_error("Point cloud undefined / empty\n");
    sizeImageFromBox(dsmImage, bounds.to2d(), dh_meters);
    if (!dsmImage.readFromPointCloud(pset, (float) dh_meters, shr3d::MAX_VALUE))
        throw pdal::pdal_error("Error filling DSM from point cloud\n");

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
}

void Shr3dder::createMIN() {
    printf("Creating MIN\n");
    if (!is_pset_empty()) {
        sizeImageFromBox(minImage, bounds.to2d(), dh_meters);
        // Now get the minimum Z values for the DTM.
        if (!minImage.readFromPointCloud(pset, (float) dh_meters, shr3d::MIN_VALUE))
            throw pdal::pdal_error("Error filling MIN from point cloud\n");

        // Median filter, replacing only points differing by more than the AGL threshold.
        minImage.quantileFilter(2, (unsigned int) (agl_meters / minImage.scale), 0.33);

    } else {
        // Copy DSM
        minImage = getDSM();

        // Min filter, replacing only points differing by more than the AGL threshold.
        minImage.minFilter(4, (unsigned int) (agl_meters / minImage.scale));
    }

    // Fill small voids in the DSM.
    minImage.fillVoidsPyramid(true, 2);
}

void Shr3dder::createDSM2() {
    const OrthoImage<unsigned short> &dsmImage = getDSM();
    const OrthoImage<unsigned short> &minImage = getMIN();
    printf("Creating DSM2\n");

    // Find many of the trees by comparing MIN and MAX. Set their values to void.
    shr3d::OrthoImage<unsigned short> varImage = dsmImage - minImage;

    // Convert threshold & tree height to scaled units.
    unsigned int maxTreeHeightScaled = (max_tree_height_meters / minImage.scale);
    unsigned short threshold = (dz_meters / dsmImage.scale);

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
}

void adjust_scale(OrthoImage<unsigned short> &im, float new_offset, float new_scale) {
    if (new_offset == im.offset && new_scale == im.scale)
        return;
    for (std::vector<unsigned short> &r : im.data) {
        for (unsigned short &z : r) {
            if (z) { // Preserve voids
                z = (z*im.scale+im.offset-new_offset)/new_scale;
            }
        }
    }
    im.offset = new_offset;
    im.scale = new_scale;
}

bool Shr3dder::setDTM0(std::string dtmFile) {
    if (!dtm0Image.read(dtmFile.c_str()))
        return false;
    
    // If a pointcloud was specified, update bounds and gsd based on the DTM
    if (!is_pset_empty()) {
        dh_meters = dtm0Image.gsd;
        bounds = pdal::Bounds(pdal::BOX2D(dtm0Image.easting,dtm0Image.northing,
                dtm0Image.width*dh_meters+dtm0Image.easting,
                dtm0Image.height*dh_meters+dtm0Image.northing));
    }
    
    // Get DSM && MIN (get MIN since it could be created from base point cloud)
    // Note while I'm calling getDSM and getMIN, I'm using the class members (as they're not const)
    getDSM();
    getMIN();
    
    // Check that dsm and dtm match in everything important
    if (dtm0Image.zone != dsmImage.zone) {
        printf("ERROR: DSM/PSET zone %d doesn't match DTM zone %d\n",dsmImage.zone,dtm0Image.zone);
        return false;
    } else if (dtm0Image.gsd != dsmImage.gsd) {
        printf("ERROR: DSM gsd %f doesn't match DTM gsd %f\n",dsmImage.gsd,dtm0Image.gsd);
        return false;
    } else if (dtm0Image.easting != dsmImage.easting ||
            dtm0Image.northing != dsmImage.northing ||
            dtm0Image.width != dsmImage.width ||
            dtm0Image.height != dsmImage.height) {
        printf("ERROR: DSM bounds don't match DTM bounds\n");
        return false;
    }
    
    // Need to match scales and offsets between DSM and DTM
    float maxImageVal = (float) (pow(2.0, int(sizeof(unsigned short) * 8)) - 1);
    float offset = std::min(dtm0Image.offset,dsmImage.offset);
    float maxz = std::max(dtm0Image.offset+maxImageVal*dtm0Image.scale,
            dsmImage.offset+maxImageVal*dsmImage.scale);
    float scale = (maxz-offset) / maxImageVal;
    
    // Fix images
    adjust_scale(dsmImage,offset,scale);
    adjust_scale(minImage,offset,scale);
    adjust_scale(dtm0Image,offset,scale);
    
    // Adjust minimum image so that it's no lower than the DTM
    for (unsigned int j = 0; j < minImage.height; j++) {
        for (unsigned int i = 0; i < minImage.width; i++) {
            if (minImage.data[j][i] && dtm0Image.data[j][i] && 
                    (minImage.data[j][i] < dtm0Image.data[j][i]))
                minImage.data[j][i] = dtm0Image.data[j][i];
        }
    }
    
    // Need to initialize label0Image
    label0Image = OrthoImage<unsigned int>(&dtm0Image);
    
    return true;
}

void Shr3dder::createDTM0() {
    const OrthoImage<unsigned short> &dsmImage = getDSM();
    const OrthoImage<unsigned short> &minImage = getMIN();
    const OrthoImage<unsigned short> &dsm2Image = getDSM2();
    printf("Creating DTM0\n");

    // Generate label image.
    label0Image = OrthoImage<unsigned int>(&dsmImage);

    if (gnd_label >= 0) {
        if (is_pset_empty())
            throw pdal::pdal_error("Point cloud undefined / empty\n");
        PointCloud ground = pset.CropToClass(gnd_label);
        dtm0Image = OrthoImage<unsigned short>(&dsmImage); // Size to match dsm
        dtm0Image.readFromPointCloud(ground,dsmImage.gsd,shr3d::MAX_VALUE); // Using max to better match DSM

        // Need to match scales and offsets between DSM and DTM
        adjust_scale(dtm0Image,dsmImage.offset,dsmImage.scale);
    
        // Fill voids
        dtm0Image.fillVoidsPyramid(false);
    
        // Replicate DSM voids in DTM
        for (unsigned int j = 0; j < dtm0Image.height; j++) {
            for (unsigned int i = 0; i < dtm0Image.width; i++) {
                if (!dsmImage.data[j][i])
                    dtm0Image.data[j][i] = 0;
            }
        }

        return;
    }

    // Convert horizontal and vertical uncertainty values to bin units.
    int dh_bins = MAX(1, (int) floor(dh_meters / dsmImage.gsd));
    unsigned int dz_short = (unsigned int) (dz_meters / dsmImage.scale);

    // Allocate a DTM image as SHORT and copy in the Min DSM values.
    dtm0Image = minImage;

    // Classify ground points.
    classifyGround(label0Image, dsm2Image, dtm0Image, dh_bins, dz_short); // Updates labelImage, dtmImage

    // For DSM voids, also set DTM value to void.
    // Note: because we've changed the DSM by this point (setting voids where all the trees are),
    //  use the minImage which will have the same voids as the original DSM
    printf("Setting DTM0 values to VOID where DSM is VOID...\n");
    for (unsigned int j = 0; j < minImage.height; j++) {
        for (unsigned int i = 0; i < minImage.width; i++) {
            if (minImage.data[j][i] == 0) dtm0Image.data[j][i] = 0;
        }
    }

    // Median filter, replacing only points differing by more than the DZ threshold.
    dtm0Image.medianFilter(1, dz_short);
}

void Shr3dder::createDTM() {
    const OrthoImage<unsigned short> &dsm2Image = getDSM2();
    const OrthoImage<unsigned short> &dtm0Image = getDTM0();
    const OrthoImage<unsigned int> &label0Image = getLBL0();
    printf("Creating DTM\n");
    
    // Convert horizontal and vertical uncertainty values to bin units.
    unsigned int dz_short = (unsigned int) (dz_meters / dsmImage.scale);
    unsigned int agl_short = (unsigned int) (agl_meters / dsmImage.scale);
    
    // Allocate a DTM image as SHORT and copy in the Pre-DTM values.
    dtmImage = dtm0Image;

    // Allocate a LBL image copy in the Pre-Label values.
    labelImage = label0Image;

    // Refine the object label image and export building outlines.
    classifyNonGround(dsm2Image, dtmImage, labelImage, dz_short, agl_short, (float) min_area_meters); // Updates labelImage

    // Fill small voids in the DTM after all processing is complete.
    dtmImage.fillVoidsPyramid(true, 2);
}

void Shr3dder::createIntensity() {
    const OrthoImage<unsigned short> &dsmImage = getDSM();
    printf("Creating INT\n");
    if (is_pset_empty())
        throw pdal::pdal_error("Point cloud undefined / empty\n");

    intImage = OrthoImage<unsigned short>(&dsmImage);

    if (!is_pset_empty())    {
        OrthoImage<double> sumImage = OrthoImage<double>(&dsmImage);
        OrthoImage<unsigned short> cntImage = OrthoImage<unsigned short>(&dsmImage);

        // Determine scalings
        double mx = 1.0/dsmImage.gsd;
        double bx = mx*(pset.xOff - dsmImage.easting) - 0.5;
        double my =-1.0/dsmImage.gsd;
        double by = dsmImage.height-1 + my*(pset.yOff - dsmImage.northing) + 0.5;
        double mz = 1/dsmImage.scale;
        double bz = mz*(pset.zOff - dsmImage.offset);
        int dz_short = (unsigned int) (dz_meters / dsmImage.scale);

        // Loop over points
        double max_i = 0;
        for (unsigned long i = 0; i < pset.numPoints; i++) {
            int x = mx*pset.x(i) + bx;
            int y = my*pset.y(i) + by;
            int z = mz*pset.z(i) + bz;

            for (int y1 = std::max(y,0); y1 <= std::min(y+1,(int) dsmImage.height-1); ++y1) {
                for (int x1 = std::max(x,0); x1 <= std::min(x+1,(int) dsmImage.width-1); ++x1) {
                    int z0 = dsmImage.data[y1][x1];
                    if (abs(z-z0) < dz_short) {
                        sumImage.data[y1][x1] += pset.i(i);
                        ++cntImage.data[y1][x1];
                        max_i = std::max(max_i, sumImage.data[y1][x1] / cntImage.data[y1][x1]);
                    }
                }
            }
        }

        // Divide & copy
        double scale = std::numeric_limits<unsigned short>::max() / max_i;
        for (unsigned int y = 0; y < intImage.height; ++y) {
            for (unsigned int x = 0; x < intImage.width; ++x) {
                intImage.data[y][x] = (unsigned short) (scale * sumImage.data[y][x] / cntImage.data[y][x]);
            }
        }

        // Fill small voids in the INT.
        intImage.fillVoidsPyramid(true, 2);
    }
}

void Shr3dder::createMinAGL() {
    const OrthoImage<unsigned short> &dsmImage = getDSM();
    const OrthoImage<unsigned short> &dtmImage = getDTM();
    printf("Creating MinAGL\n");

    unsigned int dz_short = (unsigned int) (dz_meters / dsmImage.scale);
    unsigned int agl_short = (unsigned int) (agl_meters / dsmImage.scale);

    // Create the Minimum AGL Image; image gets the minimum values that are greater than the DTM+AGL
    if (!is_pset_empty())    {
        // Initialize size and allocate zeros
        minAglImage = OrthoImage<unsigned short>(&dsmImage);

        // Determine scalings
        double mx = 1.0/minAglImage.gsd;
        double bx = mx*(pset.xOff - minAglImage.easting) - 0.5;
        double my =-1.0/minAglImage.gsd;
        double by = minAglImage.height-1 + my*(pset.yOff - minAglImage.northing) + 0.5;
        double mz = 1/minAglImage.scale;
        double bz = mz*(pset.zOff - minAglImage.offset);

        // Loop over points
        for (unsigned long i = 0; i < pset.numPoints; i++) {
            int x = mx*pset.x(i) + bx;
            int y = my*pset.y(i) + by;
            unsigned short z = mz*pset.z(i) + bz;

            for (int y1 = std::max(y,0); y1 <= std::min(y+1,(int) minAglImage.height-1); ++y1) {
                for (int x1 = std::max(x,0); x1 <= std::min(x+1,(int) minAglImage.width-1); ++x1) {
                    unsigned short& z0 = minAglImage.data[y1][x1];
                    const unsigned short& z1 = dtmImage.data[y1][x1];
                    if (z1 && (z > z1+agl_short)  && (!z0 || z < z0))
                        z0 = z;
                }
            }
        }
    } else {
        // For when the input is not a pointset
        minAglImage = dsmImage;
        for (unsigned int y=0; y<minAglImage.height; ++y)
            for (unsigned int x=0; x<minAglImage.width; ++x)
                if (minAglImage.data[y][x] <= dtmImage.data[y][x]+agl_short)
                    minAglImage.data[y][x] = 0;
    }

    // Median filter, replacing only points differing by more than the AGL threshold.
    //minAglImage.quantileFilter(1, agl_short, 0.4, 0, false);

    // Filter quantile
    float quantile = 0.5f;

    // Apply filter to the image
    minAglImage.filter([&](unsigned short* val, const unsigned short& ref, std::vector<unsigned short> &ngbrs) {
        // Find quantile
        size_t ix = std::min((size_t) floor(quantile * ngbrs.size()),ngbrs.size()-1);
        std::partial_sort(ngbrs.begin(), ngbrs.begin() + (ix + 1), ngbrs.end());
        short qValue = static_cast<short>(ngbrs[ix]);
        // Only replace if it differs by more than dz from the median
        if (ref && abs(qValue - static_cast<short>(ref)) > static_cast<long>(dz_short))
            *val = qValue;
    },2);

    // Fill small voids in the image.
    //minAglImage.fillVoidsPyramid(true, 2);
}

void Shr3dder::labelClasses() {
    const OrthoImage<unsigned short> &dsmImage = getDSM();
    const OrthoImage<unsigned short> &dtmImage = getDTM();
    const OrthoImage<unsigned short> &dsm2Image = getDSM2();
    const OrthoImage<unsigned int> &labelImage = getLBL();
    printf("Creating CLS\n");
    // Produce a classification raster image with LAS standard point classes.

    // Convert horizontal and vertical uncertainty values to bin units.
    int dh_bins = MAX(1, (int) floor(dh_meters / dsmImage.gsd));
    unsigned int maxTreeHeightScaled = (max_tree_height_meters / dtmImage.scale);
    unsigned int dz_short = (unsigned int) (dz_meters / dsmImage.scale);
    unsigned int agl_short = (unsigned int) (agl_meters / dsmImage.scale);

    classImage = OrthoImage<unsigned char>(&dsmImage);
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
}

void Shr3dder::labelBuildings() {
    const OrthoImage<unsigned char>& classImage = getCLS();
    printf("Labeling buildings\n");
    bldgImage = OrthoImage<unsigned char>(classImage);
    for (unsigned int j = 0; j < classImage.height; j++) {
        for (unsigned int i = 0; i < classImage.width; i++) {
            if (classImage.data[j][i] != LAS_BUILDING) bldgImage.data[j][i] = 0;
        }
    }
}

void Shr3dder::createOutlines() {
    printf("Creating building outlines\n");
    const shr3d::OrthoImage<unsigned char>&bldgImage = getBLDG();
    
    bldgLabels = shr3d::OrthoImage<int>(&bldgImage);
    bldgLabels.labelConnectedComponentsFrom(&bldgImage);

    bldgLabels3.nn_upsample(&bldgLabels,3);

    std::map<int,GeoPolygon<double>> bounds = GeoPolygon<double>::traceBoundaries(bldgLabels3);
    printf("Traced %lu building outlines.\n",bounds.size());

    bldgOutlines.clear();
    for (std::pair<int,GeoPolygon<double>> pr : bounds) {
        GeoPolygon<double> poly = pr.second.buildingSimplify();
        if (poly.ring.empty())
            continue;
        bldgOutlines[pr.first] = poly;
    }
    printf("After simplifying, left with %lu building outlines.\n",bldgOutlines.size());
    
    bldgProcessed = true;
}


// Extend object boundaries to capture points missed around the edges.
void extendObjectBoundaries(const OrthoImage<unsigned short> &dsmImage, OrthoImage<unsigned int> &labelImage,
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
void labelObjectBoundaries(const OrthoImage<unsigned short> &dsmImage, OrthoImage<unsigned int> &labelImage,
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

bool findObjectBoundsInColumn(const OrthoImage<unsigned int> &labelImage, ObjectType &obj, unsigned int column, int &min_row, int &max_row) {
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
void fillObjectBounds(OrthoImage<unsigned int> &newLabelImage, const OrthoImage<unsigned int> &labelImage, const OrthoImage<unsigned short> &dsmImage, ObjectType &obj,
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
bool addNeighbors(std::vector<PixelType> &neighbors, OrthoImage<unsigned int> &labelImage,
                  const OrthoImage<unsigned short> &dsmImage, ObjectType &obj, unsigned int dzShort) {
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
void groupObjects(OrthoImage<unsigned int> &labelImage, const OrthoImage<unsigned short> &dsmImage,
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
void finishLabelImage(OrthoImage<unsigned int> &labelImage) {
    for (unsigned int j = 0; j < labelImage.height; j++) {
        for (unsigned int i = 0; i < labelImage.width; i++) {
            if ((labelImage.data[j][i] != LABEL_GROUND)) {
                labelImage.data[j][i] = 1;
            }
        }
    }
}

// Classify ground points, fill the voids, and generate a bare earth terrain model. 
void Shr3dder::classifyGround(OrthoImage<unsigned int> &labelImage, const OrthoImage<unsigned short> &dsmImage,
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
        labelObjectBoundaries(dtmImage, labelImage, dhBins, dzShort); // Updates labelImage

        // Extend labels for object boundaries..
        printf("Extending object boundaries...\n");
        extendObjectBoundaries(dtmImage, labelImage, dhBins, dzShort); // Updates labelImage

        // Group the objects.
        printf("Grouping objects...\n");
        std::vector<ObjectType> objects;
        groupObjects(labelImage, dtmImage, objects, maxCount, dzShort); // Updates labelImage, objects
        printf("Number of objects = %ld\n", objects.size());

        // Generate object groups and void fill them in the DEM image.
        printf("Labeling and removing objects...\n");
        OrthoImage<unsigned int> newlabelImage(labelImage);
        for (size_t i = 0; i < objects.size(); i++) {
            fillObjectBounds(newlabelImage,labelImage, dtmImage, objects[i], dhBins); // Updates newlabelImage, objects
        }
        std::swap(newlabelImage,labelImage);

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
void Shr3dder::classifyNonGround(const OrthoImage<unsigned short> &dsmImage, const OrthoImage<unsigned short> &dtmImage,
                                 OrthoImage<unsigned int> &labelImage, unsigned int dzShort, unsigned int aglShort,
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
