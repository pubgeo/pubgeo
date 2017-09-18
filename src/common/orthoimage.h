// Copyright 2017 The Johns Hopkins University Applied Physics Laboratory.
// Licensed under the MIT License. See LICENSE.txt in the project root for full license information.

// orthoimage.h
//


#ifndef PUBGEO_ORTHO_IMAGE_H
#define PUBGEO_ORTHO_IMAGE_H

#include <stdio.h>
#include <math.h>
#include <vector>
#include <typeinfo>
#include <cstring>
#include <algorithm>

#ifdef WIN32
#include "gdal_priv.h"
#include "ogr_spatialref.h"
#else

#include <gdal/gdal_priv.h>
#include "gdal/ogr_spatialref.h"

#endif

#include "util.h"
#include "PointCloud.h"
#include "Image.h"

namespace pubgeo {
    typedef enum {
        MIN_VALUE, MAX_VALUE
    } MIN_MAX_TYPE;

//
// Ortho image template class
//
    template<class TYPE>
    class OrthoImage : public Image<TYPE> {

    public:

        double easting;
        double northing;
        int zone;
        float gsd;

        // Constructor
        OrthoImage() {
            memset(this, 0, sizeof(OrthoImage));
            this->scale = 1.0f;
        }

        // Destructor
        ~OrthoImage() {
            this->Deallocate();
        }

        // Read any GDAL-supported image.
        bool read(char *fileName) {
            // Open the image.
            GDALAllRegister();
            CPLSetConfigOption("GDAL_DATA", ".\\gdaldata");
            GDALDataset *poDataset = (GDALDataset *) GDALOpen(fileName, GA_ReadOnly);
            if (poDataset == NULL) {
                printf("Error opening %s.\n", fileName);
                return false;
            }

            // Get geospatial metadata.
            double adfGeoTransform[6];
            printf("Driver: %s/%s\n", poDataset->GetDriver()->GetDescription(),
                   poDataset->GetDriver()->GetMetadataItem(GDAL_DMD_LONGNAME));
            this->width = (unsigned int) poDataset->GetRasterXSize();
            this->height = (unsigned int) poDataset->GetRasterYSize();
            this->bands = (unsigned int) poDataset->GetRasterCount();
            printf("Width = %d\nHeight = %d\nBands = %d\n", this->width, this->height, this->bands);
            if (poDataset->GetProjectionRef() != NULL) printf("Projection is %s.\n", poDataset->GetProjectionRef());
            if (poDataset->GetGeoTransform(adfGeoTransform) == CE_None) {
                printf("GeoTransform = (%.6f,%.6f,%.6f,%.6f,%.6f,%.6f)\n",
                       adfGeoTransform[0], adfGeoTransform[1], adfGeoTransform[2], adfGeoTransform[3],
                       adfGeoTransform[4],
                       adfGeoTransform[5]);
            }
            double xscale = adfGeoTransform[1];
            double yscale = -adfGeoTransform[5];
            this->gsd = (float) ((xscale + yscale) / 2.0);
            this->easting = adfGeoTransform[0] - adfGeoTransform[2] * xscale;
            this->northing = adfGeoTransform[3] - (this->height) * yscale;
            OGRSpatialReference myOGRS(poDataset->GetProjectionRef());
            int north = 1;
            this->zone = myOGRS.GetUTMZone(&north);
            if (north == 0) {

                this->zone *= -1;
            }
            printf("UTM Easting = %lf\n", this->easting);
            printf("UTM Northing = %lf\n", this->northing);
            printf("UTM Zone = %d\n", this->zone);
            printf("GSD = %f\n", this->gsd);

            // Get band information.
            GDALRasterBand *poBand = poDataset->GetRasterBand(1);
            if (poBand == NULL) {
                printf("Error opening first band.");
                return false;
            }
            GDALDataType BandDataType = poBand->GetRasterDataType();
            int bytesPerPixel = GDALGetDataTypeSize(BandDataType) / 8;
            unsigned char *poBandBlock = (unsigned char *) CPLMalloc(bytesPerPixel * this->width);
            int ok = 0;
            double noData = poBand->GetNoDataValue(&ok);
            if (!ok) {
                // Set noData only for floating point images.
                if (strcmp(typeid(TYPE).name(), "float") == 0)
                    noData = -10000.0;
                else
                    noData = 0;
            }

            // Get scale and offset values.
            if (strcmp(typeid(TYPE).name(), "float") == 0) {
                // Do not scale if floating point values.
                this->scale = 1.0;
                this->offset = 0.0;
            } else {
                float minVal = MAX_FLOAT;
                float maxVal = MAX_FLOAT;
                for (unsigned int i = 0; i < this->bands; i++) {
                    poBand = poDataset->GetRasterBand(i + 1);
                    if (poBand == NULL) {
                        printf("Error opening band %d", i + 1);
                        return false;
                    }
                    int bGotMin, bGotMax;
                    double adfMinMax[2];
                    adfMinMax[0] = poBand->GetMinimum(&bGotMin);
                    adfMinMax[1] = poBand->GetMaximum(&bGotMax);
                    if (!(bGotMin && bGotMax)) GDALComputeRasterMinMax((GDALRasterBandH) poBand, false, adfMinMax);
                    if (minVal > adfMinMax[0]) minVal = (float) adfMinMax[0];
                    if (maxVal < adfMinMax[1]) maxVal = (float) adfMinMax[1];
                }
                minVal -= 1;    // Reserve zero for noData value
                maxVal += 1;
                float maxImageVal = (float) (pow(2.0, int(sizeof(TYPE) * 8)) - 1);
                this->offset = minVal;
                this->scale = (maxVal - minVal) / maxImageVal;
            }
            printf("Offset = %f\n", this->offset);
            printf("Scale = %f\n", this->scale);

            // Allocate memory for the image.
            this->Allocate(this->width, this->height, this->bands);

            // Read the image, one band at a time.
            for (unsigned int irow = 0; irow < this->height; irow++) {
                for (unsigned int ib = 0; ib < this->bands; ib++) {
                    // Read the next row.
                    poBand = poDataset->GetRasterBand(ib + 1);
                    poBand->RasterIO(GF_Read, 0, irow, this->width, 1, poBandBlock, this->width, 1, BandDataType, 0, 0);
                    switch (BandDataType) {
                        case GDT_Byte: {
                            for (unsigned int icols = 0; icols < this->width; icols++) {
                                long k = icols * this->bands;
                                if (poBandBlock[icols] == (unsigned char) noData)
                                    this->data[irow][k + ib] = 0;
                                else
                                    this->data[irow][k + ib] = (TYPE) ((poBandBlock[icols] - this->offset) /
                                                                       this->scale);
                            }
                            break;
                        }
                        case GDT_UInt16: {
                            unsigned int *poBandBlockShort = (unsigned int *) poBandBlock;
                            for (unsigned int icols = 0; icols < this->width; icols++) {
                                long k = icols * this->bands;
                                if (poBandBlockShort[icols] == (unsigned int) noData)
                                    this->data[irow][k + ib] = 0;
                                else
                                    this->data[irow][k + ib] = (TYPE) ((poBandBlockShort[icols] - this->offset) /
                                                                       this->scale);
                            }
                            break;
                        }
                        case GDT_Float32: {
                            float *poBandBlockFloat = (float *) poBandBlock;
                            for (unsigned int icols = 0; icols < this->width; icols++) {
                                long k = icols * this->bands;
                                if (poBandBlockFloat[icols] == (float) noData)
                                    this->data[irow][k + ib] = 0;
                                else
                                    this->data[irow][k + ib] = (TYPE) ((poBandBlockFloat[icols] - this->offset) /
                                                                       this->scale);
                            }
                            break;
                        }
                        default: {
                            // TODO: shouldn't get here, relay error
                        }
                    }
                }
            }

            // Deallocate memory.
            CPLFree(poBandBlock);
            delete poDataset;
            return true;
        }

        // Write GEOTIFF image using GDAL.
        bool write(char *fileName, bool convertToFloat = false, bool egm96 = false) {
            GDALAllRegister();
            const char *pszFormat = "GTiff";
            GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
            if (poDriver == NULL) return false;
            char **papszMetadata = poDriver->GetMetadata();
            if (!CSLFetchBoolean(papszMetadata, GDAL_DCAP_CREATE, false)) return false;

            // Get the GDAL data type to match the image data type.
            GDALDataType theBandDataType = (GDALDataType) sizeof(TYPE);
            if (strcmp(typeid(TYPE).name(), "float") == 0) theBandDataType = GDT_Float32;

            // If converting this image to FLOAT, then update the output format type.
            if (convertToFloat) theBandDataType = GDT_Float32;

            // Write geospatial metadata.
            char **papszOptions = nullptr;
            GDALDataset *poDstDS = poDriver->Create(fileName, this->width, this->height, this->bands, theBandDataType,
                                                    papszOptions);
            double adfGeoTransform[6] = {this->easting, this->gsd, 0, this->northing + this->height * this->gsd, 0,
                                         -1 * this->gsd};
            poDstDS->SetGeoTransform(adfGeoTransform);
            OGRSpatialReference oSRS;
            char *pszSRS_WKT = NULL;
            oSRS.SetProjCS("UTM (WGS84)");
            oSRS.SetUTM(abs(this->zone), (this->zone > 0));
            int theZone = 255;
            int theZoneIsNorthern = 255;
            theZone = oSRS.GetUTMZone(&theZoneIsNorthern);
            oSRS.SetWellKnownGeogCS("WGS84");
            if (!egm96)
                oSRS.SetVertCS("WGS 84", "World Geodetic System 1984");
            else
                oSRS.SetVertCS("EGM96 Geoid", "EGM96 geoid");
            oSRS.SetExtension("VERT_DATUM", "PROJ4_GRIDS", "g2009conus.gtx");
            oSRS.exportToWkt(&pszSRS_WKT);
            poDstDS->SetProjection(pszSRS_WKT);
            CPLFree(pszSRS_WKT);

            // Write the image.
            if (convertToFloat) {
                // Write the image as FLOAT and convert values.
                float *raster = new float[this->width * this->bands];
                float noData = -10000.0;
                for (unsigned int i = 1; i <= this->bands; i++) {
                    GDALRasterBand *poBand = poDstDS->GetRasterBand(i);
                    poBand->SetNoDataValue(noData);
                    for (unsigned int j = 0; j < this->height; j++) {
                        for (unsigned int k = 0; k < this->width * this->bands; k++) {
                            if (this->data[j][k] == 0)
                                raster[k] = noData;
                            else
                                raster[k] = (float(this->data[j][k]) * this->scale) + this->offset;
                        }
                        poBand->RasterIO(GF_Write, 0, j, this->width, 1, raster, this->width, 1, theBandDataType,
                                         sizeof(float) * this->bands,
                                         this->width * this->bands * sizeof(float));
                    }
                }
                delete[]raster;
            } else {
                // Write the image without conversion.
                for (unsigned int i = 1; i <= this->bands; i++) {
                    GDALRasterBand *poBand = poDstDS->GetRasterBand(i);
                    for (unsigned int j = 0; j < this->height; j++) {
                        poBand->RasterIO(GF_Write, 0, j, this->width, 1, this->data[j], this->width, 1, theBandDataType,
                                         sizeof(TYPE) * this->bands,
                                         this->width * this->bands * sizeof(TYPE));
                    }
                }
            }
            GDALClose((GDALDatasetH) poDstDS);
            return true;
        }

        // Read image from point cloud file.
        bool readFromPointCloud(char *fileName, float gsdMeters, MIN_MAX_TYPE mode = MIN_VALUE) {
            // Read a PSET file (e.g., BPF or LAS).
            PointCloud pset;
            bool ok = pset.Read(fileName);
            if (!ok) return false;

            // Calculate scale and offset for conversion to TYPE.
            float minVal = pset.bounds.zMin - 1;    // Reserve zero for noData value
            float maxVal = pset.bounds.zMax + 1;
            float maxImageVal = (float) (pow(2.0, int(sizeof(TYPE) * 8)) - 1);
            this->offset = minVal;
            this->scale = (maxVal - minVal) / maxImageVal;

            // Calculate image width and height.
            this->width = (unsigned int) ((pset.bounds.xMax - pset.bounds.xMin) / gsdMeters + 1);
            this->height = (unsigned int) ((pset.bounds.yMax - pset.bounds.yMin) / gsdMeters + 1);

            // Allocate an ortho image.
            this->Allocate(this->width, this->height);
            this->easting = pset.bounds.xMin;
            this->northing = pset.bounds.yMin;
            this->zone = pset.zone;
            this->gsd = gsdMeters;

            // Copy points into the ortho image.
            if (mode == MIN_VALUE) {
                for (unsigned long i = 0; i < pset.numPoints; i++) {
                    unsigned int x = int((pset.x(i) + pset.xOff - easting) / gsd + 0.5);
                    if ((x < 0) || (x > this->width - 1)) continue;
                    unsigned int y = this->height - 1 - int((pset.y(i) + pset.yOff - northing) / gsd + 0.5);
                    if ((y < 0) || (y > this->height - 1)) continue;
                    TYPE z = TYPE((pset.z(i) + pset.zOff - this->offset) / this->scale);
                    if ((this->data[y][x] == 0) || (z < this->data[y][x])) this->data[y][x] = z;
                }
            } else if (mode == MAX_VALUE) {
                for (unsigned long i = 0; i < pset.numPoints; i++) {
                    unsigned int x = int((pset.x(i) + pset.xOff - easting) / gsd + 0.5);
                    if ((x < 0) || (x > this->width - 1)) continue;
                    unsigned int y = this->height - 1 - int((pset.y(i) + pset.yOff - northing) / gsd + 0.5);
                    if ((y < 0) || (y > this->height - 1)) continue;
                    TYPE z = TYPE((pset.z(i) + pset.zOff - this->offset) / this->scale);
                    if ((this->data[y][x] == 0) || (z > this->data[y][x])) this->data[y][x] = z;
                }
            }
            return true;
        }

        // Read image from PDAL PointView.
        bool readFromPointView(pdal::PointViewPtr view, float gsdMeters, MIN_MAX_TYPE mode = MIN_VALUE) {
            // Read a PSET file (e.g., BPF or LAS).
            PointCloud pset;
            bool ok = pset.Read(view);
            if (!ok) return false;

            // Calculate scale and offset for conversion to TYPE.
            float minVal = pset.bounds.zMin - 1;    // Reserve zero for noData value
            float maxVal = pset.bounds.zMax + 1;
            float maxImageVal = (float) (pow(2.0, int(sizeof(TYPE) * 8)) - 1);
            this->offset = minVal;
            this->scale = (maxVal - minVal) / maxImageVal;

            // Calculate image width and height.
            this->width = (unsigned int) ((pset.bounds.xMax - pset.bounds.xMin) / gsdMeters + 1);
            this->height = (unsigned int) ((pset.bounds.yMax - pset.bounds.yMin) / gsdMeters + 1);

            // Allocate an ortho image.
            this->Allocate(this->width, this->height);
            this->easting = pset.bounds.xMin;
            this->northing = pset.bounds.yMin;
            this->zone = pset.zone;
            this->gsd = gsdMeters;

            // Copy points into the ortho image.
            if (mode == MIN_VALUE) {
               for (unsigned long i = 0; i < pset.numPoints; i++) {
                    double dx = view->getFieldAs<double>(pdal::Dimension::Id::X, i);
                    double dy = view->getFieldAs<double>(pdal::Dimension::Id::Y, i);
                    double dz = view->getFieldAs<double>(pdal::Dimension::Id::Z, i);
                    unsigned int x = int((dx - easting) / gsd + 0.5);
                    if ((x < 0) || (x > this->width - 1)) continue;
                    unsigned int y = this->height - 1 - int((dy - northing) / gsd + 0.5);
                    if ((y < 0) || (y > this->height - 1)) continue;
                    TYPE z = TYPE((dz - this->offset) / this->scale);
                    if ((this->data[y][x] == 0) || (z < this->data[y][x])) this->data[y][x] = z;
                }
            } else if (mode == MAX_VALUE) {
                for (unsigned long i = 0; i < pset.numPoints; i++) {
                    double dx = view->getFieldAs<double>(pdal::Dimension::Id::X, i);
                    double dy = view->getFieldAs<double>(pdal::Dimension::Id::Y, i);
                    double dz = view->getFieldAs<double>(pdal::Dimension::Id::Z, i);
                    unsigned int x = int((dx - easting) / gsd + 0.5);
                    if ((x < 0) || (x > this->width - 1)) continue;
                    unsigned int y = this->height - 1 - int((dy - northing) / gsd + 0.5);
                    if ((y < 0) || (y > this->height - 1)) continue;
                    TYPE z = TYPE((dz - this->offset) / this->scale);
                    if ((this->data[y][x] == 0) || (z > this->data[y][x])) this->data[y][x] = z;
                }
            }
            return true;
        }

        // Count voids in an image.
        // Note that voids are always labeled zero in this data structure.
        long countVoids() {
            long count = 0;
            for (unsigned int j = 0; j < this->height; j++) {
                for (unsigned int i = 0; i < this->width; i++) {
                    if (this->data[j][i] == 0) {
                        count++;
                    }
                }
            }
            return count;
        }

        // Fill any voids in the image using a simple multigrid scheme.
        // Note that voids are always labeled zero.
        // MaxLevel by default is the maximum value of int
        void fillVoidsPyramid(bool noSmoothing, unsigned int maxLevel = MAX_INT) {
            // Check for voids.
            long count = countVoids();
            if (count == 0) return;

            // Create image pyramid.
            std::vector<OrthoImage<TYPE> *> pyramid;
            pyramid.push_back(this);
            unsigned int level = 0;
            while ((count > 0) && (level < maxLevel)) {
                // Create next level.
                OrthoImage<TYPE> *newImagePtr = new OrthoImage<TYPE>;

                unsigned int nextWidth = pyramid[level]->width / 2;
                unsigned int nextHeight = pyramid[level]->height / 2;
                newImagePtr->Allocate(nextWidth, nextHeight, 1);

                // Fill in non-void values from level below building up the pyramid with a simple running average.
                for (unsigned int j = 0; j < nextHeight; j++) {
                    for (unsigned int i = 0; i < nextWidth; i++) {
                        unsigned int j2 = MIN(MAX(0, j * 2 + 1), pyramid[level]->height - 1);
                        unsigned int i2 = MIN(MAX(0, i * 2 + 1), pyramid[level]->width - 1);

                        // Average neighboring pixels from below.
                        float z = 0;
                        int ct = 0;
                        std::vector<TYPE> neighbors;
                        for (unsigned int jj = MAX(0, j2 - 1); jj <= MIN(j2 + 1, pyramid[level]->height - 1); jj++) {
                            for (unsigned int ii = MAX(0, i2 - 1); ii <= MIN(i2 + 1, pyramid[level]->width - 1); ii++) {
                                if (pyramid[level]->data[jj][ii] != 0) {
                                    z += pyramid[level]->data[jj][ii];
                                    ct++;
                                }
                            }
                        }
                        if (ct != 0) {
                            z = z / ct;
                            newImagePtr->data[j][i] = (TYPE) z;
                        }
                    }
                }

                pyramid.push_back(newImagePtr);
                level++;
                count = pyramid[level]->countVoids();
            }

            // Void fill down the pyramid.
            for (int k = level - 1; k >= 0; k--) {
                for (unsigned int j = 0; j < pyramid[k]->height; j++) {
                    for (unsigned int i = 0; i < pyramid[k]->width; i++) {
                        // Fill this pixel if it is currently void.
                        if (pyramid[k]->data[j][i] == 0) {
                            unsigned int j2 = MIN(MAX(0, j / 2), pyramid[k + 1]->height - 1);
                            unsigned int i2 = MIN(MAX(0, i / 2), pyramid[k + 1]->width - 1);

                            if (noSmoothing) {
                                // Just use the closest pixel from above.
                                pyramid[k]->data[j][i] = pyramid[k + 1]->data[j2][i2];
                            } else {
                                // Average neighboring pixels from above.
                                float z = 0;
                                int ct = 0;
                                for (unsigned int jj = MAX(0, j2 - 1);
                                     jj <= MIN(j2 + 1, pyramid[k + 1]->height - 1); jj++) {
                                    for (unsigned int ii = MAX(0, i2 - 1);
                                         ii <= MIN(i2 + 1, pyramid[k + 1]->width - 1); ii++) {
                                        z += pyramid[k + 1]->data[jj][ii];
                                        ct++;
                                    }
                                }
                                z = z / ct;
                                pyramid[k]->data[j][i] = (TYPE) z;
                            }
                        }
                    }
                }
            }

            // Deallocate memory for all but the input DSM.
            for (unsigned int i = 1; i <= level; i++) {
                pyramid[i]->Deallocate();
            }
        }

        // Apply a median filter to an image.
        void medianFilter(int rad, unsigned int dzShort) {
            for (unsigned int j = 0; j < this->height; j++) {
                for (unsigned int i = 0; i < this->width; i++) {
                    // Skip if void.
                    if (this->data[j][i] == 0) continue;

                    // Define bounds;
                    unsigned int i1 = MAX(0, i - rad);
                    unsigned int i2 = MIN(i + rad, this->width - 1);
                    unsigned int j1 = MAX(0, j - rad);
                    unsigned int j2 = MIN(j + rad, this->height - 1);

                    // Add valid values to the list.
                    std::vector<unsigned short> values;
                    for (unsigned int jj = j1; jj <= j2; jj++) {
                        for (unsigned int ii = i1; ii <= i2; ii++) {
                            if (this->data[jj][ii] != 0) {
                                values.push_back(this->data[jj][ii]);
                            }
                        }
                    }

                    // Find the median value.
                    unsigned long num = values.size();
                    if (num > 0) {
                        std::sort(values.begin(), values.end());
                        unsigned short medianValue = values[num / 2];
                        if (fabs(float(medianValue) - float(this->data[j][i])) > dzShort)
                            this->data[j][i] = medianValue;
                    }
                }
            }
        }

        // Apply a minimum filter to an image.
        void minFilter(int rad) {
            OrthoImage<TYPE> tempImage;
            tempImage.Allocate(this->width, this->height);
            for (unsigned int j = 0; j < this->height; j++) {
                for (unsigned int i = 0; i < this->width; i++) {
                    tempImage.data[j][i] = this->data[j][i];
                }
            }
            for (unsigned int j = 0; j < this->height; j++) {
                for (unsigned int i = 0; i < this->width; i++) {
                    // Skip if void.
                    if (this->data[j][i] == 0) continue;

                    // Define bounds;
                    unsigned int i1 = MAX(0, i - rad);
                    unsigned int i2 = MIN(i + rad, this->width - 1);
                    unsigned int j1 = MAX(0, j - rad);
                    unsigned int j2 = MIN(j + rad, this->height - 1);

                    // Check all neighbors.
                    for (unsigned int jj = j1; jj <= j2; jj++) {
                        for (unsigned int ii = i1; ii <= i2; ii++) {
                            if (this->data[jj][ii] != 0) {
                                if (this->data[jj][ii] < tempImage.data[j][i])
                                    tempImage.data[j][i] = this->data[jj][ii];
                            }
                        }
                    }
                }
            }
            for (unsigned int j = 0; j < this->height; j++) {
                for (unsigned int i = 0; i < this->width; i++) {
                    this->data[j][i] = tempImage.data[j][i];
                }
            }
        }

        void edgeFilter(int dzShort) {
            OrthoImage<TYPE> tempImage;
            tempImage.Allocate(this->width, this->height);
            for (unsigned int j = 0; j < this->height; j++) {
                for (unsigned int i = 0; i < this->width; i++) {
                    tempImage.data[j][i] = this->data[j][i];
                }
            }
            for (unsigned int j = 0; j < this->height; j++) {
                for (unsigned int i = 0; i < this->width; i++) {
                    // Skip if void.
                    if (this->data[j][i] == 0) continue;

                    // Define bounds;
                    unsigned int i1 = MAX(0, i - 1);
                    unsigned int i2 = MIN(i + 1, this->width - 1);
                    unsigned int j1 = MAX(0, j - 1);
                    unsigned int j2 = MIN(j + 1, this->height - 1);

                    // If there's an edge with any neighbor, then remove.
                    for (unsigned int jj = j1; jj <= j2; jj++) {
                        for (unsigned int ii = i1; ii <= i2; ii++) {
                            if (fabs(float(this->data[jj][ii] - this->data[j][i])) > dzShort) {
                                tempImage.data[j][i] = 0;
                            }
                        }
                    }
                }
            }
            for (unsigned int j = 0; j < this->height; j++) {
                for (unsigned int i = 0; i < this->width; i++) {
                    this->data[j][i] = tempImage.data[j][i];
                }
            }
        }
    };
}

#endif //PUBGEO_ORTHO_IMAGE_H
