// Copyright 2017 The Johns Hopkins University Applied Physics Laboratory.
// Licensed under the MIT License. See LICENSE.txt in the project root for full license information.

// shr3d.h
//

#ifndef PUBGEO_SHR3D_H
#define PUBGEO_SHR3D_H

#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include "orthoimage.h"

#define LABEL_OBJECT    1
#define LABEL_GROUND    UINT_MAX
#define LABEL_IN_ONE    UINT_MAX-1
#define LABEL_ACCEPTED  UINT_MAX-2
#define LABEL_TEMP      UINT_MAX-3
#define LABEL_BOUNDARY  UINT_MAX-5

/*
 * Las Classifications list
    0 Created, never classified
    1 Unclassified1
    2 Ground
    3 Low Vegetation
    4 Medium Vegetation
    5 High Vegetation
    6 Building
    7 Low Point (noise)
    8 Model Key-point (mass point)
    9 Water
    10 Reserved for ASPRS Definition
    11 Reserved for ASPRS Definition
    12 Overlap Points2
    13-31 Reserved for ASPRS Definition
 */
#define LAS_UNCLASSIFIED    (1*40)
#define LAS_GROUND          (2*40)
#define LAS_TREE            (5*40)
#define LAS_BUILDING        (6*40)
// For now, multiply these by 40 so we can view the images easily.

namespace shr3d {
    using namespace pubgeo;

    typedef struct {
        unsigned int i;
        unsigned int j;
    } PixelType;

    typedef struct {
        unsigned int label;
        long int xmin;
        long int xmax;
        long int ymin;
        long int ymax;
        long int count;
    } ObjectType;

    enum ImageType {
        DSM,
        MIN,
        DTM,
        CLASS,
        BUILDING,
        BUILDING_OUTLINES,
        DSM2,
        LABEL,
        LABELED_BUILDINGS,
        LABELED_BUILDINGS_3,
    };

    class Shr3dder {
    public:
        // Options
        double dh_meters;
        double dz_meters;
        double agl_meters;
        double min_area_meters;
        double max_tree_height_meters;
        bool egm96;

        // Constructor
        Shr3dder() : dh_meters(1), dz_meters(1), agl_meters(2),
                min_area_meters(50), max_tree_height_meters(40), egm96(false) {}

        // Function declarations.
        void process(const OrthoImage<unsigned short> &dsmImage, const OrthoImage<unsigned short> &minImage,
                std::map<ImageType,std::string> outputFilenames);

        bool createDSM(const PointCloud& pset, OrthoImage<unsigned short> &dsmImage);

        bool createMIN(const PointCloud& pset, OrthoImage<unsigned short> &minImage);

        void createDTM(const OrthoImage<unsigned short> &dsmImage,  const OrthoImage<unsigned short> &minImage,
                OrthoImage<unsigned short> &dtmImage, OrthoImage<unsigned short> &dsm2Image, OrthoImage<unsigned long> &labelImage);

        OrthoImage<unsigned char> labelClasses(const OrthoImage<unsigned short> &dsmImage, const OrthoImage<unsigned short> &dtmImage,
                       const OrthoImage<unsigned short> &dsm2Image, const OrthoImage<unsigned long> &labelImage);

        OrthoImage<unsigned char> labelBuildings(const OrthoImage<unsigned char> &classImage);

        void classifyGround(OrthoImage<unsigned long> &labelImage, OrthoImage<unsigned short> &dsmImage,
                                   OrthoImage<unsigned short> &dtmImage, int dhBins, unsigned int dzShort);

        void classifyNonGround(OrthoImage<unsigned short> &dsmImage, OrthoImage<unsigned short> &dtmImage,
                                      OrthoImage<unsigned long> &labelImage, unsigned int dzShort,
                                      unsigned int aglShort,
                                      float minAreaMeters);

        void fillInsideBuildings(OrthoImage<unsigned char> &classImage);
    };
}

#endif // PUBGEO_SHR3D_H
