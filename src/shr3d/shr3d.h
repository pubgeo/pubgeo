// Copyright 2017 The Johns Hopkins University Applied Physics Laboratory.
// Licensed under the MIT License. See LICENSE.txt in the project root for full license information.

// shr3d.h
//

#ifndef PUBGEO_SHR3D_H
#define PUBGEO_SHR3D_H

#include <cstdio>
#include <cstring>
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

    class Shr3dder {
    public:
        // Function declarations.
        static void classifyGround(OrthoImage<unsigned long> &labelImage, OrthoImage<unsigned short> &dsmImage,
                                   OrthoImage<unsigned short> &dtmImage, int dhBins, unsigned int dzShort);

        static void classifyNonGround(OrthoImage<unsigned short> &dsmImage, OrthoImage<unsigned short> &dtmImage,
                                      OrthoImage<unsigned long> &labelImage, unsigned int dzShort,
                                      unsigned int aglShort,
                                      float minAreaMeters);

        static void fillInsideBuildings(OrthoImage<unsigned char> &classImage);
    };
}

#endif // PUBGEO_SHR3D_H