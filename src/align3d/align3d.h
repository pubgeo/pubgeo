// Copyright 2017 The Johns Hopkins University Applied Physics Laboratory.
// Licensed under the MIT License. See LICENSE.txt in the project root for full license information.

#ifndef PUBGEO_ALIGN3D_H
#define PUBGEO_ALIGN3D_H

#include <vector>
#include "orthoimage.h"
#include "Image.h"

namespace align3d {
    using namespace pubgeo;

    typedef struct {
        float tx;
        float ty;
        float tz;
        float rms;
    }
            AlignResult;

    typedef struct {
        float gsd;
        float maxdz;
        float maxt;
    }
            AlignParameters;

    typedef struct {
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        double width;
        double height;
    }
            AlignBounds;

    bool computeRMS(float dx, float dy, long numSamples, long maxSamples, std::vector<double> &xlist,
                    std::vector<double> &ylist, OrthoImage<unsigned short> &referenceDSM,
                    OrthoImage<unsigned short> &targetDSM, float &medianDZ, float &rms, long &ndx, float &completeness);

    // Estimate 3D rigid body transform parameters to align target points with reference.
    void EstimateRigidBody(OrthoImage<unsigned short> &referenceDSM, OrthoImage<unsigned short> &targetDSM, float maxt,
                           AlignBounds &bounds, AlignResult &result);

    // Align target file to match reference file.
    bool AlignTarget2Reference(char *referenceFileName, char *targetFileName, AlignParameters params);
}

#endif // PUBGEO_ALIGN3D_H