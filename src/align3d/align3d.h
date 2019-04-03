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