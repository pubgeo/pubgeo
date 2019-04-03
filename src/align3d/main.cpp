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

#include "align3d.h"
#include <chrono>

void printArguments();

// Main program for simple 3d alignment.
int main(int argc, char **argv) {
    // If no parameters, then print command line arguments.
    if (argc < 3) {
        printf("\nNumber of arguments = %d\n", argc);
        for (int i = 1; i < argc; i++) {
            printf("  ARG[%d] = %s\n", i, argv[i]);
        }
        printArguments();
        return -1;
    }

    // Parse command line parameters.
    char referenceFileName[1024];
    char targetFileName[1024];
    strcpy(referenceFileName, argv[1]);
    strcpy(targetFileName, argv[2]);
    align3d::AlignParameters params;
    params.gsd = 1.0;
    params.maxt = 10.0;
    params.maxdz = 0.0;
    for (int i = 2; i < argc; i++) {
        if (strstr(argv[i], "maxdz=")) { params.maxdz = (float) atof(&(argv[i][6])); }
        if (strstr(argv[i], "gsd=")) { params.gsd = (float) atof(&(argv[i][4])); }
        if (strstr(argv[i], "maxt=")) { params.maxt = (float) atof(&(argv[i][5])); }
    }

    // Default MAXDZ = GSD x 2 to ensure reliable performance on steep slopes.
    if (params.maxdz == 0.0) params.maxdz = params.gsd * 2.0;

    printf("Selected Parameters:\n");
    printf("  ref   = %s\n", referenceFileName);
    printf("  tgt   = %s\n", targetFileName);
    printf("  gsd   = %f\n", params.gsd);
    printf("  maxdz = %f\n", params.maxdz);
    printf("  maxt  = %f\n", params.maxt);

    // Initialize the timer.
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();

    try {
        // Align the target point cloud to the reference.
        AlignTarget2Reference(referenceFileName, targetFileName, params);
    } catch (char const *err) {
        std::cerr << "Reporting error: " << err << std::endl;
    }
    // Report total elapsed time.
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    printf("Total time elapsed = %f seconds\n", std::chrono::duration<double>(t1-t0).count());

    return 0;
}

// Print command line arguments.
void printArguments() {
    printf("Command Line: align3d <reference> <target> <parameters>\n");
    printf("Parameters:\n");
    printf("  gsd=   Ground Sample Distance (GSD) for gridding point cloud (meters); default = 1.0\n");
    printf("  maxdz= Max local Z difference (meters) for matching; default = 2*gsd\n");
    printf("  maxt=  Maximum horizontal translation in search (meters); default = 10.0\n");
    printf("Examples:\n");
    printf("  align3d ref.las tgt.las maxt=10.0 gsd=0.5 maxdz=0.5\n\n");
}

