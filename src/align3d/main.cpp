// Copyright 2017 The Johns Hopkins University Applied Physics Laboratory.
// Licensed under the MIT License. See LICENSE.txt in the project root for full license information.

#include "align3d.h"

void printArguments();

// Main program for simple 3d alignment.
int main(int argc, char **argv) {
    // If no parameters, then print command line arguments.
    if (argc < 4) {
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
    time_t t0;
    time(&t0);

    try {
        // Align the target point cloud to the reference.
        AlignTarget2Reference(referenceFileName, targetFileName, params);
    } catch (char const *err) {
        std::cerr << "Reporting error: " << err << std::endl;
    }
    // Report total elapsed time.
    time_t t1;
    time(&t1);
    printf("Total time elapsed = %f seconds\n", double(t1 - t0));

    return 0;
}

// Print command line arguments.
void printArguments() {
    printf("Command Line: align-3d <reference> <target> <parameters>\n");
    printf("Parameters:\n");
    printf("  maxdz= Max local Z difference (meters) for matching\n");
    printf("  gsd=   Ground Sample Distance (GSD) for gridding (meters)\n");
    printf("  maxt=	 Maximum XYZ translation in search (meters); default = 10.0\n");
    printf("Examples:\n");
    printf("  align-3d ref.las tgt.las maxt=10.0 gsd=0.5 maxdz=0.5 \n\n");
}
