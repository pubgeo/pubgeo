// Copyright 2017 The Johns Hopkins University Applied Physics Laboratory.
// Licensed under the MIT License. See LICENSE.txt in the project root for full license information.

#include <PointCloud.h>

using namespace pubgeo;

int main() {
    if (!PointCloud::TransformPointCloud("someFile.las", "someOut.las", 3.0, 4.0, 5.0)) {
        return 1;
    }

    return 0;
}