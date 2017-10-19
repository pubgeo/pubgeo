// Copyright 2017 The Johns Hopkins University Applied Physics Laboratory.
// Licensed under the MIT License. See LICENSE.txt in the project root for full license information.

// Created by almess1 on 4/21/17.
//

#ifndef PUBGEO_IMAGE_H
#define PUBGEO_IMAGE_H

#include <cstring>
namespace pubgeo {
    template<class TYPE>
    class Image {

    public:

        unsigned int width;
        unsigned int height;
        unsigned int bands;
        float scale;
        float offset;
        TYPE **data;    // Order is HEIGHT (outside), WIDTH, BAND (inside).

        // Deallocate memory.
        void Deallocate() {
            if (!data) return;
            for (unsigned int y = 0; y < height; y++) delete[]data[y];
            delete[]data;
            data = nullptr;
        }

        // Allocate memory.
        void Allocate(unsigned int numColumns, unsigned int numRows, unsigned int numBands = 1) {
            if (data) Deallocate();
            width = numColumns;
            height = numRows;
            bands = numBands;
            data = new TYPE *[height];
            for (unsigned int y = 0; y < height; y++) {
                data[y] = new TYPE[width * bands];
                std::memset(data[y], 0, width * bands * sizeof(TYPE));
            }
        }

        // Perform a deep copy of an image
        void Clone(const Image<TYPE>* im) {
            Deallocate();

            // Copy members
            width = im->width;
            height = im->height;
            bands = im->bands;
            scale = im->scale;
            offset = im->offset;

            // Deep copy of data
            // Note: not using Allocate so we can skip the memset
            data = new TYPE *[height];
            for (unsigned int y = 0; y < height; y++) {
                data[y] = new TYPE[width * bands];
                std::memcpy(data[y], im->data[y], width * bands * sizeof(TYPE));
            }
        }
    };
}
#endif //PUBGEO_IMAGE_H
