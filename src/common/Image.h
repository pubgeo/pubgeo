// Copyright 2017 The Johns Hopkins University Applied Physics Laboratory.
// Licensed under the MIT License. See LICENSE.txt in the project root for full license information.

// Created by almess1 on 4/21/17.
//

#ifndef PUBGEO_IMAGE_H
#define PUBGEO_IMAGE_H

#include <cstring>
#include <thread>
#include <vector>
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

        // Default constructor
        Image() : width(0), height(0), bands(0), scale(1), offset(0), data(nullptr) {}

        // Copy constructor
        Image(const Image<TYPE>& i) :
            width(i.width),
            height(i.height),
            bands(i.bands),
            scale(i.scale),
            offset(i.offset),
            data(nullptr) {
            // Deep copy of data
            // Note: not using Allocate so we can skip the memset
            data = new TYPE *[height];
            for (unsigned int y = 0; y < height; y++) {
                data[y] = new TYPE[width * bands];
                std::memcpy(data[y], i.data[y], width * bands * sizeof(TYPE));
            }
        }

        // Move constructor
        Image(Image<TYPE>&& i) : Image() { // initialize via default constructor
            swap(*this, i);
        }

        // Destructor
        virtual ~Image() {
            Deallocate();
        }

        // Swap function
        friend void swap(Image<TYPE>& first, Image<TYPE>& second) {
            using std::swap; // Allow consideration of std::swap for following calls

            swap(first.width, second.width);
            swap(first.height, second.height);
            swap(first.bands, second.bands);
            swap(first.scale, second.scale);
            swap(first.offset, second.offset);
            swap(first.data, second.data);
        }

        // Assignment operator
        Image<TYPE>& operator=(Image<TYPE> rhs) {
            swap(*this, rhs);
            return *this;
        }

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

        /**
         * Apply a filter to this image, copying off a reference version for use by the filter algorithm
         * and overwriting this image with the new results
         *
         * filterFunc   A lambda / function pointer / functional that produces a value for a single
         *  cell in the output image, and takes the following arguments:
         *      FilterFunction(TYPE* val, const TYPE& ref, std::vector<TYPE> &ngbrs)
         *  where val is the output cell value, ref is the original value of the cell, and ngbrs are the
         *  neighboring values of the cell, including itself.
         *
         * rad  The distance in each direction that neighbors should be pulled from;
         *  e.g. rad=0 would only return 1 neighbor (itself), rad=1 would return 9, rad=2 would return 25, etc.
         *
         * voidVal  The image value that designates a void
         *
         * skipVoids    If true, cells that are void do not have the filter applied, and are also not
         *  inserted into the neighbors list
         */
        template <class FilterFunction>
        void filter(FilterFunction filterFunc, int rad = 1, TYPE voidVal = 0, bool skipVoids = true) {
            Image<TYPE> ref(*this);
            Image<TYPE>::filter(this, &ref, filterFunc, rad, voidVal, skipVoids);
        }

        /**
         * Apply a filter to the source image, storing results in the destination image
         *
         * filterFunc   See definition in non-static version of this function
         * dest Destination image; where results are stored
         * src  Source image; where reference and neighbor values are pulled from
         * rad  See definition in non-static version of this function
         * voidVal  See definition in non-static version of this function
         * skipVoids    See definition in non-static version of this function
         */
        template <class FilterFunction, class DEST_TYPE>
        static void filter(Image<DEST_TYPE>* dest, const Image<TYPE>* src, FilterFunction filterFunc, int rad = 1, TYPE voidVal = 0, bool skipVoids = true) {
            // TODO: Double check that src and dest are the same size

            std::vector<std::thread> workers;
            unsigned int N = std::min(std::thread::hardware_concurrency(),src->height);
            for (unsigned int k = 0; k < N; ++k) {
                workers.push_back(std::thread([=](){

                    // Create vector for neighbors
                    std::vector<TYPE> ngbrs;
                    ngbrs.reserve((2*rad+1)*(2*rad+1));

                    for (unsigned int j = k; j < src->height; j+=N) {
                        // Define row bounds:
                        unsigned int j1 = std::max((int) j - rad, 0);
                        unsigned int j2 = std::min(j + rad, src->height - 1);

                        for (unsigned int i = 0; i < src->width; i++) {
                            // Skip if void.
                            if (skipVoids && (src->data[j][i] == voidVal)) continue;

                            // Define bounds;
                            unsigned int i1 = std::max((int) i - rad, 0);
                            unsigned int i2 = std::min(i + rad, src->width - 1);

                            // Add valid values to the list.
                            for (unsigned int jj = j1; jj <= j2; jj++) {
                                for (unsigned int ii = i1; ii <= i2; ii++) {
                                    if (!skipVoids || (src->data[jj][ii] != voidVal)) {
                                        ngbrs.push_back(src->data[jj][ii]);
                                    }
                                }
                            }

                            // Apply filter
                            filterFunc(&(dest->data[j][i]), src->data[j][i], ngbrs);

                            // Clear neighbors list
                            ngbrs.clear();
                        }
                    }
                }));
            }

            for (unsigned int k = 0; k < N; ++k) {
                workers.at(k).join();
            }
        }

        /**
         * Apply a filter to two source images, storing results in the destination image
         *
         * dest Destination image; where results are stored
         * A    First source image; where reference and neighbor values are pulled from
         * B    Second source image; where reference and neighbor values are pulled from
         * selectFunc   Function that determines whether the filter function should be applied
         *  based soley on the reference values from images A and B.  This function must satisfy
         *  the following template: bool(std::pair<TYPE_A, TYPE_B>)
         * filterFunc   Function that operates on reference and neighbor values from images A and
         *  B and computes a single filtered value.  This function must satisfy the following
         *  template: void(TYPE&, const std::pair<TYPE_A, TYPE_B>&, std::vector<std::pair<TYPE_A, TYPE_B>>&)
         *  where the first argument is the destination value, the second is the reference values from
         *  A and B, and the third argument is the list of neighbor values from A and B
         * rad  The distance in each direction that neighbors should be pulled from;
         *  e.g. rad=0 would only return 1 neighbor (itself), rad=1 would return 9, rad=2 would return 25, etc.
         */
        template <class SelectionFunction, class ImFilterFunction, class TYPE_A, class TYPE_B>
        static void filter2(Image<TYPE>* dest, const Image<TYPE_A>* A, const Image<TYPE_B>* B,
                SelectionFunction selectFunc, ImFilterFunction filterFunc, int rad = 1) {
            std::vector<std::thread> workers;
            unsigned int N = std::min(std::thread::hardware_concurrency(), dest->height);
            for (unsigned int k = 0; k < N; ++k) {
                workers.push_back(std::thread([=](){

                    // Create vector for neighbors
                    std::vector<std::pair<TYPE_A,TYPE_B>> ngbrs;
                    ngbrs.reserve((2*rad+1)*(2*rad+1));

                    for (unsigned int j = k; j < dest->height; j+=N) {
                        // Define row bounds:
                        unsigned int j1 = std::max((int) j - rad, 0);
                        unsigned int j2 = std::min(j + rad, dest->height - 1);

                        for (unsigned int i = 0; i < dest->width; i++) {
                            // Define column bounds:
                            unsigned int i1 = std::max((int) i - rad, 0);
                            unsigned int i2 = std::min(i + rad, dest->width - 1);

                            std::pair<TYPE_A,TYPE_B> ref = std::make_pair(A->data[j][i],B->data[j][i]);
                            if (selectFunc(ref)) {

                                // Add valid values to the list.
                                for (unsigned int jj = j1; jj <= j2; jj++) {
                                    for (unsigned int ii = i1; ii <= i2; ii++) {
                                        ngbrs.emplace_back(A->data[jj][ii],B->data[jj][ii]);
                                    }
                                }

                                // Apply filter
                                filterFunc(&(dest->data[j][i]), ref, ngbrs);

                                // Clear neighbors list
                                ngbrs.clear();
                            }
                        }
                    }
                }));
            }

            for (unsigned int k = 0; k < N; ++k) {
                workers.at(k).join();
            }
        }
    };
}
#endif //PUBGEO_IMAGE_H
