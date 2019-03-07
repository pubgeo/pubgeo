// Copyright 2017 The Johns Hopkins University Applied Physics Laboratory.
// Licensed under the MIT License. See LICENSE.txt in the project root for full license information.

// Created by almess1 on 4/21/17.
//

#ifndef PUBGEO_IMAGE_H
#define PUBGEO_IMAGE_H

#include <cstring>
#include <set>
#include <thread>
#include <unordered_set>
#include <vector>

#include "disjoint_set.h"

namespace pubgeo {
    template<class TYPE>
    class Image {

    public:
        unsigned int width;
        unsigned int height;
        unsigned int bands;
        float scale;
        float offset;
        std::vector<std::vector<TYPE>> data;    // Order is HEIGHT (outside), WIDTH, BAND (inside).

        // Default constructor
        Image() : width(0), height(0), bands(1), scale(1), offset(0), data() {}

        // Alternate constructor (creates identically sized, blank image)
        template<class OTYPE>
        Image(const Image<OTYPE>* i, unsigned int nbands = 1) :
                width(i->width),
                height(i->height),
                bands(nbands),
                scale(i->scale),
                offset(i->offset),
                data(height,std::vector<TYPE>(width*bands,0)) {}

        // Destructor
        virtual ~Image() {}

        // Allocate memory.
        void Allocate(unsigned int numColumns, unsigned int numRows, unsigned int numBands = 1) {
            Deallocate();
            width = numColumns;
            height = numRows;
            bands = numBands;
            data.resize(height,std::vector<TYPE>(width*bands,0));
        }

        void Deallocate() {
            width = 0;
            height = 0;
            bands = 1;
            data.clear();
        }

        // Test if object does not contain data
        bool empty() {
            return data.empty();
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

        /**
         * Labels connected components and interior holes
         *
         * Assumes src_im is a binary image, where 0 labels background and !0 labels foreground pixels
         *
         * label_im will contain positive labels for foreground objects, and negative labels for background
         * objects.  Any background that is connected to the image edge will be labeled 0.
         *
         * Functionality is loosely based on bwconncomp, but algorithm is drawn from:
         * https://en.wikipedia.org/wiki/Connected-component_labeling#Two-pass
         */
        template <class OTYPE>
        void labelConnectedComponentsFrom(const Image<OTYPE>* src_im) {
            size_t M = src_im->height;
            size_t N = src_im->width;

            // Initial labeling
            std::set<size_t> ngbrs;
            DisjointSet fg_labels;
            DisjointSet bg_labels;
            bg_labels.add(); // Add exterior background

            for (size_t j = 0; j < M; ++j) {
                for (size_t i = 0; i < N; ++i) {
                    ngbrs.clear();
                    DisjointSet* labels;

                    // Find neighbors
                    if (src_im->data[j][i]) {
                        // Foreground
                        labels = &fg_labels;

                        // Find neighbors (8-connectivity)
                        if (j) {
                            if (i && src_im->data[j-1][i-1])
                                ngbrs.insert(data[j-1][i-1]);
                            if (src_im->data[j-1][i])
                                ngbrs.insert(data[j-1][i]);
                            if (i<N-1 && src_im->data[j-1][i+1])
                                ngbrs.insert(data[j-1][i+1]);
                        }
                        if (i && src_im->data[j][i-1])
                            ngbrs.insert(data[j][i-1]);
                    } else {
                        // Background
                        labels = &bg_labels;

                        // Find neighbors (4-connectivity)
                        if (!j || j==M-1 || !i || i==N-1) // If we're on the edge
                            ngbrs.insert(0);
                        if (j && !src_im->data[j-1][i])
                            ngbrs.insert(data[j-1][i]);
                        if (i && !src_im->data[j][i-1])
                            ngbrs.insert(data[j][i-1]);
                    }

                    if (ngbrs.empty()) {
                        // Assign new label
                        data[j][i] = (TYPE) labels->add();
                    } else {
                        // Assign pixel the smallest label
                        size_t lbl = *ngbrs.begin();
                        data[j][i] = (TYPE) lbl;

                        // Link neighbors
                        for (size_t l : ngbrs)
                            labels->merge(l,lbl);
                    }
                }
            }

            // Simplify labels
            std::vector<TYPE> bg_final_labels = bg_labels.flatten<TYPE>();
            std::vector<TYPE> fg_final_labels = fg_labels.flatten<TYPE>(1);

            // Pass 2
            for (size_t j = 0; j < M; ++j) {
                for (size_t i = 0; i < N; ++i) {
                    TYPE &lbl = data[j][i];
                    lbl = src_im->data[j][i] ? fg_final_labels[lbl] : -bg_final_labels[lbl];
                }
            }
        }

        // Sets image to the src image upsampled by scale_factor, using nearest-neighbor algorithm
        virtual void nn_upsample(const Image<TYPE>* src, unsigned int scale_factor) {
            width = src->width*scale_factor;
            height = src->height*scale_factor;
            bands = src->bands;
            scale = src->scale;
            offset = src->offset;

            // First, allocate
            for (size_t y = 0; y < data.size(); ++y) {
                data[y].resize(width * bands);
            }
            data.resize(height,std::vector<TYPE>(width * bands));

            // Second, copy
            size_t y, x;
            for (size_t sy = 0; sy < src->height; ++sy) {
                y = sy*scale_factor;
                for (size_t sx = 0; sx < src->width * bands; ++sx) {
                    x = sx*scale_factor;
                    for (size_t i = 0; i < scale_factor; ++i)
                        data[y][x+i] = src->data[sy][sx]; // Duplicate element along columns
                }

                // Duplicate rows
                for (size_t j = 1; j < scale_factor; ++j)
                    data[y+j] = data[y];
            }
        }
    };
}
#endif //PUBGEO_IMAGE_H
