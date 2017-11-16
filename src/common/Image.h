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
    typedef std::pair<size_t,size_t> Pixel; // Pixel coordinate type

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

        // Alternate constructor (creates identically sized, blank image)
        template<class OTYPE>
        Image(const Image<OTYPE>* i, unsigned int nbands = 1) :
                width(0),
                height(0),
                bands(0),
                scale(i->scale),
                offset(i->offset),
                data(nullptr) {
            Allocate(i->width,i->height,nbands);
        }

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

        /**
         * Traces all object boundaries in image
         *
         * Assumes image is a label image, as it uses equality to determine what pixels are grouped.
         *
         * All objects with non-zero labels will be traced.
         *
         * See traceBounds for details on algorithm.
         *
         * Output is map of closed, CW boundaries indexed by label.
         */
        std::map<TYPE,std::vector<Pixel>> traceBoundaries() const {
            std::map<TYPE,std::vector<Pixel>> boundaries;

            for (size_t j = 0; j < height; ++j) {
                for (size_t i = 0; i < width; ++i) {
                    const TYPE& v = data[j][i];
                    if (v && !boundaries.count(v)) {
                        boundaries.emplace(v,traceBounds(j,i));
                    }
                }
            }

            return boundaries;
        }

        /**
         * Traces object boundary starting at a location on its top
         *
         * Assumes image is a label image, as it uses equality to determine what pixels are grouped.
         *
         * Function will return a closed, CW boundary of object.  If object is foreground (pixel
         * value > 0), the algorithm will use 8-connectivity; if it's background, algorithm will use
         * 4-connectivity.
         *
         * Behavior is undefined if starting location is not on the top boundary of the object (i.e.
         * there is no row < r that contains a pixel with the same label).
         *
         * Inputs are the starting pixel row & column, and the output is a vector of the traced coordinates.
         */
        std::vector<Pixel> traceBounds(size_t r, size_t c) const {
            const TYPE& v = data[r][c]; // Label at location
            std::vector<Pixel> b; // Structure that will contain pixel bounds

            const std::array<int,8> dj = {1, 1, 0,-1,-1,-1, 0, 1}; // Delta row, indexed by direction
            const std::array<int,8> di = {0,-1,-1,-1, 0, 1, 1, 1}; // Delta column, indexed by direction

            int stride = v > 0 ? 1 : 2; // For foreground pixels, we'll have 8-connectivity, and 4-connectivity for background

            size_t m = r; // Current pixel row
            size_t n = c; // Current pixel column
            int first_dir = -1; // First direction on boundary, initialize invalid
            int last_dir = 0;   // Direction traveled from previous pixel to current pixel
            int new_dir; // Direction from current pixel to next pixel
            int fin_dir; // Last direction to check

            while (true) {
                b.emplace_back(m,n); // Add pixel location
                fin_dir = (last_dir+4) % 8; // Final direction reflects back to the pixel we just came from

                // Iterate through directions, first looking to the left relative to the last direction,
                // and then progressively looking to the right
                for (new_dir = (last_dir+6) % 8; new_dir != fin_dir; new_dir = (new_dir+stride) % 8) {
                    size_t p = m+dj[new_dir]; // Test pixel row
                    size_t q = n+di[new_dir]; // Test pixel column
                    if (!(p < 0 || p >= height || q < 0 || q >= width) && (data[p][q] == v))
                        break; // If pixel coordinate is valid and labeled correctly, break
                }

                if ((new_dir == fin_dir) && (b.size()==1)) { // If there's only one pixel with this label
                    b.emplace_back(m,n); // Close polygon
                    break;
                } else if (m==r && n==c && new_dir == first_dir) { // If we're back where we started
                    break;
                } else if (first_dir < 0) { // If this is the first segment, initialize first_dir
                    first_dir = new_dir;
                }

                // Update state
                m+=dj[new_dir];
                n+=di[new_dir];
                last_dir = new_dir;
            }

            return b;
        }
    };
}
#endif //PUBGEO_IMAGE_H
