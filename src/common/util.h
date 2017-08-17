// Copyright 2017 The Johns Hopkins University Applied Physics Laboratory.
// Licensed under the MIT License. See LICENSE.txt in the project root for full license information.

//
// Created by almess1 on 4/20/17.
//

#ifndef PUBGEO_UTIL_H
#define PUBGEO_UTIL_H

#include <limits>

namespace pubgeo {
    template<typename T>
    inline T &Max(const T &a, const T &b) {
        return a < b ? b : a;
    }

    template<typename T>
    inline T &Min(const T &a, const T &b) {
        return a > b ? b : a;
    }

    const int MAX_INT = std::numeric_limits<int>::max();
    const float MAX_FLOAT = std::numeric_limits<float>::max();
}
#endif //PUBGEO_UTIL_H
