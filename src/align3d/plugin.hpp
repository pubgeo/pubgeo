// Copyright (c) 2017, Bradley J Chambers (brad.chambers@gmail.com)
// All rights reserved.
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

// PDAL Plugin of S. Almes, S. Hagstrom, D. Chilcott, H. Goldberg, M. Brown,
// “Open Source Geospatial Tools to Enable Large Scale 3D Scene Modeling,”
// FOSS4G, 2017.

#pragma once

#include <cstdint>
#include <string>

#include <pdal/Filter.hpp>
#include <pdal/pdal_export.hpp>
#include <pdal/util/ProgramArgs.hpp>

namespace pdal
{

class PDAL_DLL Align3dFilter : public Filter
{
public:
    Align3dFilter() : Filter(), m_fixed(nullptr)
    {
    }

    static void* create();
    static int32_t destroy(void*);
    std::string getName() const;

private:
    virtual void addArgs(ProgramArgs& args);
    virtual PointViewSet run(PointViewPtr view);

    double m_gsd;
    double m_maxt;
    double m_maxdz;
    PointViewPtr m_fixed;

    Align3dFilter& operator=(const Align3dFilter&) = delete;
    Align3dFilter(const Align3dFilter&) = delete;
};
}
