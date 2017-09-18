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

#include <algorithm>

#include "orthoimage.h"
#include "plugin.hpp"
#include "shr3d.h"

#include <pdal/pdal_macros.hpp>

namespace pdal
{

static PluginInfo const s_info = PluginInfo(
    "writers.shr3d", "Shareable High Resolution 3D (Almes et al., 2017)",
    "http://pdal.io/stages/writers.shr3d.html");

CREATE_SHARED_PLUGIN(1, 0, Shr3dWriter, Writer, s_info)

std::string Shr3dWriter::getName() const
{
    return s_info.name;
}

void Shr3dWriter::addArgs(ProgramArgs& args)
{
    args.add("filename", "Output filename", m_filename).setPositional();
    args.add("dh", "Horizontal uncertainty", m_dh, 0.5);
    args.add("dz", "Vertical uncertainty", m_dz, 0.5);
    args.add("agl", "Minimum building height above ground level", m_agl, 2.0);
    args.add("area", "Minimum building area", m_area, 50.0);
    args.add("egm96", "Set vertical datum to EGM96", m_egm96, false);
}

void Shr3dWriter::write(const PointViewPtr view)
{
    shr3d::OrthoImage<unsigned short> dsmImage;
    if (!dsmImage.readFromPointView(view, m_dh, shr3d::MAX_VALUE))
        throw pdal_error("Error createing DSM\n");
    dsmImage.medianFilter(1, static_cast<unsigned int>(m_agl / dsmImage.scale));
    dsmImage.fillVoidsPyramid(true, 2);

    shr3d::OrthoImage<unsigned short> minImage;
    if (!minImage.readFromPointView(view, m_dh, shr3d::MIN_VALUE))
        throw pdal_error("Error creating minimum Z image\n");
    minImage.medianFilter(1, static_cast<unsigned int>(m_agl / minImage.scale));
    minImage.fillVoidsPyramid(true, 2);

    for (unsigned int j = 0; j < dsmImage.height; ++j)
    {
        for (unsigned int i = 0; i < dsmImage.width; ++i)
        {
            float minValue = static_cast<float>(minImage.data[j][i]);

            // This check is to avoid penalizing spurious returns under very
            // tall buildings. CAUTION: This is a hack to address an observed
            // lidar sensor issue and may not generalize well.
            if ((static_cast<float>(dsmImage.data[j][i]) - minValue) <
                (40.0 / minImage.scale))
            {
                // If this is in the trees, then set to void.
                bool found = false;
                double threshold = m_dz / dsmImage.scale;
                unsigned int i1 = std::max(0u, i - 1u);
                unsigned int i2 = std::min(i + 1u, dsmImage.width - 1u);
                unsigned int j1 = std::max(0u, j - 1u);
                unsigned int j2 = std::min(j + 1u, dsmImage.height - 1u);
                for (unsigned int jj = j1; jj <= j2; ++jj)
                {
                    for (unsigned int ii = i1; ii <= i2; ++ii)
                    {
                        float diff = static_cast<float>(dsmImage.data[jj][ii]) -
                                     minValue;
                        if (diff < threshold)
                            found = true;
                    }
                }
                if (!found)
                    dsmImage.data[j][i] = 0;
            }
        }
    }

    // Convert horizontal and vertical uncertainty values to bin units.
    int dh_bins = std::max(1, static_cast<int>(floor(m_dh / dsmImage.gsd)));
    unsigned int dz_short = static_cast<unsigned int>(m_dz / dsmImage.scale);
    unsigned int agl_short = static_cast<unsigned int>(m_agl / dsmImage.scale);

    // Generate label image.
    shr3d::OrthoImage<unsigned long> labelImage;
    labelImage.Allocate(dsmImage.width, dsmImage.height);
    labelImage.easting = dsmImage.easting;
    labelImage.northing = dsmImage.northing;
    labelImage.zone = dsmImage.zone;
    labelImage.gsd = dsmImage.gsd;

    // Allocate a DTM image as SHORT and copy in the DSM values.
    shr3d::OrthoImage<unsigned short> dtmImage;
    dtmImage.Allocate(dsmImage.width, dsmImage.height);
    dtmImage.easting = dsmImage.easting;
    dtmImage.northing = dsmImage.northing;
    dtmImage.zone = dsmImage.zone;
    dtmImage.gsd = dsmImage.gsd;
    dtmImage.scale = dsmImage.scale;
    dtmImage.offset = dsmImage.offset;
    dtmImage.bands = dsmImage.bands;
    for (unsigned int j = 0; j < dsmImage.height; ++j)
    {
        for (unsigned int i = 0; i < dsmImage.width; ++i)
        {
            dtmImage.data[j][i] = minImage.data[j][i];
        }
    }

    // Classify ground points.
    shr3d::Shr3dder::classifyGround(labelImage, dsmImage, dtmImage, dh_bins,
                                    dz_short);

    // For DSM voids, also set DTM value to void.
    for (unsigned int j = 0; j < dsmImage.height; ++j)
    {
        for (unsigned int i = 0; i < dsmImage.width; ++i)
        {
            if (dsmImage.data[j][i] == 0)
                dtmImage.data[j][i] = 0;
        }
    }

    // Median filter, replacing only points differing by more than the DZ
    // threshold.
    dtmImage.medianFilter(1, static_cast<unsigned int>(m_dz / dsmImage.scale));

    // Refine the object label image and export building outlines.
    shr3d::Shr3dder::classifyNonGround(dsmImage, dtmImage, labelImage, dz_short,
                                       agl_short, static_cast<float>(m_area));

    // Fill small voids in the DTM after all processing is complete.
    dtmImage.fillVoidsPyramid(true, 2);

    // Write the DTM image as FLOAT.
    dtmImage.write(const_cast<char*>(m_filename.c_str()), true, m_egm96);
}

} // namespace pdal
