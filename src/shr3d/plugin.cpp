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
    args.add("dh", "Horizontal uncertainty", shr3dder.dh_meters, shr3dder.dh_meters);
    args.add("dz", "Vertical uncertainty", shr3dder.dz_meters, shr3dder.dz_meters);
    args.add("agl", "Minimum building height above ground level", shr3dder.agl_meters, shr3dder.agl_meters);
    args.add("area", "Minimum building area", shr3dder.min_area_meters, shr3dder.min_area_meters);
    args.add("egm96", "Set vertical datum to EGM96", shr3dder.egm96, shr3dder.egm96);
    args.add("image_type", "Image Type (0=DSM,1=MIN,2=DTM,3=CLASS,4=BUILDING)", output, output);
}

void Shr3dWriter::write(const PointViewPtr view)
{
    if (output < shr3d::DSM || output > shr3d::LABEL)
        throw pdal_error("Unknown output image type\n");
    std::map<shr3d::ImageType,std::string> outputFilenames;
    outputFilenames[(shr3d::ImageType) output] = m_filename;

    // Read in point cloud
    shr3d::PointCloud pset;
    if (!pset.Read(view))
        throw pdal_error("Error reading from view\n");

    // Create DSM
    shr3d::OrthoImage<unsigned short> dsmImage;
    if (!shr3dder.createDSM(pset, dsmImage))
        throw pdal_error("Error creating DSM\n");

    // Create MIN
    shr3d::OrthoImage<unsigned short> minImage;
    if (!shr3dder.createMIN(pset,minImage))
        throw pdal_error("Error creating minimum Z image\n");

    // Create everything else
    shr3dder.process(dsmImage, minImage, outputFilenames);
}

} // namespace pdal
