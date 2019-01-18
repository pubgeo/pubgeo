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
// "Open Source Geospatial Tools to Enable Large Scale 3D Scene Modeling,"
// FOSS4G, 2017.

#include <algorithm>
#include <iterator>

#include <boost/algorithm/string/trim.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include "orthoimage.h"
#include "plugin.hpp"
#include "shr3d.h"

namespace pdal
{

static StaticPluginInfo const s_info = {
    "writers.shr3d",
    "Shareable High Resolution 3D (Almes et al., 2017)",
    "http://pdal.io/stages/writers.shr3d.html"
};

CREATE_SHARED_STAGE(Shr3dWriter, s_info)

Shr3dWriter::Shr3dWriter() : shr3dder(), m_filename() {
    output_types = (boost::format("%d,%d,%d,%d,%d,%d,%d")
        % shr3d::DSM
        % shr3d::MIN
        % shr3d::DTM
        % shr3d::INTENSITY
        % shr3d::CLASS
        % shr3d::BUILDING
        % shr3d::BUILDING_OUTLINES).str();
}

std::string Shr3dWriter::getName() const
{
    return s_info.name;
}

void Shr3dWriter::addArgs(ProgramArgs& args)
{
    args.add("filename", "Output filename (if multiple image types are specified, will be used as a basename)", m_filename).setPositional();
    args.add("dh", "Horizontal uncertainty", shr3dder.dh_meters, shr3dder.dh_meters);
    args.add("dz", "Vertical uncertainty", shr3dder.dz_meters, shr3dder.dz_meters);
    args.add("agl", "Minimum building height above ground level", shr3dder.agl_meters, shr3dder.agl_meters);
    args.add("area", "Minimum building area", shr3dder.min_area_meters, shr3dder.min_area_meters);
    args.add("egm96", "Set vertical datum to EGM96", shr3dder.egm96, shr3dder.egm96);
    args.add("bounds", "Bounds for cropped image, e.g.: '([xmin, xmax], [ymin, ymax])'", shr3dder.bounds);
    args.add("image_type",
        (boost::format("Output image type (%d=DSM,%d=MIN,%d=DTM,%d=INTENSITY,%d=CLASS,%d=BUILDING,%d=BUILDING SHAPEFILE)")
            % shr3d::DSM % shr3d::MIN % shr3d::DTM % shr3d::INTENSITY % shr3d::CLASS % shr3d::BUILDING % shr3d::BUILDING_OUTLINES).str(),
        output_types, output_types);
}

std::vector<shr3d::ImageType> getImageTypeList(std::string str) {
    std::vector<shr3d::ImageType> output;

    //Remove whitespace from string
    boost::algorithm::trim(str);

    // Parse and convert
    boost::tokenizer<> tok(str);
    transform(tok.begin(), tok.end(), std::back_inserter(output),
        [](std::string s) -> shr3d::ImageType {
            return (shr3d::ImageType) boost::lexical_cast<unsigned int>(s);
        });
    return output;
}

void Shr3dWriter::write(const PointViewPtr view)
{
    // Set output filenames
    std::map<shr3d::ImageType,std::string> outputFilenames;
    
    std::vector<shr3d::ImageType> types = getImageTypeList(output_types);
    
    if (types.empty()) {
        throw pdal_error("No image type specified");
    } else if (types.size() == 1) {
        outputFilenames[types.front()] = m_filename;
    } else {
        for (shr3d::ImageType t : types) {
            switch (t) {
                case shr3d::DSM: outputFilenames[t] = m_filename + "_DSM.tif"; break;
                case shr3d::MIN: outputFilenames[t] = m_filename + "_MIN.tif"; break;
                case shr3d::DTM: outputFilenames[t] = m_filename + "_DTM.tif"; break;
                case shr3d::INTENSITY: outputFilenames[t] = m_filename + "_INT.tif"; break;
                case shr3d::CLASS: outputFilenames[t] = m_filename + "_class.tif"; break;
                case shr3d::BUILDING: outputFilenames[t] = m_filename + "_buildings.tif"; break;
                case shr3d::BUILDING_OUTLINES: outputFilenames[t] = m_filename + "_buildings.shp"; break;
                case shr3d::MINAGL: outputFilenames[t] = m_filename + "_MINAGL.tif"; break;
                default: throw pdal_error("Unknown image type");
            }
        }
    }

    // Read in point cloud
    if (!shr3dder.setPSET(view))
        throw pdal_error("Error reading from view\n");

    // Run process
    shr3dder.process(outputFilenames);

    // Reset for subsequent uses
    shr3dder.reset();
}

} // namespace pdal
