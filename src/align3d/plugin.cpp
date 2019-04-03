// Original work Copyright (c) 2017, Bradley J Chambers (brad.chambers@gmail.com)
// Modified work Copyright (c) 2018-2019, The Johns Hopkins University /
//   Applied Physics Laboratory (JHU/APL)
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

#include "align3d.h"
#include "orthoimage.h"
#include "plugin.hpp"

namespace pdal
{

static StaticPluginInfo const s_info = {
    "filters.align3d",
    "Shareable High Resolution 3D (Almes et al., 2017)",
    "http://pdal.io/stages/filters.align3d.html"
};

CREATE_SHARED_STAGE(Align3dFilter, s_info)

std::string Align3dFilter::getName() const
{
    return s_info.name;
}

void Align3dFilter::addArgs(ProgramArgs& args)
{
    args.add("gsd", "Ground Sample Distance (GSD) for gridding", m_gsd, 1.0);
    args.add("maxt", "Maximum XYZ translation in search", m_maxt, 10.0);
    args.add("maxdz", "Max local Z difference for matching", m_maxdz, 0.0);
}

PointViewSet Align3dFilter::run(PointViewPtr view)
{
    PointViewSet viewSet;
    if (!view->size())
        return viewSet;

    if (!m_fixed)
    {
        m_fixed = view;
        return viewSet;
    }

    // Default MAXDZ = GSD x 2 to ensure reliable performance on steep
    // slopes.
    if (m_maxdz == 0.0)
        m_maxdz = m_gsd * 2.0;

    // Align the target point cloud to the reference.
    // AlignTarget2Reference(referenceFileName, targetFileName, params);
    // Read the reference LAS file as a DSM.
    // Fill small voids.
    // Remove points along edges which are difficult to match.
    log()->get(LogLevel::Debug) << "Reading reference point cloud\n";
    pubgeo::OrthoImage<unsigned short> referenceDSM;
    if (!referenceDSM.readFromPointView(m_fixed, m_gsd, pubgeo::MAX_VALUE))
        throw pdal_error("Error reading reference PointView\n");
    referenceDSM.fillVoidsPyramid(true, 2);
    log()->get(LogLevel::Debug) << "Filtering reference point cloud.\n";
    referenceDSM.edgeFilter((long)(m_maxdz / referenceDSM.scale));

    // Read the target LAS file as a DSM.
    // Fill small voids.
    // Remove points along edges which are difficult to match.
    log()->get(LogLevel::Debug) << "Reading target point cloud\n";
    pubgeo::OrthoImage<unsigned short> targetDSM;
    if (!targetDSM.readFromPointView(view, m_gsd, pubgeo::MAX_VALUE))
        throw pdal_error("Error reading target PointView\n");
    targetDSM.fillVoidsPyramid(true, 2);
    log()->get(LogLevel::Debug) << "Filtering target point cloud.\n";
    targetDSM.edgeFilter((long)(m_maxdz / targetDSM.scale));

    // Get overlapping bounds.
    align3d::AlignBounds bounds;
    bounds.xmin = std::max(referenceDSM.easting, targetDSM.easting);
    bounds.ymin = std::max(referenceDSM.northing, targetDSM.northing);
    bounds.xmax =
        std::min(referenceDSM.easting + (referenceDSM.width * referenceDSM.gsd),
                 targetDSM.easting + (targetDSM.width * targetDSM.gsd));
    bounds.ymax = std::min(
        referenceDSM.northing + (referenceDSM.height * referenceDSM.gsd),
        targetDSM.northing + (targetDSM.height * targetDSM.gsd));
    bounds.width = bounds.xmax - bounds.xmin;
    bounds.height = bounds.ymax - bounds.ymin;
    double overlap_km = bounds.width / 1000.0 * bounds.height / 1000.0;
    log()->get(LogLevel::Debug)
        << "Overlap = " << bounds.width << " m x " << bounds.height
        << " m = " << overlap_km << " km\n";
    if (overlap_km == 0.0)
        throw pdal_error("PointViews do not overlap");

    // Estimate rigid body transform to align target points to reference.
    align3d::AlignResult result;
    log()->get(LogLevel::Debug) << "Estimating rigid body transformation.\n";
    EstimateRigidBody(referenceDSM, targetDSM, m_maxt, bounds, result);

    for (PointId i = 0; i < view->size(); ++i)
    {
        double x = view->getFieldAs<double>(Dimension::Id::X, i);
        double y = view->getFieldAs<double>(Dimension::Id::Y, i);
        double z = view->getFieldAs<double>(Dimension::Id::Z, i);
        view->setField(Dimension::Id::X, i, x + result.tx);
        view->setField(Dimension::Id::Y, i, y + result.ty);
        view->setField(Dimension::Id::Z, i, z + result.tz);
    }
    viewSet.insert(view);

    MetadataNode root = getMetadata();
    root.add("x_offset", result.tx);
    root.add("y_offset", result.ty);
    root.add("z_offset", result.tz);
    root.add("z_rms", result.rms);

    return viewSet;
}

} // namespace pdal
