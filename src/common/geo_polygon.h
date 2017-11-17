/*
 * geo_polygon.h
 *
 *  Created on: Nov 17, 2017
 *      Author: fosterkh
 *
 * Copyright 2017 The Johns Hopkins University Applied Physics Laboratory.
 * Licensed under the MIT License. See LICENSE.txt in the project root for full license information.
 */

#ifndef PUBGEO_GEO_POLYGON_H_
#define PUBGEO_GEO_POLYGON_H_

#include <array>
#include <vector>

#include <ogrsf_frmts.h>

#include "orthoimage.h"

namespace pubgeo {

/**
 * Implements a geo-referenced polygon, with optional holes
 *
 * Points are stored relative to a given easting and northing (and zone),
 * with scale factors so integer types can be used for point coordinates.
 */
template <class TYPE = size_t>
class GeoPolygon {
public:
    typedef std::array<TYPE,2> Point;
    std::vector<Point> ring;
    std::vector<GeoPolygon> inner_rings;

    double easting;
    double northing;
    int zone;
    double xscale;
    double yscale;

    // Default constructor
    GeoPolygon() : easting(0), northing(0), zone(0), xscale(1), yscale(1) {}

    /**
     * Creates a geo-polygon by tracing object boundary starting at a location on its top
     *
     * Assumes input image is a label image, as it uses equality to determine what pixels are grouped.
     *
     * Function will create a closed, CW boundary of object.  If object is foreground (pixel
     * value > 0), the algorithm will use 8-connectivity; if it's background, algorithm will use
     * 4-connectivity.
     *
     * Behavior is undefined if starting location is not on the top boundary of the object (i.e.
     * there is no row < r that contains a pixel with the same label).
     *
     * Inputs are the starting pixel row & column, and the output is a vector of the traced coordinates.
     */
    template <class LBL_TYPE>
    GeoPolygon(const OrthoImage<LBL_TYPE>& labels, size_t r, size_t c) :
            ring(),
            inner_rings(),
            easting(labels.easting+labels.gsd*0.5),
            northing(labels.northing+labels.gsd*(labels.height-0.5)),
            zone(labels.zone),
            xscale(labels.gsd),
            yscale(-labels.gsd) {
        const LBL_TYPE& v = labels.data[r][c]; // Label at start location

        const std::array<int,8> dj = {1, 1, 0,-1,-1,-1, 0, 1}; // Delta row, indexed by direction
        const std::array<int,8> di = {0,-1,-1,-1, 0, 1, 1, 1}; // Delta column, indexed by direction

        int stride = v > 0 ? 1 : 2; // For foreground pixels, we'll have 8-connectivity, and 4-connectivity for background

        TYPE m = r; // Current pixel row
        TYPE n = c; // Current pixel column
        int first_dir = -1; // First direction on boundary, initialize invalid
        int last_dir = 0;   // Direction traveled from previous pixel to current pixel
        int new_dir; // Direction from current pixel to next pixel
        int fin_dir; // Last direction to check

        while (true) {
            ring.push_back({n,m}); // Add pixel location
            fin_dir = (last_dir+4) % 8; // Final direction reflects back to the pixel we just came from

            // Iterate through directions, first looking to the left relative to the last direction,
            // and then progressively looking to the right
            for (new_dir = (last_dir+6) % 8; new_dir != fin_dir; new_dir = (new_dir+stride) % 8) {
                TYPE p = m+dj[new_dir]; // Test pixel row
                TYPE q = n+di[new_dir]; // Test pixel column
                if (!(p < 0 || p >= labels.height || q < 0 || q >= labels.width) && (labels.data[p][q] == v))
                    break; // If pixel coordinate is valid and labeled correctly, break
            }

            if ((new_dir == fin_dir) && (ring.size()==1)) { // If there's only one pixel with this label
                ring.push_back({n,m}); // Close polygon
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
    }

    /**
     * Traces all object boundaries in image
     *
     * Assumes input image is a label image, as it uses equality to determine what pixels are grouped.
     *
     * All objects with non-zero labels will be traced; holes (objects with labels < 0) will be
     * assigned as an inner ring of their exterior polygon.
     *
     * See GeoPolygon constructor for details on algorithm.
     *
     * Output is map of closed, CW boundaries indexed by label.
     */
    template <class LBL_TYPE>
    static std::map<LBL_TYPE,GeoPolygon> traceBoundaries(const OrthoImage<LBL_TYPE>& labels) {
        std::map<LBL_TYPE,GeoPolygon> boundaries;

        // Create all polygons (exterior and interior)
        for (size_t j = 0; j < labels.height; ++j) {
            for (size_t i = 0; i < labels.width; ++i) {
                const LBL_TYPE& v = labels.data[j][i];
                if (v && !boundaries.count(v))
                    boundaries.emplace(v,GeoPolygon(labels,j,i));
            }
        }

        // Assign inner rings to exterior polygons
        auto end = boundaries.lower_bound(0);
        for (auto iter = boundaries.begin(); iter!=end; boundaries.erase(iter++)) { // Remove after incrementing iterator
            const Point& pt = iter->second.ring.front();
            LBL_TYPE ext_label = labels.data[pt[1]-1][pt[0]];
            boundaries[ext_label].inner_rings.push_back(iter->second);
        }

        return boundaries;
    }

    // Save collection of polygons to a shapefile
    template <class LBL_TYPE>
    static bool write(const std::string& filename, const std::map<LBL_TYPE,GeoPolygon>& bounds) {
        // Need to delete file if it already exists
        if (FILE *file = fopen(filename.c_str(), "r")) {
            fclose(file);
            std::remove(filename.c_str());
        }

        // Don't create a file if it's going to be empty
        if (bounds.empty())
            return false;

        // GDAL setup
        OGRRegisterAll();
        OGRSFDriver* poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("ESRI Shapefile");
        if(poDriver == NULL) {
            printf("ERROR: ESRI Shapefile driver not available.\n");
            return false;
        }

        OGRDataSource* poDS = poDriver->CreateDataSource(filename.c_str(), NULL);
        if(poDS == NULL) {
            printf("ERROR: Creation of output shapefile '%s' failed.\n", filename.c_str());
            return false;
        }

        // Set coordinate system
        int zone = bounds.empty() ? 1 : bounds.begin()->second.zone;
        OGRSpatialReference oSRS;
        oSRS.SetProjCS("UTM (WGS84)");
        oSRS.SetUTM(abs(zone), (zone > 0));
        int theZoneIsNorthern = 255;
        oSRS.GetUTMZone(&theZoneIsNorthern);
        oSRS.SetWellKnownGeogCS("WGS84");

        OGRLayer* poLayer = poDS->CreateLayer("buildings", &oSRS, wkbMultiPolygon, NULL);
        if(poLayer == NULL) {
            printf("ERROR: Layer creation failed.\n");
            return false;
        }

        OGRFieldDefn oField("Label", OFTInteger);
        if(poLayer->CreateField(&oField) != OGRERR_NONE) {
            printf("ERROR: Creating label field failed.\n");
            return false;
        }

        // Populate shapefile with data
        for (auto iter = bounds.lower_bound(0); iter != bounds.end(); ++iter) {
            OGRFeature* poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
            poFeature->SetField("Label", iter->first);

            // Add geometry
            OGRPolygon polygon = iter->second.makeOGRPolygon();
            poFeature->SetGeometry(&polygon);

            // Add feature to file
            if(poLayer->CreateFeature(poFeature) != OGRERR_NONE) {
                printf("ERROR: Failed to create feature in shapefile.\n");
                return false;
            }
            OGRFeature::DestroyFeature(poFeature);
        }

        GDALClose(poDS);
        return true;
    }

    // Convert polygon to an OGRPolygon
    OGRPolygon makeOGRPolygon() const {
        OGRPolygon polygon;
        OGRLinearRing ogrRing = makeOGRLinearRing();
        polygon.addRing(&ogrRing);

        for (const GeoPolygon& innerRing : inner_rings) {
            OGRLinearRing ogrRing2 = innerRing.makeOGRLinearRing();
            polygon.addRing(&ogrRing2);
        }

        return polygon;
    }

private:
    // Convert points to OGRLinearRing
    OGRLinearRing makeOGRLinearRing() const {
        OGRLinearRing ogrRing;
        for (const Point& pt : ring)
            ogrRing.addPoint(easting+xscale*pt[0],northing+yscale*pt[1]);
        return ogrRing;
    }
};

} // pubgeo

#endif /* PUBGEO_GEO_POLYGON_H_ */
