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

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include <ogrsf_frmts.h>

#include "orthoimage.h"

namespace pubgeo {

/**
 * 2-D Representation of a point, with associated convenience functions
 */
template <class TYPE = size_t>
struct GeoPoint {
    TYPE x;
    TYPE y;

    GeoPoint() : x(0), y(0) {}
    GeoPoint(const TYPE& x, const TYPE& y) : x(x), y(y) {}

    // Provide for numeric indexing; caution, this does not check bounds
    TYPE& operator[](size_t idx) {
        return *((&x)+idx);
    }
    const TYPE operator[](size_t idx) const {
        return *((&x)+idx);
    }

    double abs() const {
        return sqrt((double) (x*x+y*y));
    }

    double dot(const GeoPoint& rhs) {
        return x*rhs.x + y*rhs.y;
    }

    double cross(const GeoPoint& rhs) {
        return x*rhs.y - y*rhs.x;
    }

    GeoPoint& operator+=(const GeoPoint& rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }
    GeoPoint& operator-=(const GeoPoint& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }
    GeoPoint& operator*=(const double& rhs) {
        x *= rhs;
        y *= rhs;
        return *this;
    }
    GeoPoint& operator/=(const double& rhs) {
        x /= rhs;
        y /= rhs;
        return *this;
    }
};
template <class T>
inline GeoPoint<T> operator+(GeoPoint<T> lhs, const GeoPoint<T>& rhs) { lhs += rhs; return lhs; }
template <class T>
inline GeoPoint<T> operator-(GeoPoint<T> lhs, const GeoPoint<T>& rhs) { lhs -= rhs; return lhs; }
template <class T>
inline GeoPoint<T> operator*(GeoPoint<T> lhs, const double& rhs) { lhs *= rhs; return lhs; }
template <class T>
inline GeoPoint<T> operator*(const double& lhs, const GeoPoint<T>& rhs) { return rhs * lhs; }
template <class T>
inline GeoPoint<T> operator/(GeoPoint<T> lhs, const double& rhs) { lhs /= rhs; return lhs; }
template <class T>
std::ostream& operator<<(std::ostream& os, const GeoPoint<T>& pt) { os << '(' << pt.x << ',' << pt.y << ')'; return os; }
template <class T>
inline bool operator==(const GeoPoint<T>& lhs, const GeoPoint<T>& rhs) { return lhs.x == rhs.x && lhs.y == rhs.y; }
template <class T>
inline bool operator!=(const GeoPoint<T>& lhs, const GeoPoint<T>& rhs) { return !operator==(lhs,rhs); }

/**
 * Implements a geo-referenced polygon, with optional holes
 *
 * Points are stored relative to a given easting and northing (and zone),
 * with scale factors so integer types can be used for point coordinates.
 */
template <class TYPE = size_t>
class GeoPolygon {
public:
    typedef GeoPoint<TYPE> Point;
    typedef GeoPoint<double> DblPoint;
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
     * Function will create an unclosed, CW boundary of object.  If object is foreground (pixel
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

        size_t m = r; // Current pixel row
        size_t n = c; // Current pixel column
        int first_dir = -1; // First direction on boundary, initialize invalid
        int last_dir = 0;   // Direction traveled from previous pixel to current pixel
        int new_dir; // Direction from current pixel to next pixel
        int fin_dir; // Last direction to check

        while (true) {
            ring.push_back({(TYPE) n, (TYPE) m}); // Add pixel location
            fin_dir = (last_dir+4) % 8; // Final direction reflects back to the pixel we just came from

            // Iterate through directions, first looking to the left relative to the last direction,
            // and then progressively looking to the right
            for (new_dir = (last_dir+6) % 8; new_dir != fin_dir; new_dir = (new_dir+stride) % 8) {
                size_t p = m+dj[new_dir]; // Test pixel row
                size_t q = n+di[new_dir]; // Test pixel column
                if (!(p < 0 || p >= labels.height || q < 0 || q >= labels.width) && (labels.data[p][q] == v))
                    break; // If pixel coordinate is valid and labeled correctly, break
            }

            if ((new_dir == fin_dir) && (ring.size()==1)) { // If there's only one pixel with this label
                return;
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

        // Unclose-polygon
        ring.pop_back();
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
        std::map<LBL_TYPE,LBL_TYPE> hole_mapping; // Mapping between interior holes to exterior polygons

        // Create all polygons (exterior and interior)
        for (size_t j = 0; j < labels.height; ++j) {
            for (size_t i = 0; i < labels.width; ++i) {
                const LBL_TYPE& v = labels.data[j][i];
                if (v && !boundaries.count(v)) {
                    boundaries.emplace(v,GeoPolygon(labels,j,i));
                    if (v < 0)
                        hole_mapping[v] = labels.data[j-1][i];
                }
            }
        }

        // Assign inner rings to exterior polygons
        for (const auto& mp : hole_mapping) {
            auto iter = boundaries.find(mp.first);
            boundaries[mp.second].inner_rings.push_back(iter->second);
            boundaries.erase(iter);
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
        GDALAllRegister();
        GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
        if(poDriver == NULL) {
            printf("ERROR: ESRI Shapefile driver not available.\n");
            return false;
        }

        GDALDataset* poDS = poDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL);
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

    /**
     * Applies several algorithms applicable to buildings to simplify polygon
     *
     * SCALE parameter controls the size of features to retain; units are pixels
     */
    GeoPolygon<double> buildingSimplify(unsigned int scale = 36) const {
        // Initialize new polygon
        GeoPolygon<double> poly;
        poly.easting = easting;
        poly.northing = northing;
        poly.zone = zone;
        poly.xscale = xscale;
        poly.yscale = yscale;

        // Find local corners in the polygon
        std::list<size_t> indices = findCorners(scale);
        if (indices.size() < 3) return poly; // Insufficient # of points to proceed


        // Add points we might have missed because of gently curving lines, and because
        // of the omission we would be missing significant portions of the building
        addMissingPoints(indices, 0.5 * scale);

        // Adjust point locations to fractions of a pixel, and add points as necessary to keep all
        // angles greater than 60 degrees
        std::list<DblPoint> points = improvePoints(indices, scale / 3);

        // Transfer points to new polygon
        for (DblPoint p : points)
            poly.ring.push_back(p);

        // Remove points that add miniscule amounts of area to object
        // This function also will close the polygon.
        poly.vwReduce(scale * 0.25);

        // Recursively apply to sub polygons
        poly.inner_rings.reserve(inner_rings.size());
        for (const GeoPolygon<TYPE>& inner : inner_rings) {
            GeoPolygon<double> new_inner = inner.buildingSimplify(scale);
            if (!new_inner.ring.empty())
                poly.inner_rings.push_back(new_inner);
        }

        return poly;
    }


    /**
     * Apply Visvalingam-Whyatt algorithm to polygon, removing points whose
     * contribution to the polygon's area is less than MIN_AREA
     *
     * See http://archive.is/Tzq2#selection-63.0-63.27 for algorithm details
     */
    void vwReduce(double min_area) {
        // Map of point indices and associated area
        std::map<size_t,double> effectiveAreas;

        // Define lambdas to help simplify code
        auto getNextIter = [&](std::map<size_t,double>::iterator iter) {
            if (++iter == effectiveAreas.end())
                return effectiveAreas.begin();
            else
                return iter;
        };
        auto getPrevIter = [&](std::map<size_t,double>::iterator iter) {
            if (iter == effectiveAreas.begin())
                return --effectiveAreas.end();
            else
                return --iter;
        };
        auto updateArea = [&](const std::map<size_t,double>::iterator& iter) {
            iter->second = area(getPrevIter(iter)->first,iter->first,getNextIter(iter)->first);
        };

        // Add indices to list
        size_t max_ix = ring.front() == ring.back() ? ring.size()-1 : ring.size();
        for (size_t ix = 0; ix < max_ix; ++ix)
            effectiveAreas[ix] = 0; // Initialize area with 0

        // Calculate areas
        for (std::map<size_t,double>::iterator iter = effectiveAreas.begin(); iter !=  effectiveAreas.end(); ++iter)
            updateArea(iter);

        // Remove points until we hit the area limit, or only have three points left
        for (size_t n = max_ix; n > 3; --n) {
            auto iter = std::min_element(effectiveAreas.begin(),effectiveAreas.end(),
                    [](const std::pair<size_t,double>& l, const std::pair<size_t,double>& r) {
                return l.second < r.second; // Sort by area
            });
            if (iter->second > min_area)
                break;
            else {
                auto prev = getPrevIter(iter);
                auto next = getNextIter(iter);
                effectiveAreas.erase(iter);
                updateArea(prev);
                updateArea(next);
            }
        }

        // Create new polygon, and copy over points
        std::vector<Point> new_ring;
        new_ring.reserve(ring.size());
        for (auto pair : effectiveAreas)
            new_ring.push_back(ring[pair.first]);
        // Close polygon
        new_ring.push_back(new_ring.front());
        // Swap
        ring.swap(new_ring);
    }
private:
    // Convert points to OGRLinearRing
    OGRLinearRing makeOGRLinearRing() const {
        OGRLinearRing ogrRing;
        for (const Point& pt : ring)
            ogrRing.addPoint(easting+xscale*pt.x,northing+yscale*pt.y);
        return ogrRing;
    }

    // Collection of functions that advance/decrement an index into ring, providing
    // for looping back to the start at the end
    size_t prev(size_t ix) const {
        return (ix+ring.size()-1) % ring.size();
    }
    size_t next(size_t ix) const {
        return (ix+1) % ring.size();
    }
    size_t incr(size_t ix, size_t delta) const {
        return (ix + delta) % ring.size();
    }
    size_t decr(size_t ix, size_t delta) const {
        delta = delta % ring.size();
        return delta > ix ? (ix + ring.size() - delta) : ix - delta;
    }

    // Calculate distance between points of indices i and j
    double dist(size_t i, size_t j) const {
        return (ring[i]-ring[j]).abs();
    }

    // Calculate angle at the jth ring point
    double angle(size_t j) const {
        return angle(prev(j),j,next(j));
    }

    // Calculate angle between ring points i, j, and k
    double angle(size_t i, size_t j, size_t k) const {
        return angle(ring[i],ring[j],ring[k]);
    }

    // Calculate angle between points A, B, and C
    double angle(const Point& A, const Point& B, const Point& C) const {
        Point dAB = B-A;
        Point dCB = B-C;
        double xp = dAB.dot(dCB)/(dAB.abs()*dCB.abs());
        if (xp > 1)
            return 0;
        else if (xp < -1)
            return M_PI;
        else
            return acos(xp);

    }

    // Calculate area of triangle defined by ring points a, b, and c
    double area(size_t a, size_t b, size_t c) const {
        return 0.5*abs((ring[a].x-ring[c].x)*(ring[b].y-ring[a].y)-(ring[a].x-ring[b].x)*(ring[c].y-ring[a].y));
    }

    /**
     * Determine points that are local corners, based on between a given point
     * and those points that are +/- scale in the boundary away from it.
     *
     * The corner metric is (pi-theta)*l_AB*l_BC, where l_XY is the distance
     * between two points and theta is the angle between A, B, and C.
     * The minimum threshold for the corner metric is equivalent to an angle of
     * 180-22.5 degrees between two otherwise straight lines; if the boundary
     * points between A, B and C are not straight, the corner metric is reduced
     * to account for that.
     */
    std::list<size_t> findCorners(unsigned int scale) const {
        std::vector<double> areas(ring.size());
        std::list<size_t> possible_peaks;
        std::vector<bool> suppressed_peaks(ring.size());

        double min_area = scale*scale / 8.0 * M_PI;
        size_t min_pk_dist = scale*5/6;

        // Calculate areas [corner metric] for each point
        for (size_t j=0; j<ring.size(); ++j) {
            size_t i = decr(j,scale);
            size_t k = incr(j,scale);
            areas[j] = (M_PI-angle(i,j,k))*dist(i,j)*dist(j,k);
            if (areas[j] >= min_area)
                possible_peaks.push_back(j);
        }

        // Sort indices by decreasing area
        possible_peaks.sort([&](const size_t& l, const size_t& r) {
            return areas[l] >= areas[r]; // Reverse sort
        });

        // Non-maximal suppression
        std::list<size_t>::iterator iter = possible_peaks.begin();
        while (iter != possible_peaks.end()) {
            size_t i = *iter;
            if (areas[i] >= areas[next(i)] && areas[i] >= areas[prev(i)] && !suppressed_peaks[i]) {
                size_t j = i;
                size_t k = i;
                for (size_t l=0; l<min_pk_dist; ++l) {
                    j = next(j);
                    k = prev(k);
                    suppressed_peaks[j] = true;
                    suppressed_peaks[k] = true;
                }
                ++iter;
            } else {
                possible_peaks.erase(iter++);
            }
        }

        // Sort indices by index
        possible_peaks.sort();

        return possible_peaks;
    }


    /**
     * Add point indices to the list to make sure all boundary (ring) points are within
     * MIN_DISTANCE of the polygon defined by POINTS
     */
    void addMissingPoints(std::list<size_t>& points, double min_distance) const {
        std::list<size_t>::iterator curr_iter = points.begin();

        // Helpful lambda
        auto getNextIter = [&](std::list<size_t>::iterator iter) {
            if (++iter == points.end())
                return points.begin();
            else
                return iter;
        };

        bool first_point = true;
        size_t start_ix = points.front();
        while (first_point || *curr_iter != start_ix) {
            std::list<size_t>::iterator next_iter = getNextIter(curr_iter);
            if (next_iter == points.end())
                next_iter = points.begin();
            size_t i = *curr_iter;
            size_t k = *next_iter;
            double l = dist(i,k);

            // For each ring point between i and k, find the maximum distance
            // between the point and the straight line between i and k
            double max_d = 0;
            size_t max_ix = 0;
            for (size_t j = next(i); j != k; j=next(j)) {
                double d = 2*area(i,j,k)/l;
                if (d > max_d) {
                    max_d = d;
                    max_ix = j;
                }
            }

            if (max_d >= min_distance) {
                // If point is too far away, add it to the list
                points.insert(next_iter,max_ix);
            } else {
                // Otherwise, move on to the next line segment
                curr_iter = next_iter;
                if (first_point)
                    first_point = false;
            }
        }
    }

    /**
     * Create a new polygon, where corner points are adjusted to fractions of a pixel based on the
     * weighted best fit of to boundary points between each corner, and corners are added until
     * all angles are greater than 60 degrees.
     */
    std::list<DblPoint> improvePoints(std::list<size_t>& points, unsigned int scale) const {
        std::vector<bool> verified(ring.size()); // List of point indices, and whether they're improved point has been checked
        std::list<size_t>::iterator curr_iter = points.begin(); // Current iterator into list of point indices

        std::map<size_t,DblPoint> pt1s; // Map of boundary pixel index to improved point location
        double eps = 1e-6; // Determines how parallel line segments must be before they are considered co-linear

        //----  Series of helpful lambdas ----//
        // Gets an incremented iterator, wrapping at end
        auto getNextIter = [&](std::list<size_t>::iterator iter) {
            if (++iter == points.end())
                return points.begin();
            else
                return iter;
        };
        // Gets a decremented iterator, wrapping back at start
        auto getPrevIter = [&](std::list<size_t>::iterator iter) {
            if (iter == points.begin())
                return --points.end();
            else
                return --iter;
        };
        // Invalidate improved point at index i
        auto voidPoint1 = [&](size_t i) {
            if (pt1s.count(i))
                pt1s.erase(i);
            verified[i]=false;
        };
        // Tries to add a boundary point between indices *I and *K; returns true if list of points has been modified
        auto insertBetween = [&](std::list<size_t>::iterator I, std::list<size_t>::iterator K) {
            // We want to insert a point; find the point between I and K that will add the most area
            double max_a = 0;
            size_t max_ix = 0;
            for (size_t j = next(*I); j!= *K; j=next(j)) {
                double a = area(*I,j,*K);
                if (a > max_a) {
                    max_a = a;
                    max_ix = j;
                }
            }
            if (max_a > eps) {
                // We're going to modify either *I, *K, or add a point in-between -- either way, the
                // calculated improved points are invalidated.
                voidPoint1(*I);
                voidPoint1(*K);

                if (M_PI - angle(*getPrevIter(I),*I,max_ix) < eps)
                    *I = max_ix; // New point is co-linear with *I and the preceding point; replace *I with new point
                else if (M_PI - angle(max_ix,*K,*getNextIter(K)) < eps)
                    *K = max_ix; // New point is co-linear with *K and the following point; replace *K with new point
                else
                    points.insert(K,max_ix);
                return true;
            } else {
                return false;
            }
        };


        // Loop through corner indices
        while (!verified[*curr_iter]) {
            std::list<size_t>::iterator prev_iter = getPrevIter(curr_iter);
            std::list<size_t>::iterator next_iter = getNextIter(curr_iter);
            // Calculate improved points (if not already calculated)
            if (!pt1s.count(*prev_iter))
                pt1s[*prev_iter] = calcNewPoint(*getPrevIter(prev_iter),*prev_iter,*curr_iter,scale);
            if (!pt1s.count(*curr_iter))
                pt1s[*curr_iter] = calcNewPoint(*prev_iter,*curr_iter,*next_iter,scale);
            if (!pt1s.count(*next_iter))
                pt1s[*next_iter] = calcNewPoint(*curr_iter,*next_iter,*getNextIter(next_iter),scale);

            bool pts_changed = false;

            // Test if angle between improved points is less than 60 degrees
            if (angle(pt1s[*prev_iter],pt1s[*curr_iter],pt1s[*next_iter]) < M_PI/3) {
                // We want to insert a point; find the point between prev and curr that will add the most area
                if ( (ring.size() + *curr_iter - *prev_iter) % ring.size() > 1)
                    pts_changed |= insertBetween(prev_iter,curr_iter);
                // Find the point between curr and next that will add the most area
                if ( (ring.size() + *next_iter - *curr_iter) % ring.size() > 1)
                    pts_changed |= insertBetween(curr_iter,next_iter);
            }

            if (!pts_changed) {
                // We didn't have to change points; mark the improved point as ok and move on
                verified[*curr_iter] = true;
                curr_iter = getNextIter(curr_iter);
            } else if (!verified[*prev_iter]) {
                // Depending on what changed, we might have to back up a step
                curr_iter = prev_iter;
            }
        }

        // Copy points from map to list & return
        std::list<DblPoint> newPoints;
        for (const size_t& i : points) {
            newPoints.push_back(pt1s[i]);
        }
        return newPoints;
    }

    /**
     * Calculate a new corner point based on the intersections between the weighted best fit lines
     * of the boundary points between i and j, and j and k.
     *
     * For lines that are essentially parallel, the returned point is the average of the closest point
     * on each line to the jth boundary point.
     */
    DblPoint calcNewPoint(size_t i, size_t j, size_t k, unsigned int scale) const {
        std::pair<DblPoint,DblPoint> AB = calcBestFitLine(i,j,scale);
        std::pair<DblPoint,DblPoint> BC = calcBestFitLine(j,k,scale);

        if (AB.first.dot(BC.first) > 1 - 1e-6) {
            // Lines are basically parallel
            return 0.5*(closestLinePointToPoint(AB,j)+closestLinePointToPoint(BC,j));
        } else {
            // Intersect lines
            double t = (BC.second-AB.second).cross(BC.first) / AB.first.cross(BC.first);
            return t*AB.first+AB.second;
        }
    }

    /**
     * Calculate the weighted best fit line for boundary [ring] points between i and k.
     * Points are weighted such that points within scale of the corner have a linearly
     * increasing weight; this is designed to help square up corners.
     *
     * Returned best fit line consists of a pair of DBlPoints; the first represents the
     * line unit vector and the second is a point the line intersects.  This description
     * of the line allows for a more parametric description, avoiding issues where the
     * slope is infinite.
     */
    std::pair<DblPoint,DblPoint> calcBestFitLine(size_t i, size_t k, unsigned int scale) const {
        std::pair<DblPoint,DblPoint> rv;

        // Assign weights
        std::vector<double> wts(ring.size(),0);
        unsigned int w = 0;
        double c = 1.0/scale;
        for (size_t j = i; j!=next(k); j=next(j))
            wts[j]=std::min(++w,scale)*c;
        w = 0;
        for (size_t j = k; j!=prev(i); j=prev(j))
            wts[j]=wts[j]*std::min(++w,scale)*c;

        // Compute least squares of parametric x and y
        for (size_t d = 0; d < 2; ++d) {
            double t = 0;
            double SW = 0;
            double SWt = 0;
            double SWz = 0;
            double SWtt = 0;
            double SWtz = 0;

            for (size_t j = i; j!=next(k); j=next(j)) {
                SW  += wts[j];
                SWt += wts[j]*t;
                SWz += wts[j]*ring[j][d];
                SWtt+= wts[j]*t*t;
                SWtz+= wts[j]*t*ring[j][d];
                ++t;
            }

            rv.first[d] = (SW*SWtz - SWt*SWz) / (SWtt*SW - SWt*SWt);
            rv.second[d] = (SWz-rv.first[d]*SWt) / SW;
        }
        return rv;
    }

    // Determine the closest point on line to the ith ring point
    DblPoint closestLinePointToPoint(const std::pair<DblPoint,DblPoint>& line, size_t i) const {
        return line.second + (ring[i]-line.second).dot(line.first)*line.first;
    }
};

} // pubgeo

#endif /* PUBGEO_GEO_POLYGON_H_ */
