// Copyright (c) 2017-2019, The Johns Hopkins University /
// Applied Physics Laboratory (JHU/APL)
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

// Please reference the following when reporting any results using this software:
//
// M. Bosch, A. Leichtman, D. Chilcott, H. Goldberg, M. Brown, “Metric
// Evaluation Pipeline for 3D Modeling of Urban Scenes,” ISPRS Archives, 2017.
//
// S. Almes, S. Hagstrom, D. Chilcott, H. Goldberg, M. Brown, “Open Source
// Geospatial Tools to Enable Large Scale 3D Scene Modeling,” FOSS4G, 2017.
//
// For more information, please see: http://www.jhuapl.edu/pubgeo.html

#ifndef PUBGEO_SHR3D_H
#define PUBGEO_SHR3D_H

#include <cstdio>
#include <cstring>
#include <map>
#include <memory>
#include <string>

#include <pdal/util/Bounds.hpp>

#include "geo_polygon.h"
#include "orthoimage.h"

#define LABEL_OBJECT    1
#define LABEL_GROUND    UINT_MAX
#define LABEL_IN_ONE    UINT_MAX-1
#define LABEL_ACCEPTED  UINT_MAX-2
#define LABEL_TEMP      UINT_MAX-3
#define LABEL_BOUNDARY  UINT_MAX-5

/*
 * Las Classifications list
    0 Created, never classified
    1 Unclassified1
    2 Ground
    3 Low Vegetation
    4 Medium Vegetation
    5 High Vegetation
    6 Building
    7 Low Point (noise)
    8 Model Key-point (mass point)
    9 Water
    10 Reserved for ASPRS Definition
    11 Reserved for ASPRS Definition
    12 Overlap Points2
    13-31 Reserved for ASPRS Definition
 */
#define LAS_UNCLASSIFIED    (1*40)
#define LAS_GROUND          (2*40)
#define LAS_TREE            (5*40)
#define LAS_BUILDING        (6*40)
// For now, multiply these by 40 so we can view the images easily.

namespace shr3d {
    using namespace pubgeo;

    typedef struct {
        unsigned int i;
        unsigned int j;
    } PixelType;

    typedef struct {
        unsigned int label;
        long int xmin;
        long int xmax;
        long int ymin;
        long int ymax;
        long int count;
    } ObjectType;

    enum ImageType {
        DSM,
        MIN,
        DTM,
        INTENSITY,
        CLASS,
        BUILDING,
        BUILDING_OUTLINES,
        DSM2,
        MINAGL,
        DTM0,
        LABEL0,
        LABEL,
        LABELED_BUILDINGS,
        LABELED_BUILDINGS_3,
    };

    class Shr3dder {
    public:
        // Options
        double dh_meters;
        double dz_meters;
        double agl_meters;
        double min_area_meters;
        double max_tree_height_meters;
        bool egm96;
        pdal::Bounds bounds;
        int gnd_label;

        // Constructor
        Shr3dder() : dh_meters(1), dz_meters(1), agl_meters(2),
                min_area_meters(50), max_tree_height_meters(40),
                egm96(false), gnd_label(-1), bldgProcessed(false) {}

        // Function declarations.
        void process(std::map<ImageType,std::string> outputFilenames);

        bool setDSM(std::string dsmFile)    { return dsmImage.read(dsmFile.c_str()); }
        bool setDTM0(std::string dtmFile);
        bool setPSET(std::string psetFile)  { return pset.Read(psetFile.c_str()); }
        bool setPSET(const pdal::PointViewPtr view)   { return pset.Read(view); }
        
        bool is_pset_empty() { return !pset.numPoints; }
        
        void reset() {
            dsmImage.Deallocate();
            minImage.Deallocate();
            dtmImage.Deallocate();
            intImage.Deallocate();
            classImage.Deallocate();
            bldgImage.Deallocate();
            dsm2Image.Deallocate();
            minAglImage.Deallocate();
            labelImage.Deallocate();
            bldgLabels.Deallocate();
            bldgLabels3.Deallocate();
            bldgOutlines.clear();
            bldgProcessed = false;
        }

#define IMGET(im,func) if(im.empty()) func(); return im
        const OrthoImage<unsigned short> &getDSM()  {IMGET(dsmImage,createDSM);}
        const OrthoImage<unsigned short> &getMIN()  {IMGET(minImage,createMIN);}
        const OrthoImage<unsigned short> &getDSM2() {IMGET(dsm2Image,createDSM2);}
        const OrthoImage<unsigned short> &getDTM0() {IMGET(dtm0Image,createDTM0);}
        const OrthoImage<unsigned int>   &getLBL0() {IMGET(label0Image,createDTM0);}
        const OrthoImage<unsigned short> &getDTM()  {IMGET(dtmImage,createDTM);}
        const OrthoImage<unsigned int>   &getLBL()  {IMGET(labelImage,createDTM);}
        const OrthoImage<unsigned short> &getINT()  {IMGET(intImage,createIntensity);}
        const OrthoImage<unsigned short> &getMINAGL()   {IMGET(minAglImage,createMinAGL);}
        const OrthoImage<unsigned char> &getCLS()   {IMGET(classImage,labelClasses);}
        const OrthoImage<unsigned char> &getBLDG()  {IMGET(bldgImage,labelBuildings);}
        const OrthoImage<int> &getBLDGLBL() {IMGET(bldgLabels,createOutlines);}
        const OrthoImage<int> &getBLDGLBL3(){IMGET(bldgLabels3,createOutlines);}
        const std::map<int,GeoPolygon<double>>& getBLDGPOLY() {
            if (!bldgProcessed)
                createOutlines();
            return bldgOutlines;
        }
#undef IMGET
        
    protected:
        PointCloud pset;
        OrthoImage<unsigned short> dsmImage;
        OrthoImage<unsigned short> minImage;
        OrthoImage<unsigned short> dtm0Image;
        OrthoImage<unsigned short> dtmImage;
        OrthoImage<unsigned short> intImage;
        OrthoImage<unsigned char> classImage;
        OrthoImage<unsigned char> bldgImage;
        OrthoImage<unsigned short> dsm2Image;
        OrthoImage<unsigned short> minAglImage;
        OrthoImage<unsigned int> label0Image;
        OrthoImage<unsigned int> labelImage;
        OrthoImage<int> bldgLabels;
        OrthoImage<int> bldgLabels3;
        std::map<int,GeoPolygon<double>> bldgOutlines;
        bool bldgProcessed;
        
        void createDSM();

        void createMIN();

        void createDSM2();

        void createDTM0();

        void createDTM();

        void createMinAGL();

        void createIntensity();

        void labelClasses();

        void labelBuildings();
        
        void createOutlines();

        void classifyGround(OrthoImage<unsigned int> &labelImage, const OrthoImage<unsigned short> &dsmImage,
                                   OrthoImage<unsigned short> &dtmImage, int dhBins, unsigned int dzShort);

        void classifyNonGround(const OrthoImage<unsigned short> &dsmImage, const OrthoImage<unsigned short> &dtmImage,
                                      OrthoImage<unsigned int> &labelImage, unsigned int dzShort,
                                      unsigned int aglShort,
                                      float minAreaMeters);

        void fillInsideBuildings(OrthoImage<unsigned char> &classImage);
    };
}

#endif // PUBGEO_SHR3D_H
