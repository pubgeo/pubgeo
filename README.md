# JHU/APL pubgeo
Open source software and public data for geospatial computer vision.

## SHR3D

Shareable High Resolution 3D Point Cloud Classification
>SHR3D classifies a 3D point cloud to generate simple geospatial products, including: a Digital Surface Model (DSM); a Digital 
Terrain Model (DTM); a classification image labeling ground, buildings, and trees; and a simplified building outlines vector product.

### SHR3D Usage
    ./shr3d <Input File (LAS|TIF)> DH=<Horizontal uncertainty(m)> DZ=<Vertical uncertainty(m)> AGL=<Min building height(m)> <Options>
    Options:
        AREA=<Min building area(m^2)>   Buildings smaller than this size will not be labeled (default: 50)
        EGM96                           Set this flag to write vertical datum
        BOUNDS=<MINX,MAXX,MINY,MAXY>    Set to define image bounds
        DTM=<DTM File (TIF)>            Path to optional DTM file
        GND_LABEL=<LABEL>               Set to the value of the point cloud ground classification label (typically 2) if the point cloud has already been partially classified.

#### Example:
    ./shr3d shr3dMe.las DH=0.5 DZ=0.5 AGL=2 AREA=50
This will produce multiple files as a result:
* Raster files (*.tif):
  * **Classification image**- 'Colored' by classification for each pixel: buildings, ground, and trees
  * **Building image**- Binary mask of all detected buildings
  * **DSM**- Digital surface model (reflective surface)
  * **DTM**- Digital terrain model (bare earth, or ground)
  * **INT**- Intensity image
* Vector files (shapefiles: *.shp, *.prj, *.shx, *.dbf):
  * **Building outlines**- Simplified polygons of detected building perimeters

## ALIGN3D
Align 3D Point Cloud Registration Tool
>ALIGN3D estimates and applies a transform to align 3D point clouds. This algorithm was developed 
for use with airborne lidar, multiple view satellite imagery, and synthetic aperture radar 
derived point clouds.

### ALIGN3D Usage
    ./align3d <Reference point cloud(LAS)> <Target point cloud(LAS)> maxdz=<Maximum local z difference(m)> gsd=<Ground sample distance(m)> maxt=<Maximum XY translation for search(m)>

#### Example:
    ./align3d reference.las target.las maxt=10.0 gsd=1.0
This will produce multiple output files:
* Aligned.las- A point cloud representation of the target, aligned to the reference
* Aligned.tif- A DSM of the target, aligned to the reference
* Offsets.txt- A list of the xyz offsets to translate the target to align with the reference 
