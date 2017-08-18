# JHU/APL pubgeo
JHU/APL is committed to advancing the state of the art in geospatial computer vision by advocating for public data sets 
and open source software. For more information on this and other efforts, please visit [JHU/APL](http://www.jhuapl.edu/pubgeo.html).

## SHR3D

Shareable High Resolution 3D Point Cloud Classification
>SHR3D classifies a 3D point cloud to generate simple raster products, including a Digital Surface Model (DSM), Digital 
Terrain Model (DTM), and a classification image labeling ground, buildings, and trees.

### SHR3D Usage
    ./shr3d <Input File(LAS|TIF)> DH=<Horizontal uncertainty(m)> DZ=<Vertical uncertainty(m)> AGL=<Min building height(m)> AREA=<Min building area(m^2)>

#### Example:
    ./shr3d shr3dMe.las DH=0.5 DZ=0.5 AGL=2 AREA=50
This will produce a multiple tif files as a result:
* **Classification image**- 'Colored' by classification for each pixel: buildings, ground, and trees
* **Building image**- Binary mask of all detected buildings
* **DSM**- Digital surface model (reflective surface)
* **DTM**- Digital terrain model (bare earth, or ground)

## ALIGN3D
Align 3D Point Cloud Registration Tool
>ALIGN3D estimates and applies a transform to align 3D point clouds. This algorithm was developed 
for use with airborne lidar, multiple view satellite imagery, and synthetic aperture radar 
derived point clouds.

### ALIGN3D Usage
    ./align3d <Reference point cloud(LAS)> <Target point cloud(LAS)> maxdz=<Maximum local z difference(m)> gsd=<Ground sample distance(m)> maxt=<Maximum XY translation for search(m)>

#### Example:
    ./align-3d reference.las target.las maxt=10.0 gsd=1.0
This will produce multiple output files:
* Aligned.las- A point cloud representation of the target, aligned to the reference
* Aligned.tif- A DSM of the target, aligned to the reference
* Offsets.txt- A list of the xyz offsets to translate the target to align with the reference 
