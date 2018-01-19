FROM ubuntu:17.10
MAINTAINER JHUAPL <pubgeo@jhuapl.edu>

RUN apt update && apt upgrade -y && apt install -y --fix-missing --no-install-recommends\
    build-essential \
    ca-certificates \
	cmake \
	curl \
	gdal-bin \
	git \
	libgdal-dev \
	libpdal-dev \
	pdal \
&& rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/pubgeo/pubgeo

WORKDIR /build/
# Make a directory to work out of, and change to it
RUN cmake ../pubgeo && make -j 10

# cleanup
RUN rm -rf /pubgeo
RUN apt purge -y \
    build-essential \
    libgdal-dev \
    libpdal-dev \
    cmake \
    git
CMD echo "Please run a valid executable:" && \
    echo "docker run --rm -v <path to 3D data>:<mount point (mp)> jhuapl/pubgeo ./shr3d <mount point>/<3D file> DH=2 DZ=1 AGL=2 AREA=50" && \
    echo "docker run --rm -v <path to point clouds (pc)>:<mount point (mp)> jhuapl/pubgeo ./align3d <mp>/<reference pc> <mp>/<target pc> gsd=1.0 maxt=10.0"

# Examples
# docker run -v /home/ubuntu/pointclouds:pc jhuapl/pubgeo ./shr3d /pc/shr3dMe.las DH=2 DZ=1 AGL=2 AREA=50
# docker run -v /home/ubuntu/pointclouds:pc jhuapl/pubgeo ./align3d /pc/reference.las pc/alignMe.las gsd=1.0 maxt=10.0
# Reminders:
#   -t makes output go to terminal instantly (instead of bulk dump at end)
#   Data coming out will be in <path to 3D data> or <path to point cloud>
#   --rm will auto delete the container when complete (no need to take up disk space)
#   The files will be owned by root when they come out (known Docker behaviour)
