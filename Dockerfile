FROM ubuntu:17.04
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

# Make a directory to work out of, and change to it
WORKDIR /build
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
    echo "docker run -v <path to 3D data> shr3d <3D file> DH=2 DZ=1 AGL=2 AREA=50" && \
    echo "docker run -v <path to point cloud> align3d <reference point cloud> <target pc> gsd=1.0 maxt=10.0"
