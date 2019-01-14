FROM ubuntu:18.04 as base
RUN apt update && apt upgrade -y && apt install -y --fix-missing --no-install-recommends\
    ca-certificates \
    curl \
    gdal-bin \
    pdal \
&& rm -rf /var/lib/apt/lists/*

FROM base as build
RUN apt update && apt upgrade -y && apt install -y --fix-missing --no-install-recommends\
    cmake \
    libgdal-dev \
    libpdal-dev \
    build-essential \
    git \
&& rm -rf /var/lib/apt/lists/*

RUN cd / && git clone https://github.com/pubgeo/pubgeo

# Make a directory to work out of, and change to it
WORKDIR /pubgeo/build/
RUN cmake .. && make -j 10 && make install && mv shr3d align3d /usr/local/bin
WORKDIR /

FROM base as release
MAINTAINER JHUAPL <pubgeo@jhuapl.edu>

WORKDIR /
COPY --from=build /usr/local /usr/local

CMD echo "Please run a valid executable:" && \
    echo "docker run --rm -v <path to 3D data>:<mount point (MP)> jhuapl/pubgeo shr3d <MP>/<3D file> DH=2 DZ=1 AGL=2 AREA=50" && \
    echo "docker run --rm -v <path to point clouds (PC)>:<mount point (MP)> jhuapl/pubgeo align3d <MP>/<reference pc> <MP>/<target pc> gsd=1.0 maxt=10.0"


# Reminders:
#   Data coming out will be in <path to 3D data> or <path to point clouds (pc)>
#   --rm will auto delete the container when complete (no need to take up disk space)
#   The files will be owned by root when they come out (known Docker behaviour)
