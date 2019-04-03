# Copyright (c) 2017-2019, The Johns Hopkins University /
# Applied Physics Laboratory (JHU/APL)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Please reference the following when reporting any results using this software:
#
# M. Bosch, A. Leichtman, D. Chilcott, H. Goldberg, M. Brown, “Metric
# Evaluation Pipeline for 3D Modeling of Urban Scenes,” ISPRS Archives, 2017.
#
# S. Almes, S. Hagstrom, D. Chilcott, H. Goldberg, M. Brown, “Open Source
# Geospatial Tools to Enable Large Scale 3D Scene Modeling,” FOSS4G, 2017.
#
# For more information, please see: http://www.jhuapl.edu/pubgeo.html

FROM pdal/pdal:1.8 as builder
MAINTAINER JHUAPL <pubgeo@jhuapl.edu>

RUN apt update && apt upgrade -y && apt install -y --fix-missing --no-install-recommends\
    build-essential \
    ca-certificates \
    cmake \
    libgdal-dev \
&& rm -rf /var/lib/apt/lists/*

# Make a directory to work out of, and change to it
RUN cd / && mkdir /pubgeo
COPY src /pubgeo/src
COPY CMakeLists.txt /pubgeo/
COPY CMakeSettings.json /pubgeo/
WORKDIR /pubgeo/build/
RUN cmake .. && make -j 10 && make install && mv shr3d align3d /usr/local/bin


FROM pdal/pdal:1.8 as runner
MAINTAINER JHUAPL <pubgeo@jhuapl.edu>

WORKDIR /
COPY --from=builder /usr/local/bin /usr/local/bin
COPY --from=builder /usr/lib/libpdal_plugin_filter_align3d* /usr/lib/
COPY --from=builder /usr/lib/libpdal_plugin_writer_shr3d* /usr/lib/

CMD echo "Please run a valid executable:" && \
    echo "docker run --rm -v <path to 3D data>:<mount point (MP)> jhuapl/pubgeo shr3d <MP>/<3D file> DH=2 DZ=1 AGL=2 AREA=50" && \
    echo "docker run --rm -v <path to point clouds (PC)>:<mount point (MP)> jhuapl/pubgeo align3d <MP>/<reference pc> <MP>/<target pc> gsd=1.0 maxt=10.0"

# Reminders:
#   Data coming out will be in <path to 3D data> or <path to point clouds (pc)>
#   --rm will auto delete the container when complete (no need to take up disk space)
#   The files will be owned by root when they come out (known Docker behaviour)
