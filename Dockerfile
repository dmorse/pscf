# This is the Dockerfile for PSCF.
#
# To build the container, you need docker to be installed on
# your machine (https://www.docker.com) and run the following command
#
#   docker build -t pscf .
#
# It will create a Debian-based container that contains everything
# needed to compile pscf
#
# Everything will be installed in /pscf
#
# To run the container interactively
#
#   docker run -it pcsf /bin/bash

FROM debian

RUN apt-get update && \
    apt-get install -y \
        cmake \
        gfortran \
        liblapack3 \
        libfftw3-dev \
        python

WORKDIR .
ADD . /pcsf
RUN cd /pcsf_SRC && \
    mkdir -p build && \
    cd build && \
    export FC=`which gfortran` && \
    cmake .. -DCMAKE_BUILD_TYPE=Release \
             -DCMAKE_INSTALL_PREFIX=/pcsf && \
    make && \
    make install
