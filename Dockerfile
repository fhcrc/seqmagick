# Stage 1: Build and test the package
FROM docker.io/python:3.6-slim-stretch AS build

ADD . /seqmagick
WORKDIR /seqmagick

RUN pip3 install numpy
RUN pip3 install biopython
RUN python3 setup.py install
RUN python3 setup.py test

# Stage 2: Only install the stuff needed for runtime
FROM docker.io/python:3.6-alpine

COPY --from=build /usr/local/bin/seqmagick /usr/local/bin/
COPY --from=build /usr/local/lib/python3.6/site-packages /usr/local/lib/python3.6/site-packages
