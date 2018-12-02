FROM ubuntu:16.04
# Update the image with the latest packages (recommended)
# and install missing packages

ENV DEBIAN_FRONTEND=noninteractive
ENV IRAP_VERSION=1.0.2d1

LABEL iRAP.version="$IRAP_VERSION" maintainer="nuno.fonseca at gmail.com"

ADD build/irap_docker_setup.sh build
RUN bash build ubuntu_16 $IRAP_VERSION full 

#ENTRYPOINT ["irap"]
