FROM ubuntu:18.04
# Initially created by Natalja Kurbatova
# Update the image with the latest packages (recommended)
# and install missing packages

ENV DEBIAN_FRONTEND=noninteractive
ENV IRAP_VERSION=1.0.5

LABEL iRAP.version="$IRAP_VERSION" maintainer="nuno.fonseca at gmail.com"


ADD build/irap_docker_setup.sh build
RUN bash build ubuntu $IRAP_VERSION full 

#ENTRYPOINT ["irap"]
