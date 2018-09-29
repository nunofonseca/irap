FROM amazonlinux:latest

ENV IRAP_VERSION=devel

LABEL iRAP.version="$IRAP_VERSION" maintainer="nuno.fonseca at gmail.com"
ADD build/irap_docker_setup.sh build
RUN bash build aws $IRAP_VERSION light 
