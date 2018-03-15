FROM fedora:27

ENV IRAP_VERSION=1.0.0a
LABEL iRAP.version="$IRAP_VERSION" maintainer="nuno.fonseca at gmail.com"

# docker build -f Dockerfile -t irap/latest:v0 ..
COPY build/irap_docker_setup.sh build 
RUN bash build fedora_27 master full 


#ENTRYPOINT ["irap"]

