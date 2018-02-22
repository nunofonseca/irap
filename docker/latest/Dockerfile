FROM fedora:27

LABEL iRAP.version="pre_release" maintainer="nuno.fonseca at gmail.com"

# docker build -f Dockerfile -t irap/latest:v0 ..
COPY build/irap_docker_setup.sh build 
RUN bash build fedora_27 master full -R


#ENTRYPOINT ["irap"]

