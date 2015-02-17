#!/bin/bash
sop=$1
foldername=$2
if [ "$foldername-" == "-" ]; then
   echo "Usage: irap_pcawg_docker_setup SOP foldername" > /dev/stderr
   echo "This script will generate the indexes" > /dev/stderr
   exit 1
fi
if [ ! -e $foldername ]; then
   echo "Folder $foldername not found" > /dev/stderr
   exit 1
fi

if [ "$sop-" == "-" ]; then
    echo "$sop not defined"
    exit 1
fi

if [ "$IRAP_DOCKER_IMAGE-" == "-" ]; then
   echo "ERROR: Please define the environment variable IRAP_DOCKER_IMAGE with the ID of docker image to use"
   exit 1
fi

# download data
docker run -v `pwd`:/irap_data --entrypoint="env" -i -t $IRAP_DOCKER_IMAGE bash -c "cd /irap_data && irap conf=test_data/pcawg.conf sop=$sop data_dir=$foldername"
