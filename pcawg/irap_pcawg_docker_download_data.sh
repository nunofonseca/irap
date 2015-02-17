#!/bin/bash
foldername=$1
if [ "$foldername-" == "-" ]; then
   echo "Usage: irap_pcawg_setup foldername" > /dev/stderr
   echo "This script will download the genome and annotation" > /dev/stderr
   exit 1
fi
if [ -e $foldername ]; then
   echo "Folder $foldername already exists...aborting" > /dev/stderr
   exit 1
fi

if [ "$IRAP_DOCKER_IMAGE-" == "-" ]; then
   echo "ERROR: Please define the environment variable IRAP_DOCKER_IMAGE with the ID of docker image to use"
   exit 1
fi

# download data
docker run -v `pwd`:/irap_data --entrypoint="env" -i -t $IRAP_DOCKER_IMAGE bash -c "cd /irap_data && /opt/irap_clone/pcawg/pcawg_download_data.sh /irap_data/$foldername"
