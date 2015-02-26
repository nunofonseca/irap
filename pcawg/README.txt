iRAP docker wrapper for PanCancer (http://pancancer.info/)

Setup

1. Create a docker image. There are two docker files in the docker
   folder related to PanCancer (pcawg): 
     a) one docker file will install the devel branch of iRAP and
     b) the other the latest release.

   Examples:
    a)
     cd docker; cat iRAP-devel-pcawg.docker| sudo docker build --rm --tag fedora/irap-pcawg:devel -
    b)
     cd docker; cat iRAP-pcawg.docker| sudo docker build --rm --tag fedora/irap:devel -

2.  The scripts will look for the docker image to use (the image ID
    can be found by running the command 'docker images') in the
    environment variable IRAP_DOCKER_IMAGE.

3. Download the reference, annotation, and other data necessary for
   the analysis. The script pcawg/irap_pcawg_docker_download_data.sh
   will install all data in the user defined folder. 
   Usage:    <path2irap_clone>/pcawg/irap_pcawg_docker_download_data.sh <folderName>

   Example:
   <path2irap_clone>/pcawg/irap_pcawg_docker_download_data.sh test_data

4. Generate the indices and other files needed by iRAP
   Usage: <path2irap_clone>/pcawg/irap_pcawg_docker_setup.sh sop=<star or tophat2 sop's> data_dir=<folderName>
   Example: <path2irap_clone>/pcawg/irap_pcawg_docker_setup.sh pawg3_th2_mapping test_data
            <path2irap_clone>/pcawg/irap_pcawg_docker_setup.sh pawg3_star_mapping test_data

Align

The irap_pcawg_align script can be used to run iRAP to align the fastq
files contained in a tarball and obtain a BAM file. The produced BAM
will be named as <ID>.<Mapper>.v1.bam and placed on the current
folder.

Usage: <path2irap_clone>/pcawg/irap_pcawg_align <star or tophat2 SOP> <tarball> <sampleId> [max_threads=N max_memory=M]

Example: <path2irap_clone>/pcawg/irap_pcawg_align pawg3_th2_mapping <path2irap_clone>/pcawg/examples/test-SE-1Ls_rnaseq_fastq.tar 6664e287-c7a6-417d-85d5-4893dcb1ec2b "max_threads=2"

