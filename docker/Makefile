
IRAP_VERSION=0.4.3d

# Create a docker container
iRAP.docker.tar.gz: ../Makefile iRAP.docker
	IMAGE=`cat iRAP.docker| sudo docker build --rm --tag fedora/irap:v$(IRAP_VERSION) -  | grep "Successfully built" | cut -f 3 -d\ ` && \
	docker save $$IMAGE |gzip> $@.tmp && mv $@.tmp $@

iRAP-pcawg.docker.tar.gz: ../Makefile iRAP-pcawg.docker
	IMAGE=`cat iRAP-pcawg.docker| sudo docker build --rm --tag fedora/irap_pcawg:v$(IRAP_VERSION) -  | grep "Successfully built" | cut -f 3 -d\ ` && \
        docker save $$IMAGE |gzip> $@.tmp && mv $@.tmp $@

# run irap

