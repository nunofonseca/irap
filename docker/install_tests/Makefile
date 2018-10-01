
IRAP_VERSION=31092018

#texlive-texinfo
#cp /sbin/install-info /usr/bin/install-info
#evince automake	ksh r_cv_prog_makeinfo_v4

gen_images:
	IMAGE5=`docker build --tag fedora28/irapl:v$(IRAP_VERSION) -f  fedora28_light.docker ..  | tee - /dev/stderr  | grep "Successfully built" | cut -f 3 -d\ `
	IMAGE3=`docker build --tag ubuntu/irapl:v$(IRAP_VERSION) -f ubuntu_latest_light.docker  .. | tee - /dev/stderr  | grep "Successfully built" | cut -f 3 -d\ `
	IMAGE1=`docker build --tag fedora27/irapf:v$(IRAP_VERSION) -f fedora27_full.docker  .. | tee - /dev/stderr  | grep "Successfully built" | cut -f 3 -d\ `
	IMAGE2=`docker build --tag fedora27/irapl:v$(IRAP_VERSION) -f  fedora27_light.docker ..  | tee - /dev/stderr  | grep "Successfully built" | cut -f 3 -d\ `
	IMAGE4=`docker build --tag ubuntu/irapf:v$(IRAP_VERSION) -f  ubuntu_latest_full.docker  ..  | tee - /dev/stderr | grep "Successfully built" | cut -f 3 -d\ ` && \
	echo "All done!!!!"
	echo "Images: $(IMAGE1) $(IMAGE2) $(IMAGE3) $(IMAGE4)"
