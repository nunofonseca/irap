
IRAP_VERSION=08122016

#texlive-texinfo
#cp /sbin/install-info /usr/bin/install-info
#evince automake	ksh r_cv_prog_makeinfo_v4

gen_images:
	IMAGE1=`cat fedora20_full.docker | sudo docker build --tag fedora20/irapf:v$(IRAP_VERSION) -  | grep "Successfully built" | cut -f 3 -d\ ` && \
	IMAGE2=`cat fedora20_light.docker | sudo docker build --tag fedora20/irapl:v$(IRAP_VERSION) -  | grep "Successfully built" | cut -f 3 -d\ ` && \
	IMAGE3=`cat fedora22_light.docker | sudo docker build --tag fedora22/irapl:v$(IRAP_VERSION) -  | grep "Successfully built" | cut -f 3 -d\ ` && \
	IMAGE4=`cat ubuntu_latest.docker | sudo docker build --tag ubuntu/irapl:v$(IRAP_VERSION) -  | grep "Successfully built" | cut -f 3 -d\ ` && \
	IMAGE5=`cat ubuntu14_light.docker | sudo docker build --tag ubuntu14/irapl:v$(IRAP_VERSION) -  | grep "Successfully built" | cut -f 3 -d\ ` && \
	IMAGE6=`cat fedora22_full.docker | sudo docker build --tag fedora22/irapf:v$(IRAP_VERSION) -  | grep "Successfully built" | cut -f 3 -d\ ` && \
	echo "All done!!!!"
	echo "Images: $(IMAGE1) $(IMAGE2) $(IMAGE3) $(IMAGE4) $(IMAGE5)"
