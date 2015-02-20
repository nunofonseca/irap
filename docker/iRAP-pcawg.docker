FROM fedora:20
MAINTAINER Nuno Fonseca email: nuno.fonseca at gmail.com

# Update the image with the latest packages (recommended)
# and install missing packages
RUN sed -i "s/failovermethod=priority/failovermethod=roundrobin/" /etc/yum.repos.d/fedora.repo  && yum update -y && yum install -y zlib-devel python-devel bzip2-devel python readline-devel libgfortran gcc-gfortran gcc-c++ libX11-devel libXt-devel numpy gd-devel libxml2-devel libxml2 libpng libcurl-devel expat-devel  libpangocairo db-devel java python gcc g++ gfortran curl-config git which make bzip2 bison gettext-devel  unzip make wget sqlite sqlite-devel db4-devel libdb-devel graphviz tar && yum clean all


# 
WORKDIR /opt
RUN git clone http://code.google.com/p/irap/ irap_clone  
WORKDIR /opt/irap_clone
RUN git checkout devel && ./pcawg/irap_install_pcawg.sh -a /irap_install -s . && rm -rf /irap_install/tmp 
RUN echo source /irap_install/irap_setup.sh >> ~/.bash_profile && echo source /irap_install/irap_setup.sh >> ~/.bashrc

RUN echo '#!/bin/env bash' > /usr/bin/irap
RUN echo 'source /irap_install/irap_setup.sh' >> /usr/bin/irap
RUN echo '/irap_install/scripts/irap "$@"' >> /usr/bin/irap
RUN chmod u+x /usr/bin/irap
#ENTRYPOINT ["irap"]

