FROM fedora:22
MAINTAINER Nuno Fonseca email: nuno.fonseca at gmail.com

# Update the image with the latest packages (recommended)
# and install missing packages
RUN dnf update -y && dnf install -y boost-devel R zlib-devel python-devel bzip2-devel python readline-devel libgfortran gcc-gfortran gcc-c++ libX11-devel libXt-devel numpy gd-devel libxml2-devel libxml2 libpng texi2html libcurl-devel expat-devel  pango-devel cairo-devel  java python gcc gcc-c++ gcc-objc++  gcc-gfortran curl git which make bzip2 bison gettext-devel  unzip make wget sqlite sqlite-devel db4-devel libdb-devel graphviz texlive tar java-devel texinfo  xorg-x11-server-Xvfb texlive-incons*  automake evince  texlive-collection-latex firefox tk-devel tcl-devel libtool && dnf clean all

WORKDIR /opt

# Install and clean in a single layer
# Full iRAP installation (subset of supported mappers and quantification tools)
# GCC 4.8 will be installed to compile mono (required by osa)
RUN git clone http://github.com/nunofonseca/irap.git  irap_clone  && cd /opt/irap_clone && git checkout irap_new_release && ./scripts/irap_install.sh -a /irap_install -s . -G && cd / && rm -rf /irap_install/tmp /opt/irap_clone  /root/.cpan*

WORKDIR /
RUN echo source /irap_install/irap_setup.sh >> ~/.bash_profile && echo source /irap_install/irap_setup.sh >> ~/.bashrc

RUN echo '#!/usr/bin/env bash' > /usr/bin/irap
RUN echo 'source /irap_install/irap_setup.sh' >> /usr/bin/irap
RUN echo '/irap_install/scripts/irap "$@"' >> /usr/bin/irap
RUN chmod u+x /usr/bin/irap

#ENTRYPOINT ["irap"]

