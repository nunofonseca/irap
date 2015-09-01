FROM ubuntu:14.04
MAINTAINER Nuno Fonseca email: nuno.fonseca at gmail.com
# Initially created by Natalja Kurbatova
# Update the image with the latest packages (recommended)
# and install missing packages
# docker build --no-cache -t iRAP:latest
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \ 
	build-essential \     
	texlive \
	xvfb \
	zlibc \
	zlib1g \
	zlib1g-dev \
	libncurses5-dev \
	sqlite \ 
	sqlite3 \
        libsqlite3-dev \
	gettext \
	python \
	python-dev  \
	gfortran \
	bzip2 \ 
        libbz2-1.0 \
	libbz2-dev \ 
	libreadline-dev \
	libx11-dev \
	libxt-dev \
	python-numpy \
	libgd-dev \
	libxml2-dev \ 
	libxml2 \
	libpng12-0 \
	curl \
	texinfo \
	libcurl4-openssl-dev \
	libexpat1 \
	libexpat1-dev \ 
	libpangocairo-1.0-0 \ 
	libdb-dev \
	openjdk-7-jre \
	openjdk-7-jre-lib \
	git \
	bison \ 
	poxml \
	unzip \
	wget \
	graphviz && apt-get clean
#readline-devel 
#libgfortran  
#numpy 
#gd-devel 
#libxml2-devel 
#libxml2 
#libpng 
#libcurl-devel 
#db-devel 
#java 
#gettext-devel   
#db4-devel 
#libdb-devel 


# 
WORKDIR /opt
RUN git clone http://github.com/nunofonseca/irap.git  irap_clone  && cd /opt/irap_clone && ./scripts/irap_install.sh -a /irap_install -s . && cd / && rm -rf /irap_install/tmp /opt/irap_clone  /root/.cpan*

WORKDIR /
RUN echo source /irap_install/irap_setup.sh >> ~/.bash_profile && echo source /irap_install/irap_setup.sh >> ~/.bashrc

RUN echo '#!/usr/bin/env bash' > /usr/bin/irap
RUN echo 'source /irap_install/irap_setup.sh' >> /usr/bin/irap
RUN echo '/irap_install/scripts/irap "$@"' >> /usr/bin/irap
RUN chmod u+x /usr/bin/irap

#ENTRYPOINT ["irap"]

