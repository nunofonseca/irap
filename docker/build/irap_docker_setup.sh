#!/bin/bash

## 1 - OS (ubuntu|fedora)_version or skip_os_install
## 2 - irap version
## 3 - full/minimal installation (full/minimal)
## 4 - irap extra install options
set -eux
OS=$1
IRAP_VERSION=$2
INSTALL_TYPE=$3
shift 3
INSTALL_OPTIONS=$*

os_install_done=n
if [ "$OS-" == "fedora_27-" ]; then
    echo "Updating and installing Fedora packages"
    dnf update -y
    # R>=3.4.3
    dnf install -y --nodocs zlib-devel python-devel bzip2-devel python readline-devel libgfortran gcc-gfortran gcc-c++ libX11-devel libXt-devel numpy gd-devel libxml2-devel libxml2 libpng texi2html libcurl-devel expat-devel  pango-devel cairo-devel  java python gcc gcc-c++ gcc-objc++  gcc-gfortran curl git which make bzip2 bison gettext-devel  unzip make wget sqlite sqlite-devel db4-devel libdb-devel graphviz texlive tar java-devel texinfo texinfo-tex xorg-x11-server-Xvfb openssl-devel hdf5-devel bc
    # texlive-incons*
    dnf clean all
    os_install_done=y
    INSTALL_OPTIONS="$INSTALL_OPTIONS -R"
fi

if [ "$OS-" == "fedora_28-" ]; then
    echo "Updating and installing Fedora packages"
    dnf update -y
    # R>=3.4.3
    dnf install -y --nodocs R zlib-devel python-devel bzip2-devel python readline-devel libgfortran gcc-gfortran gcc-c++ libX11-devel libXt-devel numpy gd-devel libxml2-devel libxml2 libpng texi2html libcurl-devel expat-devel  pango-devel cairo-devel  java python gcc gcc-c++ gcc-objc++  gcc-gfortran curl git which make bzip2 bison gettext-devel  unzip make wget sqlite sqlite-devel db4-devel libdb-devel graphviz texlive tar java-devel texinfo texinfo-tex xorg-x11-server-Xvfb  openssl-devel hdf5-devel bc
    # texlive-incons*
    dnf clean all
    os_install_done=y
fi

# fedora 22
if [ "$OS-" == "fedora_22-" ]; then
    echo "Updating and installing Fedora packages"
    dnf update -y
    dnf install -y --nodocs zlib-devel python-devel bzip2-devel python readline-devel libgfortran gcc-gfortran gcc-c++ libX11-devel libXt-devel numpy gd-devel libxml2-devel libxml2 libpng texi2html libcurl-devel expat-devel  pango-devel cairo-devel  java python gcc gcc-c++ gcc-objc++  gcc-gfortran curl git which make bzip2 bison gettext-devel  unzip make wget sqlite sqlite-devel db4-devel libdb-devel graphviz texlive tar java-devel texinfo  xorg-x11-server-Xvfb automake evince  texlive-collection-latex firefox tk-devel tcl-devel libtool
    # hdf5 is not installed
    # texlive-incons*  
    dnf clean all
    os_install_done=y
    INSTALL_OPTIONS="$INSTALL_OPTIONS -R"
fi

if [ "$OS-" == "ubuntu_14-" ]; then
    echo "Updating and installing Ubuntu packages"
    apt-get update
    apt-get install -y build-essential texlive xvfb zlibc zlib1g zlib1g-dev libncurses5-dev sqlite sqlite3 libsqlite3-dev gettext python python-dev gfortran bzip2 libbz2-1.0 libbz2-dev libreadline-dev libx11-dev libxt-dev python-numpy libgd-dev libxml2-dev libxml2 libpng12-0 curl texinfo libcurl4-openssl-dev libexpat1 libexpat1-dev libpangocairo-1.0-0 libdb-dev openjdk-7-jre openjdk-7-jre-lib git bison poxml unzip wget graphviz libpcre3-dev libpcre++-dev bc libtbb-dev  xfonts-75dpi
    apt-get clean
    os_install_done=y
    # 
    INSTALL_OPTIONS="$INSTALL_OPTIONS -R"
fi
if [ "$OS-" == "ubuntu_16-" ] || [ "$OS-" == "ubuntu-" ]; then
    ## ubuntu
    echo "Updating and installing Ubuntu packages"
    #apt-get install -y
    ## R >= 3.4.2
    if [ "$IRAP_VERSION-" != "devel-" ]; then
	apt-get update
    fi
    if [ "$OS-" == "ubuntu-" ]; then
	apt-get update
	apt-get install -yq  --no-install-recommends software-properties-common ca-certificates wget curl ssh
	apt-get install -yq --no-install-recommends make
	make --version
    fi
    
    apt-get install -yq --no-install-recommends build-essential
    apt-get install -yq --no-install-recommends unzip
    apt-get install -yq --no-install-recommends libpcre3-dev
    apt-get install -yq --no-install-recommends libssl-dev
    apt-get install -yq --no-install-recommends lsb-release
    apt-get install -yq --no-install-recommends curl
    apt-get install -yq --no-install-recommends bison bzip2 curl gettext gfortran git graphviz libboost-all-dev libbz2-1.0  libbz2-dev libcurl4-openssl-dev libdb-dev libexpat1 libexpat1-dev libgd-dev libncurses5-dev libpangocairo-1.0-0 libpcre++-dev libpng-dev libreadline-dev libsqlite3-dev  libx11-dev libxml2 libxml2-dev libxt-dev lsb-release openjdk-8-jre openjdk-8-jdk poxml python-dev python-numpy sqlite sqlite3 texinfo texlive  wget xvfb zlib1g zlib1g-dev zlibc libpthread-stubs0-dev  libhdf5-dev bc libtbb-dev  xfonts-75dpi
    echo "deb http://cran.rstudio.com/bin/linux/ubuntu/ $(lsb_release -cs)/" >> /etc/apt/sources.list
    if [ "$OS-" == "ubuntu_16-" ]; then
	gpg --keyserver hkp://keyserver.ubuntu.com --recv-key E084DAB9    0
	gpg -a --export E084DAB9 | apt-key add -
	## update needed before installing R
	apt-get update
	apt-get install -y --no-install-recommends  r-base r-base-dev
    else
	INSTALL_OPTIONS="$INSTALL_OPTIONS -R"
    fi
    # clean up
    apt-get clean
    os_install_done=y
fi
if [ "$OS-" == "aws-" ]; then
    #yum update -y
    yum install -y make zlib-devel python-devel bzip2-devel python readline-devel libgfortran gcc-gfortran gcc-c++ libX11-devel libXt-devel numpy gd-devel libxml2-devel libxml2 libpng texi2html libcurl-devel expat-devel  pango-devel cairo-devel  java python gcc gcc-c++ gcc-objc++  gcc-gfortran curl git which make bzip2 bison gettext-devel  unzip make wget sqlite sqlite-devel db4-devel libdb-devel graphviz texlive tar java-devel texinfo  xorg-x11-server-Xvfb automake evince  texlive-collection-latex firefox tk-devel tcl-devel libtool python3-devel python3 openssl-devel
    ## Note: When you install the EPEL rpm package on Amazon Linux 2, the EPEL repository is already enabled.
    yum install -y https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm
    yum install -y hdf5-devel
    # texlive-incons*  
    yum clean all
    os_install_done=y
    INSTALL_OPTIONS="$INSTALL_OPTIONS -R"
fi
if [ "$OS-" == "skip_os_install-" ]; then
    os_install_done=y
fi

if [ $os_install_done == n ]; then
    echo ERROR: invalid OS parameter $OS > /dev/stderr
    exit 1
fi
mkdir -p /opt/
cd opt

#
CLONE=n
if [ "$IRAP_VERSION-" == "devel-" ]; then
    CLONE=y
    CLONE_BRANCH=devel
fi
if [ "$IRAP_VERSION-" == "master-" ]; then
    CLONE=y
    CLONE_BRANCH=master
fi
if [ "$IRAP_VERSION-" == "irap_new_release-" ]; then
    CLONE=y
    CLONE_BRANCH=irap_new_release
fi
if [ $CLONE == y ]; then
    git clone http://github.com/nunofonseca/irap.git -b $CLONE_BRANCH irap_clone --depth 50
    rm -rf irap_clone/.git
else
    wget https://github.com/nunofonseca/irap/archive/v$IRAP_VERSION.tar.gz -O irap.tar.gz
    tar xzvf irap.tar.gz
    rm -f irap.tar.gz
    mv irap-$IRAP_VERSION /opt/irap_clone
fi

cd /opt/irap_clone
INSTALL_PARAM=
if [ "$INSTALL_TYPE-" =="minimal" ]; then
    INSTALL_PARAM=-l
else
    INSTALL_PARAM=-a
fi

IRAP_DEST_DIR=/opt/irap

if [ "$INSTALL_TYPE-" == "pcawg-" ]; then
    ./pcawg/irap_install_pcawg.sh -a $IRAP_DEST_DIR -s . 
else
## install irap
# full path
    ./scripts/irap_install.sh -s . $INSTALL_PARAM $IRAP_DEST_DIR $INSTALL_OPTIONS
fi

cd /
## cleanup
rm -rf $IRAP_DEST_DIR/tmp /opt/irap_clone  /root/.cpan* 

##
cat <<EOF > /usr/bin/irap
#!/usr/bin/env bash
source $IRAP_DEST_DIR/irap_setup.sh
xvfb-run $IRAP_DEST_DIR/scripts/irap "\$@"
EOF
chmod u+x /usr/bin/irap

cat <<EOF > /usr/bin/irap_sc
#!/usr/bin/env bash
source $IRAP_DEST_DIR/irap_setup.sh
xvfb-run $IRAP_DEST_DIR/scripts/irap_sc "\$@"
EOF
chmod u+x /usr/bin/irap_sc
exit 0
