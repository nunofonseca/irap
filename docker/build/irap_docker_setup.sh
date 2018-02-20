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
    dnf install -y zlib-devel python-devel bzip2-devel python readline-devel libgfortran gcc-gfortran gcc-c++ libX11-devel libXt-devel numpy gd-devel libxml2-devel libxml2 libpng texi2html libcurl-devel expat-devel  pango-devel cairo-devel  java python gcc gcc-c++ gcc-objc++  gcc-gfortran curl git which make bzip2 bison gettext-devel  unzip make wget sqlite sqlite-devel db4-devel libdb-devel graphviz texlive tar java-devel texinfo texinfo-tex xorg-x11-server-Xvfb texlive-incons*
    dnf clean all
    os_install_done=y
fi

# fedora 20
if [ "$OS-" == "fedora_22-" ]; then
    echo "Updating and installing Fedora packages"
    dnf update -y
    dnf install -y R zlib-devel python-devel bzip2-devel python readline-devel libgfortran gcc-gfortran gcc-c++ libX11-devel libXt-devel numpy gd-devel libxml2-devel libxml2 libpng texi2html libcurl-devel expat-devel  pango-devel cairo-devel  java python gcc gcc-c++ gcc-objc++  gcc-gfortran curl git which make bzip2 bison gettext-devel  unzip make wget sqlite sqlite-devel db4-devel libdb-devel graphviz texlive tar java-devel texinfo  xorg-x11-server-Xvfb texlive-incons*  automake evince  texlive-collection-latex firefox tk-devel tcl-devel libtool
    dnf clean all
    os_install_done=y
fi

if [ "$OS-" == "ubuntu_14-" ]; then
    echo "Updating and installing Ubuntu packages"
    apt-get update
    apt-get install -y build-essential texlive xvfb zlibc zlib1g zlib1g-dev libncurses5-dev sqlite sqlite3 libsqlite3-dev gettext python python-dev gfortran bzip2 libbz2-1.0 libbz2-dev libreadline-dev libx11-dev libxt-dev python-numpy libgd-dev libxml2-dev libxml2 libpng12-0 curl texinfo libcurl4-openssl-dev libexpat1 libexpat1-dev libpangocairo-1.0-0 libdb-dev openjdk-7-jre openjdk-7-jre-lib git bison poxml unzip wget graphviz
    apt-get clean
    os_install_done=y

fi
if [ "$OS-" == "ubuntu_16-" ]; then
    ## ubuntu
    echo "Updating and installing Ubuntu packages"
    apt-get update
    apt-get install -y  build-essential texlive xvfb zlibc zlib1g zlib1g-dev  libncurses5-dev  sqlite  sqlite3 libsqlite3-dev  gettext  python  python-dev   gfortran bzip2 libbz2-1.0 libbz2-dev libreadline-dev  libx11-dev libxt-dev python-numpy r-base libgd-dev  libxml2-dev  libxml2 libpng12-0 curl texinfo libcurl4-openssl-dev libexpat1 libexpat1-dev libpangocairo-1.0-0  libdb-dev openjdk-8-jre openjdk-8-jdk git bison  poxml  unzip  libboost-all-dev wget graphviz
    apt-get clean
    os_install_done=y
fi

if [ "$OS-" == "skip_install_os-" ]; then
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
if [ "$IRAP_VERSION-" == "devel" ]; then
    CLONE=y
    CLONE_BRANCH=devel
fi
if [ "$IRAP_VERSION-" == "irap_new_release-" ]; then
    CLONE=y
    CLONE_BRANCH=irap_new_release
fi
if [ $CLONE == y ]; then
    git clone http://github.com/nunofonseca/irap.git -b $CLONE_BRANCH irap_clone --depth 50
else
    wget https://github.com/nunofonseca/irap/archive/v$IRAP_VERSION.tar.gz -o irap.tar.gz
    tar xzvf irap.tar.gz
    rm -f irap.tar.gz
    mv irap-$IRAP_VERSION irap_clone
fi

cd /opt/irap_clone
INSTALL_PARAM=
if [ "$INSTALL_TYPE-" =="minimal" ]; then
    INSTALL_PARAM=-l
else
    INSTALL_PARAM=-a
fi
## install irap
# full path
IRAP_DEST_DIR=/opt/irap
./scripts/irap_install.sh -s . $INSTALL_PARAM $IRAP_DEST_DIR -R $INSTALL_OPTIONS
cd /
## cleanup
rm -rf $IRAP_DEST_DIR/tmp /opt/irap_clone  /root/.cpan* 

##
cat <<EOF > /usr/bin/irap
#!/usr/bin/env bash
source $IRAP_DEST_DIR/irap_setup.sh
$IRAP_DEST_DIR/scripts/irap "$@"
EOF
chmod u+x /usr/bin/irap
