#!/bin/bash

set -eux

IRAP_OPT=$1

apt-get update
apt-get install -yq --no-install-recommends build-essential
apt-get install -yq --no-install-recommends texlive
apt-get install -yq --no-install-recommends xvfb
apt-get install -yq --no-install-recommends zlibc
apt-get install -yq --no-install-recommends zlib1g
apt-get install -yq --no-install-recommends zlib1g-dev
apt-get install -yq --no-install-recommends libncurses5-dev
apt-get install -yq --no-install-recommends sqlite
apt-get install -yq --no-install-recommends libsqlite3-dev
apt-get install -yq --no-install-recommends gettext
apt-get install -yq --no-install-recommends python-dev
apt-get install -yq --no-install-recommends gfortran
apt-get install -yq --no-install-recommends libbz2-dev
apt-get install -yq --no-install-recommends libreadline-dev
apt-get install -yq --no-install-recommends libx11-dev
apt-get install -yq --no-install-recommends libxt-dev
apt-get install -yq --no-install-recommends python-numpy
apt-get install -yq --no-install-recommends libgd-dev
apt-get install -yq --no-install-recommends libxml2-dev
apt-get install -yq --no-install-recommends texinfo
apt-get install -yq --no-install-recommends libcurl4-openssl-dev
apt-get install -yq --no-install-recommends libpangocairo-1.0-0
apt-get install -yq --no-install-recommends libdb-dev
apt-get install -yq --no-install-recommends openjdk-8-jre
apt-get install -yq --no-install-recommends git
apt-get install -yq --no-install-recommends bison
apt-get install -yq --no-install-recommends poxml
apt-get install -yq --no-install-recommends wget
apt-get install -yq --no-install-recommends graphviz
apt-get install -yq --no-install-recommends unzip
apt-get install -yq --no-install-recommends libpcre3-dev
apt-get install -yq --no-install-recommends libssl-dev
apt-get install -yq --no-install-recommends lsb-release
apt-get install -yq --no-install-recommends curl

# install R
UBUNTU_VER=`lsb_release -cs`

echo "deb http://cran.rstudio.com/bin/linux/ubuntu $UBUNTU_VER/" >> /etc/apt/sources.list
gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
gpg -a --export E084DAB9 | apt-key add -
apt-get install -yq --no-install-recommends r-base r-base-dev

apt-get clean

curl -sSL https://github.com/nunofonseca/irap/archive/v0.8.5.p8.tar.gz > irap.tar.gz
tar xzf irap.tar.gz
./irap-0.8.5.p8/scripts/irap_install.sh -a $IRAP_OPT -s irap-0.8.5.p8
# clean??
rm -r irap-0.8.5.p8 irap.tar.gz $IRAP_OPT/tmp .cpan*
