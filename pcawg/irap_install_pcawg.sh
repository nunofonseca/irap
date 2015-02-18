#!/usr/bin/env bash
# =========================================================
# Copyright 2014,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
#
# This file is part of iRAP.
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with iRAP.  If not, see <http://www.gnu.org/licenses/>.
#
#
#   $Id: 0.1.3 Nuno Fonseca Fri Dec 21 17:18:21 2012$
# =========================================================
# Install step by step (assumes that dependencies are already installed in docker)



DIR=
SRC_DIR=.
OPTERR=0
while getopts "s:a:"  Option
do
    case $Option in
# update/reinstall
        a ) install=all;DIR=$OPTARG;;# send all output to a log file
        s ) SRC_DIR=$OPTARG;;# send all output to a log file
    esac
done

if [ "$DIR-" = "-" ]; then
    echo "You need to provide the instalation directory" > /dev/stderr
    exit 1
fi

if [ ! -e  scripts/irap_install.sh ]; then
    echo "This script should be executed from iRAP's installation directory"
    exit 1
fi
echo DIR=$DIR
set -e 
export IRAP_DIR=$DIR
./scripts/irap_install.sh -s $SRC_DIR  -x make
./scripts/irap_install.sh -s $SRC_DIR  -x zlib
./scripts/irap_install.sh -s $SRC_DIR  -x perl
./scripts/irap_install.sh -s $SRC_DIR  -x R3
./scripts/irap_install.sh -s $SRC_DIR -x gnuplot
./scripts/irap_install.sh -s $SRC_DIR -x YAP
./scripts/irap_install.sh -s $SRC_DIR -x samtools1
./scripts/irap_install.sh -s $SRC_DIR -x samtools
./scripts/irap_install.sh -s $SRC_DIR -x bedtools

./scripts/irap_install.sh -s $SRC_DIR -x core
source $DIR/irap_setup.sh
./scripts/irap_install.sh -s $SRC_DIR  -x perl_packages
./scripts/irap_install.sh -s $SRC_DIR  -x R3_packages
./scripts/irap_install.sh -s $SRC_DIR -x tophat2
./scripts/irap_install.sh -s $SRC_DIR -x bowtie2
./scripts/irap_install.sh -s $SRC_DIR -x bowtie1
./scripts/irap_install.sh -s $SRC_DIR -x star
./scripts/irap_install.sh -s $SRC_DIR -x fastq_qc
./scripts/irap_install.sh -s $SRC_DIR -x htseq
./scripts/irap_install.sh -s $SRC_DIR -v
# iRAP is now installed...
# install all the data needed by pcawg
#./examples/setup_pcawg.sh
exit 0

