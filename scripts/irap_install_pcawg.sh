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

DIR=$1
if [ "$DIR-" = "-" ]; then
    echo "You need to provide the instalation directory" > /dev/stderr
    exit 1
fi

if [ ! -e  scritps/irap_install.sh ]; then
    echo "This script should be executed from iRAP's installation directory"
    exit 1
fi

export IRAP_DIR=$DIR
./scripts/irap_install -s . -x make
./scripts/irap_install -s . -x gnuplot
./scripts/irap_install -s . -x YAP
./scripts/irap_install -s . -x samtools
./scripts/irap_install -s . -x bedtools
./scripts/irap_install -s . -x picard
./scripts/irap_install -s . -x core
source $DIR/irap_setup.sh
./scripts/irap_install -s . -x tophat2
./scripts/irap_install -s . -x bowtie2
./scripts/irap_install -s . -x star
./scripts/irap_install -s . -x fastq_qc
./scripts/irap_install -s . -v
exit 0

