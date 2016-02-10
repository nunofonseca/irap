#!/usr/bin/env bash
# =========================================================
# Copyright 2012-2016,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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
#    $Id$
# =========================================================
set -e

old_exp=$1
new_exp=$2
bam_op=$3

if [ "$bam_op-" = "-" ]; then
    echo "ERROR: usage: irap_clone.sh exp_name new_exp_name bam_operation@{addNH,rmNH,none}"
    exit 1
fi

if [ -e $new_exp ]; then
    echo "ERROR: target $new_exp already exists"
    exit 1
fi

if [ ! -e $old_exp ]; then
    echo "ERROR: source $old_exp not found"
    exit 1
fi

if [ "$bam_op" = "none" ]; then
    bam_op_cmd="ln -s"
else
if [ "$bam_op" = "addNH" ]; then
    bam_op_cmd=bam_add_NH
else
    bam_op_cmd=bam_rm_NH
fi
fi

pushd $old_exp
old_exp_path=`pwd`
popd
echo "Cloning $old_exp to $new_exp  ($bam_op)"
mkdir -p $new_exp
mkdir -p $new_exp/report

ln -s $old_exp_path/data $new_exp

# report
ln -s $old_exp_path/report/riq $new_exp/report/
ln -s $old_exp_path/report/read_filtering_plot.png $new_exp/report/

# Mappers
MAPPERS=`ls $old_exp_path/*/*.bam|less|sed -E "s|/[^/]*$||"|sed -E "s|.*/||"|sort -u`

for m in $MAPPERS; do
    echo "$m"
    mkdir -p $new_exp/$m
    for f in `ls $old_exp_path/$m/*.hits.bam`; do
	$bam_op_cmd $f $new_exp/$m/`basename $f`
    done
done




