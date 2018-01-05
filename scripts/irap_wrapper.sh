#!/bin/bash
# =========================================================
# Copyright 2012-2017,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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
#    $Id: irap.txt Nuno Fonseca Wed Jun 26 00:29:38 2013 +0100$
# =========================================================
# This wrapper sets up the environment before running an application
# it necessary to avoid program name conflicts (different versions of the same program)
set -e
APP=$1
shift 1
if [ "$IRAP_DIR-" = "-" ]; then
    echo "Missing env variable IRAP_DIR" >&2
    exit 1
fi

if [ "$APP-" = "-" ]; then
    echo "Error: Usage: irap_wrapper.sh app app_executable app_params..." >&2
    exit 1
fi

if [ ! -e $IRAP_DIR/bin/$APP ]; then
    echo "Unknown application $APP" >&2
    exit 1
fi

export LD_LIBRARY_PATH=$IRAP_DIR/bin/$APP/lib:$LD_LIBRARY_PATH
export PATH=$IRAP_DIR/bin/$APP:$IRAP_DIR/bin/$APP/bin:$IRAP_DIR/bin/:$PATH:$IRAP_DIR/bin/$APP

##########################################################
# special cases 
#if [ $APP = "rna_mate" ]; then
#    PERL5LIB=$PERL5LIB:$IRAP_DIR/bin/$APP/bin/perl
#fi
##########################################################
if [ -e /dev/stderr ]; then
    echo "Wrapper: $*" >&2
fi
$*
