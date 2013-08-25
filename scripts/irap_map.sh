#!/bin/bash
# =========================================================
# Copyright 2012,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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
#    $Id: irap.txt Nuno Fonseca Tue Jan 8 12:29:27 2013$
# =========================================================
# This wrapper sets up the environment before running a mapper
# it necessary to avoid program name conflicts (different versions of the same program)
set -e
MAPPER=$1
shift 1
CMD=$*
if [ "$IRAP_DIR-" = "-" ]; then
    echo "Missing env variable IRAP_DIR"
    exit 1
fi

if [ "$MAPPER-" = "-" ]; then
    echo "Error: Usage: irap_map mapper_name mapper_executable mapper_params..."
    exit 1
fi

if [ ! -e $IRAP_DIR/bin/$MAPPER ]; then
    echo "Unknown mapper $MAPPER"
    exit 1
fi

export LD_LIBRARY_PATH=$IRAP_DIR/bin/$MAPPER/lib:$LD_LIBRARY_PATH
export PATH=$IRAP_DIR/bin/$MAPPER:$IRAP_DIR/bin/$MAPPER/bin:$IRAP_DIR/bin/:$PATH:$IRAP_DIR/bin/$MAPPER

##########################################################
# special cases (dependencies between mappers
if [ $MAPPER = "tophat1" ]; then
    export PATH=$IRAP_DIR/bin/bowtie1/bin:$PATH
fi
if [ $MAPPER = "bwa2" ]; then
    export LD_LIBRARY_PATH=$IRAP_DIR/bin/bwa/lib:$LD_LIBRARY_PATH
    export PATH=$IRAP_DIR/bin/bwa/bin:$PATH
fi
if [ $MAPPER = "tophat2" ]; then
    export PATH=$IRAP_DIR/bin/bowtie2/bin:$PATH
fi
if [ $MAPPER = "bismark" ]; then
    export PATH=$IRAP_DIR/bin/bowtie2/bin:$IRAP_DIR/bin/bowtie1/bin:$PATH
fi
if [ $MAPPER = "rna_mate" ]; then
    PERL5LIB=$PERL5LIB:$IRAP_DIR/bin/$MAPPER/bin/perl
#perl -MCPAN -e shell
#install ForkManager
#KWILLIAMS/Path-Class-0.26.tar.gz
fi
if [ $MAPPER = "gsnap" ]; then
    export PATH=$IRAP_DIR/bin/gmap/bin:$IRAP_DIR/bin/gmap/bin:$PATH
fi
if [ $MAPPER = "stampy" ]; then
    export PATH=$IRAP_DIR/bin/bwa/bin:$IRAP_DIR/bin/bwa/bin:$PATH
fi
if [ $MAPPER = "osa" ]; then
    CMD="mono $IRAP_DIR/bin/osa/bin/$CMD"
fi
if [ $MAPPER = "gnumap" ]; then
# ADD PATH2MPI
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/openmpi/lib/
    export PATH=$PATH:/usr/lib64/openmpi/bin/
fi
if [ $MAPPER = "passion" ]; then
    export PATH=$IRAP_DIR/bin/smalt/bin:$PATH
    alias smalt='smalt_x86_64'
fi

if [ $MAPPER = "bs_seeker" ]; then
    export PATH=$IRAP_DIR/bin/bowtie1/bin:$PATH:$IRAP_DIR/bin/bowtie1/bin
    if [ ! -e $IRAP_DIR/bin/bs_seeker/bin/bs_seeker_index ]; then
	cat <<EOF  >  $IRAP_DIR/bin/bs_seeker/bin/bs_seeker_index
#!/bin/bash 
python $IRAP_DIR/bin/bs_seeker/bin/Preprocessing_genome.py \$*
EOF
	chmod +x $IRAP_DIR/bin/bs_seeker/bin/bs_seeker_index
    fi
    if [ ! -e $IRAP_DIR/bin/bs_seeker/bin/bs_seeker ]; then
	cat <<EOF  >  $IRAP_DIR/bin/bs_seeker/bin/bs_seeker
#!/bin/bash 
python $IRAP_DIR/bin/bs_seeker/bin/BS_Seeker.py \$*
EOF
	chmod +x $IRAP_DIR/bin/bs_seeker/bin/bs_seeker
    fi
    if [ ! -e $IRAP_DIR/bin/bs_seeker/bin/bs_seeker_out2sam ]; then
	cat <<EOF  >  $IRAP_DIR/bin/bs_seeker/bin/bs_seeker_out2sam
#!/bin/bash 
python $IRAP_DIR/bin/bs_seeker/bin/BSSout2SAM.py \$*
EOF
	chmod +x $IRAP_DIR/bin/bs_seeker/bin/bs_seeker_out2sam
    fi

fi
##########################################################
if [ -e /dev/stderr ]; then
    echo "Wrapper: $CMD" >&2
fi
$CMD
