#!/usr/bin/env bash
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
#   $Id: 0.1.3 Nuno Fonseca Fri Dec 21 17:18:21 2012$
# =========================================================
install=all
IRAP_DIR1=
SRC_DIR=

IRAP_VERSION=0.8.5.p7


#
USE_CACHE=y
#############
function pinfo {
    echo "[INFO] $*"
}

function usage {
    echo "Usage: irap_install.sh  -s irap_toplevel_directory [ -c dir -a dir -a -u -m -r -p -q -b -K -R -B]  ";
    echo " -s dir : toplevel irap clone directory";
    echo " -c dir : install/update IRAP core files only to directory 'dir'. If dir is not given the value of IRAP_DIR will be used (if available).";
    echo " -a dir : install/update all files (core and 3rd party software) to directory 'dir' (default)";
    echo " -l dir : lightweight/minimal installation of iRAP (a minimum set of tools will be installed).";    
    echo " -u : update IRAP_core files (alias to -c $IRAP_DIR).";
    echo " -m : update mappers.";
    echo " -r : update R packages.";
    echo " -p : update Perl packages.";
    echo " -q : update quantifiers.";
    echo " -b : update jbrowser.";
    echo " -j : install jbrowser (with -a).";
    echo " -v : collect software versions.";
    echo " -G : install gcc 4.8 before installing Mono (GCC will be installed in \$IRAP_DIR/gcc).";
    echo " -B : install boost libraries.";
    echo " -K : use ksh instead of bash (due to an issue trapping signals) while installing some components (R).";
    echo " -R : install R.";
    echo " Advanced options:";
    echo " -f : check/fix file permissions"
    echo " -d : download all software and libraries (except R and Perl packages) but do not install.";
    echo " -x software: install/update software.";
}

function check_for_irap_env {

    if [ "$IRAP_DIR-" == "-" ]; then
	echo "ERROR: IRAP_DIR environment variable not defined"
	exit 1
    fi
    # check if IRAP_DIR is in the path
    IP="$IRAP_DIR/scripts" 
    if [[ ":$PATH:" != *":$IP:"* ]]; then
	echo "ERROR: IRAP_DIR ($IP) not in PATH."
	exit 1
    fi
}

function download {
    FILE2DOWNLOAD=$1
    FILE2=$2
    
    if [ "$FILE2-" == "-" ]; then
	FILE2=`basename $FILE2DOWNLOAD`
    fi
    
    if [ "$USE_CACHE-" == "y-" ]; then
	# avoid copying
	rm -f $FILE2
	ln -s $SRC_DIR/download/$FILE2 .
    else
	set +e
	if [ $OS == "linux" ]; then	    
	    wget  --no-check-certificate -c -nc -L "$FILE2DOWNLOAD" -O $FILE2
	else
	    curl $FILE2DOWNLOAD
	fi
	set -e
    fi
}
############################################
# download everything to the download folder
function download2cache {

    USE_CACHE=n
    mkdir -p $SRC_DIR/download
    pushd $SRC_DIR/download
    PACKAGES2DOWNLOAD=`set  | grep _URL=|sed "s/_URL.*//"|grep -v "PACKAGES"|uniq`
    echo $PACKAGES2DOWNLOAD > /dev/stderr
    for p in $PACKAGES2DOWNLOAD; do
	URL=${p}_URL
	FILE=${p}_FILE
	pinfo "Downloading ($p)  ${!URL}"
	download_software $p
	if [ ! -e ${!FILE} ]; then
	    echo "Failed downloading $p ${!URL}" > /dev/stderr
	    exit 1
	else
	    ls -lh ${!FILE}
	fi
    done
    pinfo "Downloaded files in $SRC_DIR/download"
}

function check_dependencies {
    DEVEL_LIBRARIES_REQUIRED="zlib-devel python-devel bzip2-devel python readline-devel libgfortran gcc-gfortran gcc-c++ libX11-devel libXt-devel numpy gd-devel libxml2-devel libxml2 libpng libcurl-devel expat-devel  libpangocairo bison gettext-devel  sqlite-devel sqlite [db-devel|db4-devel|libdb-devel] R"
    MISSING=0
    pinfo "If installation fails then please check if the following libraries are installed:"
    pinfo "$DEVEL_LIBRARIES_REQUIRED"
    # Binaries that should be available
    # make is required to...compile make
    BINARIES="java python gcc g++ gfortran curl-config git which make bzip2 unzip bash R"
    pinfo "Checking dependencies..."
    for bin in $BINARIES; do
	PATH2BIN=`which $bin 2> /dev/null`
	if [ "$PATH2BIN-" == "-" ]; then
	    pinfo " $bin not found!"
	    #
	    if [ "$bin" == "R" ]; then
		pinfo "WARNING:Please install the R package (and update all packages) or run irap_install.sh with -R to install R (version 3.2 or above)"
	    else
		MISSING=1
	    fi
	else
	    pinfo " $bin found: $PATH2BIN"
	fi
    done
    pinfo "Checking dependencies...done."
    if [ $MISSING == 1 ]; then
	pinfo "ERROR: Unable to proceed"
	exit 1
    fi

}
#################################
# VERSIONS, SRC file and URL
BFAST_VERSION=0.7.0a
BFAST_FILE=bfast-$BFAST_VERSION.tar.gz
BFAST_URL=http://sourceforge.net/projects/bfast/files/bfast/0.7.0/$BFAST_FILE
# current 1.1.2 - minor changes, no need to upgrade
bowtie1_VERSION=1.1.1
bowtie1_FILE=bowtie-${bowtie1_VERSION}-linux-x86_64.zip
bowtie1_URL=http://sourceforge.net/projects/bowtie-bio/files/bowtie/$bowtie1_VERSION/$bowtie1_FILE
# current - 2.2.6->2.2.9
bowtie2_VERSION=2.2.9
bowtie2_FILE=bowtie2-${bowtie2_VERSION}-linux-x86_64.zip
bowtie2_URL=http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$bowtie2_VERSION/$bowtie2_FILE

#
GEM_VERSION=20130406-045632
GEM_FILE=GEM-binaries-Linux-x86_64-core_2-${GEM_VERSION}.tbz2
GEM_URL=http://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%203/$GEM_FILE/download
#GEM_VERSION=1.6-i3
#GEM_FILE=GEM-gemtools-${GEM_VERSION}.tar.gz
#GEM_URL=http://barnaserver.com/gemtools/$GEM_FILE
#GEM_VERSION=20121106-022124
#GEM_FILE=GEM-binaries-Linux-x86_64-core_i3-$GEM_VERSION.tbz2
#GEM_URL=http://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%202/$GEM_FILE/download
#
tophat1_VERSION=1.4.1
tophat1_FILE=tophat-${tophat1_VERSION}.Linux_x86_64.tar.gz
tophat1_URL=http://ccb.jhu.edu/software/tophat/downloads/$tophat1_FILE

#
#tophat2_VERSION=2.0.13
#tophat2_FILE=tophat-${tophat2_VERSION}.Linux_x86_64.tar.gz
#tophat2_URL=http://tophat.cbcb.umd.edu/downloads/$tophat2_FILE
#tophat2_URL=http://ccb.jhu.edu/software/tophat/downloads/$tophat2_FILE

# version used in the pcawg SOP - 2.0.12
# current version: -> 2.0.12->2.1.0
tophat2_VERSION=2.1.1
tophat2_FILE=tophat-${tophat2_VERSION}.Linux_x86_64.tar.gz
#tophat2_URL=http://tophat.cbcb.umd.edu/downloads/$tophat2_FILE
tophat2_URL=http://ccb.jhu.edu/software/tophat/downloads/$tophat2_FILE

# 
HISAT2_VERSION=2.0.5
HISAT2_FILE=hisat2-${HISAT2_VERSION}-Linux_x86_64.zip
HISAT2_URL=ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/$HISAT2_FILE


# 0.7.6 - now in sf
SMALT_VERSION=0.7.4
SMALT_FILE=smalt-$SMALT_VERSION.tgz
SMALT_URL=ftp://ftp.sanger.ac.uk/pub4/resources/software/smalt/$SMALT_FILE
#
SOAPsplice_VERSION=1.10
SOAPsplice_FILE=SOAPsplice-v$SOAPsplice_VERSION.tar.gz
SOAPsplice_URL=http://soap.genomics.org.cn/down/$SOAPsplice_FILE
#
SOAP2_VERSION=2.21
SOAP2_FILE=soap${SOAP2_VERSION}release.tar.gz
SOAP2_URL=http://soap.genomics.org.cn/down/$SOAP2_FILE
# 2.4.0i->2.5.0c
STAR_VERSION=2.5.0c
STAR_FILE=${STAR_VERSION}.tar.gz
STAR_URL=https://github.com/alexdobin/STAR/archive/$STAR_FILE

# 
#GSNAP_VERSION=2013-11-27->gmap-gsnap-2015-12-31.v3.tar.gz
GSNAP_VERSION=2015-12-31
GSNAP_FILE=gmap-gsnap-${GSNAP_VERSION}.v3.tar.gz
GSNAP_URL=http://research-pub.gene.com/gmap/src/$GSNAP_FILE

# 2.2.0
mapsplice_VERSION=2.2.0
mapsplice_FILE=MapSplice-v$mapsplice_VERSION.zip
mapsplice_URL=http://protocols.netlab.uky.edu/~zeng/$mapsplice_FILE
# 0.7.4->0.7.12
bwa_VERSION=0.7.12
bwa_FILE=bwa-${bwa_VERSION}.tar.bz2
bwa_URL=http://sourceforge.net/projects/bio-bwa/files/$bwa_FILE
# 4.0.2.1->4.1.1.1
osa_VERSION=4.0.2.1
osa_FILE=OSAv$osa_VERSION.zip
osa_URL=http://omicsoft.com/osa/Software/$osa_FILE

# deprecated
#EMBAM_VERSION=0.1.14
#EMBAM_FILE=emBAM_${EMBAM_VERSION}.tar.gz
#EMBAM_URL=http://embam.googlecode.com/files/$EMBAM_FILE

RUBY_VERSION=1.9.3-p484
RUBY_FILE=ruby-${RUBY_VERSION}.tar.gz
RUBY_URL=http://ftp.ruby-lang.org/pub/ruby/1.9/$RUBY_FILE

# 
#PERL_VERSION=5.20.1 -> 5.22.1
PERL_VERSION=5.20.3
PERL_FILE=perl-$PERL_VERSION.tar.gz
PERL_URL=http://www.cpan.org/src/5.0/$PERL_FILE

# previous: 1.55
BOOST_VERSION=1.55.0
BOOST_FILE=boost_`echo $BOOST_VERSION|sed "s/\./_/g"`.tar.bz2
BOOST_URL=http://sourceforge.net/projects/boost/files/boost/$BOOST_VERSION/$BOOST_FILE

gnuplot_VERSION=4.6.4   
gnuplot_FILE=gnuplot-$gnuplot_VERSION.tar.gz
gnuplot_URL=http://sourceforge.net/projects/gnuplot/files/gnuplot/$gnuplot_VERSION/$gnuplot_FILE

# deprecated
R2_VERSION=2.15.2
R2_FILE=R-${R2_VERSION}.tar.gz 
R2_URL=http://cran.r-project.org/src/base/R-2/$R2_FILE

# 
R3_VERSION=3.2.5
R3_FILE=R-${R3_VERSION}.tar.gz 
R3_URL=http://cran.r-project.org/src/base/R-3/$R3_FILE

# 
SAMTOOLS_VERSION=0.1.18
SAMTOOLS_FILE=samtools-$SAMTOOLS_VERSION.tar.bz2
SAMTOOLS_URL=http://sourceforge.net/projects/samtools/files/samtools/$SAMTOOLS_VERSION/$SAMTOOLS_FILE

# new samtools: 1.3
# prev: 1.1
SAMTOOLS1_VERSION=1.3.1
SAMTOOLS1_FILE=samtools-$SAMTOOLS1_VERSION.tar.bz2
SAMTOOLS1_URL=http://sourceforge.net/projects/samtools/files/samtools/$SAMTOOLS1_VERSION/$SAMTOOLS1_FILE

VCFTOOLS_VERSION=0.1.14
VCFTOOLS_FILE=vcftools-$VCFTOOLS_VERSION.tar.gz
VCFTOOLS_URL=https://github.com/vcftools/vcftools/releases/download/v$VCFTOOLS_VERSION/$VCFTOOLS_FILE


# 1.3
BCFTOOLS_VERSION=1.3
BCFTOOLS_FILE=bcftools-$BCFTOOLS_VERSION.tar.bz2
BCFTOOLS_URL=https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VERSION/$BCFTOOLS_FILE


ZLIB_VERSION=1.2.8
ZLIB_FILE=zlib-$ZLIB_VERSION.tar.gz
ZLIB_URL=http://zlib.net/$ZLIB_FILE

#
BEDTOOLS_VERSION=2.25.0
BEDTOOLS_FILE=bedtools-$BEDTOOLS_VERSION.tar.gz
BEDTOOLS_URL=https://github.com/arq5x/bedtools2/releases/download/v$BEDTOOLS_VERSION/$BEDTOOLS_FILE

# 1.0.3->1.2.0
stringtie_VERSION=1.2.0
stringtie_FILE=stringtie-${stringtie_VERSION}.Linux_x86_64.tar.gz
stringtie_URL=http://ccb.jhu.edu/software/stringtie/dl/$stringtie_FILE

cufflinks1_VERSION=1.3.0
cuffdiff1_VERSION=1.3.0
cufflinks1_FILE=cufflinks-${cufflinks1_VERSION}.Linux_x86_64.tar.gz
cufflinks1_URL=http://cufflinks.cbcb.umd.edu/downloads/$cufflinks1_FILE

cufflinks2_VERSION=2.2.1
cuffdiff2_VERSION=2.2.1
cufflinks2_FILE=cufflinks-${cufflinks2_VERSION}.Linux_x86_64.tar.gz
cufflinks2_URL=http://cufflinks.cbcb.umd.edu/downloads/$cufflinks2_FILE

#IReckon_VERSION=1.0.6
#IReckon_FILE=IReckon-$IReckon_VERSION.jar
#IReckon_URL=http://compbio.cs.toronto.edu/ireckon/$IReckon_FILE

#SAVANT_VERSION=2.0.3
#SAVANT_FILE=Savant-${SAVANT_VERSION}-Linux-x86_64-Install 
#SAVANT_URL=http://genomesavant.com/savant/dist/v`echo $SAVANT_VERSION|sed "s/./_/g"`/$SAVANT_FILE

# new: 0.7.5 - now in github
BitSeq_VERSION=0.7.0
BitSeq_FILE=BitSeq-$BitSeq_VERSION.tar.gz
BitSeq_URL=http://bitseq.googlecode.com/files/$BitSeq_FILE

#MMSEQ_VERSION=1.0.0-beta2
#MMSEQ_FILE=mmseq_${MMSEQ_VERSION}.zip
#MMSEQ_URL=http://www.bgx.org.uk/software/$MMSEQ_FILE

#htseq_VERSION=0.5.4p5
htseq_VERSION=0.6.1p1
htseq_FILE=HTSeq-${htseq_VERSION}.tar.gz
htseq_URL=http://pypi.python.org/packages/source/H/HTSeq/$htseq_FILE

# new: 1.6.1
# latest (1.5.2) fails (same error since 1.2.4)
# Parameter COVERAGE_STATS not found. Check the spelling!
FLUX_CAPACITOR_VERSION=1.2.3
FLUX_CAPACITOR_FILE=flux-capacitor-$FLUX_CAPACITOR_VERSION.tgz
#FLUX_CAPACITOR_URL=http://sammeth.net/artifactory/barna-nightly/barna/barna.capacitor/$FLUX_CAPACITOR_VERSION/$FLUX_CAPACITOR_FILE
FLUX_CAPACITOR_URL=http://sammeth.net/artifactory/barna/barna/barna.capacitor/$FLUX_CAPACITOR_VERSION/$FLUX_CAPACITOR_FILE

#
NURD_VERSION=1.1.1
NURD_FILE=NURD_v${NURD_VERSION}.tar.gz
NURD_URL=http://bioinfo.au.tsinghua.edu.cn/software/NURD/share/$NURD_FILE

#
IsoEM_VERSION=1.1.4
IsoEM_FILE=IsoEM-${IsoEM_VERSION}.zip
IsoEM_URL=http://dna.engr.uconn.edu/software/IsoEM/$IsoEM_FILE

# 0.6.2-> 0.9.0
Sailfish_VERSION=0.9.0
Sailfish_FILE=v${Sailfish_VERSION}.tar.gz
Sailfish_URL=http://github.com/kingsfordgroup/sailfish/archive/$Sailfish_FILE

SALMON_VERSION=0.7.2
SALMON_FILE=Salmon-${SALMON_VERSION}_linux_x86_64.tar.gz
SALMON_URL=https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/${SALMON_FILE}

# kallisto - 0.42.1 -> 0.42.4 (multithreads)
kallisto_VERSION=0.42.4
kallisto_FILE=kallisto_linux-v$kallisto_VERSION.tar.gz
kallisto_URL=https://github.com/pachterlab/kallisto/releases/download/v$kallisto_VERSION/$kallisto_FILE

# 1.2.21->1.2.22
rsem_VERSION=1.2.22
rsem_FILE=rsem-${rsem_VERSION}.tar.gz
rsem_URL=http://deweylab.biostat.wisc.edu/rsem/src/$rsem_FILE

FUSIONMAP_VERSION=2015-03-31
FUSIONMAP_FILE=FusionMap_${FUSIONMAP_VERSION}.zip
FUSIONMAP_URL=http://omicsoft.com/fusionmap/Software/$FUSIONMAP_FILE


SCRIPTURE_VERSION=beta2
SCRIPTURE_FILE=scripture-${SCRIPTURE_VERSION}.jar 
SCRIPTURE_URL=ftp://ftp.broadinstitute.org/pub/papers/lincRNA/$SCRIPTURE_FILE

IGV_VERSION=2.1.24
IGV_FILE=IGV_$IGV_VERSION.zip
IGV_URL=http://www.broadinstitute.org/igv/projects/downloads/$IGV_FILE

IGV_TOOLS_VERSION=2.1.24
IGV_TOOLS_FILE=igvtools_nogenomes_$IGV_TOOLS_VERSION.zip
IGV_TOOLS_URL=http://www.broadinstitute.org/igv/projects/downloads/$IGV_TOOLS_FILE

# new: FASTX_VERSION=0.0.13->0.0.13.2
FASTX_VERSION=0.0.13
FASTX_FILE=fastx_toolkit_${FASTX_VERSION}_binaries_Linux_2.6_amd64.tar.bz2
FASTX_URL=http://hannonlab.cshl.edu/fastx_toolkit/$FASTX_FILE

# 0.10.1->0.11.4
FASTQC_VERSION=0.11.4
FASTQC_FILE=fastqc_v${FASTQC_VERSION}.zip
FASTQC_URL=http://www.bioinformatics.babraham.ac.uk/projects/fastqc/$FASTQC_FILE

# 1.12.0
JBROWSE_VERSION=1.7.5
JBROWSE_FILE=download.php?id=35
JBROWSE_URL=http://jbrowse.org/wordpress/wp-content/plugins/download-monitor/$JBROWSE_FILE
JBROWSE_EXTRA_UTILS="hgGcPercent bedGraphToBigWig wigCorrelate bigWigInfo bigWigSummary faToNib faToTwoBit hgWiggle bigWigMerge"
JBROWSE_EXTRA_UTILSURL=http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/

# 1.12.0 -- 85
NEW_JBROWSE_VERSION=1.12.0
NEW_JBROWSE_FILE=download.php?id=101
NEW_JBROWSE_URL=http://jbrowse.org/wordpress/wp-content/plugins/download-monitor/$NEW_JBROWSE_FILE
NEW_JBROWSE_EXTRA_UTILS="hgGcPercent bedGraphToBigWig wigCorrelate bigWigInfo bigWigSummary faToNib faToTwoBit hgWiggle bigWigMerge"
NEW_JBROWSE_EXTRA_UTILSURL=http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/

# 2.10.8 osa does not work with 2.10.9 up to 2.11
# 2.10.8 ---> 4.2.2
MONO_VERSION=2.10.8
MONO_FILE=mono-${MONO_VERSION}.tar.bz2    
#MONO_VERSION=4.2.2
#MONO_FILE=mono-${MONO_VERSION}.tar.bz2    
#MONO_FILE=mono-${MONO_VERSION}.30.tar.bz2    
MONO_URL=http://download.mono-project.com/sources/mono/$MONO_FILE


MAKE_VERSION=4.1
MAKE_FILE=make-${MAKE_VERSION}.tar.gz
MAKE_URL=http://ftp.gnu.org/gnu/make/$MAKE_FILE
##################################
function download_software {
    name=$1
    case $name in
	bfast ) download $BFAST_URL $BFAST_FILE;;
	star ) download $STAR_URL $STAR_FILE;;
	gsnap ) download $GSNAP_URL $GSNAP_FILE;;
	smalt ) download $SMALT_URL $SMALT_FILE;;
	soap_splice ) download $SOAPsplice_URL;;
	soap2 ) download $SOAP2_URL;download http://soap.genomics.org.cn/down/soap2sam.tar.gz;;
	SOAP2 ) download $SOAP2_URL;download http://soap.genomics.org.cn/down/soap2sam.tar.gz;;
	JBROWSE ) download $JBROWSE_URL $JBROWSE_FILE; for f in $JBROWSE_EXTRA_UTILS; do  download $JBROWSE_EXTRA_UTILSURL/$f $f; done ;;
	NEW_JBROWSE ) download $NEW_JBROWSE_URL $NEW_JBROWSE_FILE; for f in $NEW_JBROWSE_EXTRA_UTILS; do  download $NEW_JBROWSE_EXTRA_UTILSURL/$f $f; done ;;
	* ) url=${name}_URL; file=${name}_FILE; download ${!url} ${!file};;
    esac
}

# install all program files (binaries, libraries,...) in a folder of the bin directory
# 
function install_binary {
    PROGNAME=$1
    SRCDIR=$2
    TDIR=$BIN_DIR/$PROGNAME/bin
    shift 2
    FILES=$*

    pushd $SRCDIR
    mkdir -p $TDIR    
    chmod +x $FILES
    cp -rf $FILES $TDIR
    popd
}

function gen_setup_irap {
    cat <<EOF > $1
export IRAP_DIR=$IRAP_DIR
export PATH=\$IRAP_DIR/bin/bowtie1/bin:\$IRAP_DIR/bin:\$IRAP_DIR/scripts:\$IRAP_DIR/python/bin/:\$PATH
export LD_LIBRARY_PATH=\$IRAP_DIR/lib:\$LD_LIBRARY_PATH:/usr/local/lib
export CFLAGS="-I\$IRAP_DIR/include -I\$IRAP_DIR/include/bam -I\$IRAP_DIR/include/boost  \$CFLAGS"
export R_LIBS_USER=$IRAP_DIR/Rlibs
export R_LIBS=$IRAP_DIR/Rlibs
export CXXFLAGS="-I\$IRAP_DIR/include -I\$IRAP_DIR/include/bam -I\$IRAP_DIR/include/boost -L\$IRAP_DIR/lib \$CXXFLAGS"
export PERL5LIB=\$IRAP_DIR/perl/lib/perl5:\$IRAP_DIR/lib/perl5:\$IRAP_DIR/lib/perl5/x86_64-linux:\$IRAP_DIR/lib/perl5/$PERL_VERSION
export PYTHONUSERBASE=\$IRAP_DIR/python
# Adjust to your needs and uncomment the following lines
# in order to use iRAP with the LSF job scheduler
# note: memory values are in MB 
#export IRAP_LSF_GROUP=/irap
#export QUEUE=research-rh6
export MEM=10000
export THREADS=8
#export JOB_MEM_INCR 8000
#export JOB_MAX_MEM 32000
#export IRAP_LSF_PARAMS=
EOF
    if [ "$INSTALL_R3-" == "y-" ]; then
	echo "export R_LIBS=	" >> $1
	echo "PATH=\$IRAP_DIR/R3/bin:\$PATH" >> $1
    fi
    mkdir -p $IRAP_DIR/Rlibs
}


####################################
# Mappers
function bfast_install {
    MAPPER=bfast
    pinfo "Starting $MAPPER source installation..."
    download_software $MAPPER
    tar xzvf $BFAST_FILE
    pushd bfast-$BFAST_VERSION
    ./configure --prefix=$BIN_DIR/$MAPPER
    make clean
    make
    make install
    pinfo "$MAPPER installation complete."
    popd
}

function bowtie1_install {
    MAPPER=bowtie1
    pinfo "Starting $MAPPER source installation..."
    download_software $MAPPER
    unzip $bowtie1_FILE
    pushd bowtie-$bowtie1_VERSION
    #export BITS=64
    #make clean
    #make -j 4 all
    FILES="bowtie bowtie-*"
    install_binary $MAPPER . $FILES 
    install_binary $MAPPER scripts \*
    pinfo "$MAPPER installation complete."
    popd
}

function bowtie2_install {
    MAPPER=bowtie2
    pinfo "Starting $MAPPER source installation..."
    download_software $MAPPER
    unzip $bowtie2_FILE
    pushd `echo $bowtie2_FILE|sed "s/-linux-x86_64.zip//"`
    #export BITS=64
    #make clean
    #make -j $J all
    FILES="bowtie2 bowtie2-build* bowtie2-inspect* bowtie2-align*"
    install_binary $MAPPER . $FILES 
    install_binary $MAPPER scripts \*
    pinfo "$MAPPER installation complete."
    popd
}

function tophat1_install {
    MAPPER=tophat1
    pinfo  "Starting $MAPPER binary installation..."
    download_software $MAPPER
    tar xzvf $tophat1_FILE
    mkdir -p $BIN_DIR/tophat1
    pushd `echo $tophat1_FILE|sed "s/.tar.gz//"`
    install_binary $MAPPER . \*
    pinfo "$MAPPER installation complete."    
    popd
}

function tophat2_install {
    MAPPER=tophat2
    pinfo  "Starting $MAPPER binary installation..."
    download_software $MAPPER
    tar xzvf $tophat2_FILE
    pushd `echo $tophat2_FILE|sed "s/.tar.gz//"`
    install_binary $MAPPER . \*
    cp $IRAP_DIR/bin/tophat2/bin/gtf_juncs $IRAP_DIR/bin/tophat2_gtf_juncs
    pinfo "$MAPPER installation complete."    
    popd
}

function hisat2_install {
    MAPPER=HISAT2
    pinfo  "Starting $MAPPER binary installation..."
    download_software $MAPPER
    unzip $HISAT2_FILE
    pushd `echo $HISAT2_FILE|sed "s/.Linux.*//"`
    install_binary $MAPPER . \*
    cp -rf scripts/* $IRAP_DIR/bin/$MAPPER/bin
    pinfo "$MAPPER installation complete."    
    popd
}


function smalt_install {
    MAPPER=smalt
    pinfo "Starting $MAPPER binary installation..."
    download_software $MAPPER
    tar xzvf $SMALT_FILE
    pushd `echo $SMALT_FILE|sed "s/.tgz//"`
    install_binary $MAPPER .  smalt_x86_64
    pinfo "$MAPPER installation complete."    
    popd
}

function soap_splice_install {
    MAPPER=soap_splice
    SRCDIR=$2    
    pinfo "Starting $MAPPER binary installation..."
    download_software $MAPPER
    tar xzvf $SOAPsplice_FILE
    pushd `echo $SOAPsplice_FILE|sed "s/.tar.gz//"`
    install_binary $MAPPER bin  \*
    pinfo "$MAPPER installation complete."    
    popd
}

function soap2_install {
    MAPPER=soap2
    pinfo "Starting $MAPPER binary installation..."
    download_software $MAPPER
    tar xzvf $SOAP2_FILE
    pushd `echo $SOAP2_FILE|sed "s/.tar.gz//"`
    install_binary $MAPPER .  2bwt-builder  soap
    popd
    #download http://soap.genomics.org.cn/down/soap2sam.tar.gz
    tar xzvf soap2sam.tar.gz
    chmod +x soap2sam.pl
    sed "s.^#\!/usr/bin/perl.*.#\!$ENV_FP perl." -i soap2sam.pl
    cp soap2sam.pl $BIN_DIR
    pinfo "$MAPPER installation complete."    
}

function gem_install {
    MAPPER=GEM
    pinfo "Starting $MAPPER binary installation..."
    download_software $MAPPER
    tar xjvf $GEM_FILE
    # deps: requires ruby
    ruby_install
    pushd `echo $GEM_FILE|sed "s/.tbz2//"`
    install_binary $MAPPER bin \*
#    sed -i "s/^#!/.*ruby/#!$ENV_FP ruby/" $BIN_DIR/$MAPPER/bin/gem*
    pinfo "$MAPPER  installation complete."    
    popd
}


function gem2_install {
    MAPPER=GEM
    pinfo "Starting $MAPPER binary installation..."
    download_software $MAPPER
    # do not install gemtools - only the gem binaries
    tar xjvf $GEM_FILE
    install_binary $MAPPER . gem-* gtf* splits* compute* transcri*
    pinfo "$MAPPER  installation complete."    
    popd
}

function star_install {
    MAPPER=star
    pinfo  "Starting $MAPPER binary installation..."
    download_software $MAPPER    
    #gunzip -c $STAR_FILE > $EXECF
    tar xzvf $STAR_FILE
    pushd STAR-$STAR_VERSION
    mkdir -p $BIN_DIR/star/bin
    rm -rf $BIN_DIR/star/bin/star/* $BIN_DIR/star/bin/STAR/*
    cp bin/Linux_x86_64_static/STAR $BIN_DIR/star/bin/star
    cp bin/Linux_x86_64_static/STAR $BIN_DIR/star/bin/STAR
#    cp STARstatic $BIN_DIR/star/bin/star
#    cp STARstatic $BIN_DIR/star/bin/
    pinfo "$MAPPER installation complete."    
    popd
}

function gsnap_install {
    # gmap/gsnap
    MAPPER=gsnap
    pinfo "Starting $MAPPER source installation..."
    # latest stable version
    download_software $MAPPER
    tar xzvf $GSNAP_FILE
    pushd `echo $GSNAP_FILE|sed "s/.tar.gz//"|sed "s/-gsnap//;s/.v3.*//"`
    ./configure --prefix=$BIN_DIR/$MAPPER --with-samtools=$IRAP_DIR
    make clean
    make -j $J
    make install
    pinfo "$MAPPER installation complete."    
    popd
}

function bwa_install {
    MAPPER=bwa
    pinfo "Starting $MAPPER source installation..."
    download_software $MAPPER
    tar xjvf $bwa_FILE
    pushd `echo $bwa_FILE|sed "s/.tar.bz2//"`
    make clean
    make -j $J
    install_binary $MAPPER . bwa
    pinfo "$MAPPER installation complete."    
    popd
}



function osa_install {
    MAPPER=osa
    pinfo "Starting $MAPPER binary installation..."
    OLD_PATH=$PATH
    OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH	
    if [ $INSTALL_GCC != "n" ]; then
	gcc4_install
	export PATH=$IRAP_DIR/gcc/bin:$PATH
	export LD_LIBRARY_PATH=$IRAP_DIR/gcc/lib64:$LD_LIBRARY_PATH
    fi
    mono_install
    download_software $MAPPER
    unzip $osa_FILE
    pushd OSAv$osa_VERSION
    install_binary $MAPPER . \*
    PATH=$PATH
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH	
    popd
    pinfo "$MAPPER installation complete."    
}

function mapsplice_install {
    MAPPER=mapsplice
    pinfo "Starting $MAPPER source installation..."
    download_software $MAPPER
    unzip $mapsplice_FILE
    pushd MapSplice-v$mapsplice_VERSION
    make clean
    # do not recompile bowtie and samtools
    sed -i -E "1s/all:(.*)bowtie/all:\1 /" Makefile
    sed -i -E "1s/all:(.*)samtools/all:\1 /" Makefile
    make
    install_binary $MAPPER bin \*
    cp mapsplice.py $BIN_DIR/$MAPPER/
    popd
    pinfo "$MAPPER installation complete."    
}

function mappers_install {
   bwa_install
   bowtie1_install
   bowtie2_install
   tophat1_install
   tophat2_install
   #bfast_install
   smalt_install
   soap_splice_install
   soap2_install
#   gem2_install
   gsnap_install 
   ## osa_install
   star_install
   pinfo "To install MapSplice run: irap_install.sh -s . -x mapsplice"
   pinfo "To install GEM run: irap_install.sh -s . -x gem"
   pinfo "To install HISAT2 run: irap_install.sh -s . -x hisat2"
   #mapsplice_install
}

########################
function embam_install {
        
    R=$1
    if [ "$R-" == "-" ];  then
	R=$IRAP_DIR/bin/R
    fi
    pinfo "Starting emBAM source installation..."
    download_software EMBAM
    $R  CMD INSTALL $EMBAM_FILE
    pinfo "Starting emBAM source installation...done."
}
function  collect_software_versions {
    pinfo "Collecting software versions..."
    irap_gen_version_mk.sh
    pinfo "Collecting software versions...done."
}

######################################################
#                   Dependencies
######################################################

######################################################
# Make
# Utility determining automatically which pieces of a large program need to be recompiled, and issues commands to recompile them
function make_install {
    pinfo "Installing make..."
    download_software MAKE
    tar xzvf $MAKE_FILE
    pushd `echo $MAKE_FILE|sed "s/.tar.gz//"`
    ./configure 
    make -j $J
    cp make $IRAP_DIR/bin
    popd
    pinfo "Installing make...done."
}

####################################################################
# It may be necessary-mono 2.x does not compile successfully in gcc
# 4.9 or higher
function gcc4_install {

    pinfo "Installing gcc 4.8.5..."
    wget -c ftp://ftp.mirrorservice.org/sites/sourceware.org/pub/gcc/releases/gcc-4.8.5/gcc-4.8.5.tar.gz
    tar xzvf gcc-4.8.5.tar.gz
    pushd gcc-4.8.5
    ./contrib/download_prerequisites
    cd ..
    mkdir objdir
    cd objdir
    CFLAGS=  $PWD/../gcc-4.8.5/configure --prefix=$IRAP_DIR/gcc --enable-languages=c,c++,fortran --disable-multilib
    make
    make install
    popd
    pinfo "Installing gcc 4.8.5...done."
}
 
######################################################
# 
function mono_install {
    pinfo "Installing mono..."
    download_software MONO
    tar xjvf $MONO_FILE
    pushd mono-$MONO_VERSION
    ./configure --with-large-heap=yes --prefix=$IRAP_DIR
    make
    make install
    popd    
    pinfo "Installing mono...done."
}
    
######################################################
# Ruby    
# Programming language (mmseq includes several scripts in ruby)
function ruby_install {
    pinfo "Installing ruby..."
    # previous 1.9.3-p125
    download_software RUBY
    tar xzvf $RUBY_FILE
    pushd ruby-${RUBY_VERSION}
    ./configure --prefix=$IRAP_DIR
    make 
    make install
    popd
    pinfo "Installing ruby...done."
}
######################################################
# Perl 
# To ensure that there are no version issues...
function perl_install {
    pinfo "Installing perl..."
    if [ -e  ~/.cpan ]; then
	if [ "$INIT_CPAN-" != "no-" ]; then
	    pinfo "Moving ~/.cpan to ~/.cpan.bak"
	    cp -arb ~/.cpan ~/.cpan.bak
	    rm -fr ~/.cpan
	else
	    pinfo "Skipping CPAN initialization (you may enable it by doing export INIT_CPAN=yes)"
	fi
    fi
    download_software PERL
    tar -xzf $PERL_FILE
    pushd perl-$PERL_VERSION
    ./Configure -des -Dprefix=$IRAP_DIR
    make -j $J
    # skip tests...
    set +e
    make test
    set -e
    make install 
    popd
    pinfo "Installing perl...done."
}
######################################################
# BOOST
function boost_install {
    pinfo "Installing boost..."
    download_software BOOST
    tar xvjf $BOOST_FILE
    pushd `echo $BOOST_FILE|sed "s/.tar.bz2//"`
    ./bootstrap.sh --prefix=$IRAP_DIR
    ./b2 -q
    ./b2 install
    popd
    pinfo "Installing boost...done."
}

######################################################
# gnuplot 4.4
# Utility the for visual display of scientific data
function gnuplot_install {
    pinfo "Installing gnuplot..."
    download_software gnuplot
    tar xzvf $gnuplot_FILE
    pushd `echo $gnuplot_FILE|sed "s/.tar.gz//"`
    ./configure --prefix=`pwd` --without-readline --without-lisp-files
    make
    make install
    cp bin/gnuplot $IRAP_DIR/bin
    popd
    pinfo "Installing gnuplot...done."
}

######################################################
# R
function R2_install {
    pinfo "Installing R2..."
    download_software R2
    tar xzvf $R2_FILE
    pushd R-${R2_VERSION}
    export R_LIBS=
    export R_LIBS_USER=$IRAP_DIR/Rlibs
    # assume that makeinfo is installed - configure does not work on 5.2
    #sed -i "s/r_cv_prog_makeinfo_v4=no/r_cv_prog_makeinfo_v4=yes/" configure  
    CFLAGS_noboost=`echo $CFLAGS|sed -E "s|\-I[^ ]*boost||g"`    
    # clean up - delete packages previously installed
    rm -rf $IRAP_DIR/Rlibs/*
    CFLAGS=$CFLAGS_noboost $SPECIAL_SH_TO_USE ./configure --prefix=$IRAP_DIR
    #fedora 23:--disable-nls
    CFLAGS=$CFLAGS_noboost make clean
    CFLAGS=$CFLAGS_noboost make
    CFLAGS=$CFLAGS_noboost make check
    CFLAGS=$CFLAGS_noboost make install
    popd
    pinfo "Installing R2...done."
}

# install R-3.x
function R_install {
    pinfo "Installing R-3.x..."
    download_software R3
    tar xzvf $R3_FILE
    pushd R-${R3_VERSION}
    export R_LIBS=
    export R_LIBS_USER=$IRAP_DIR/Rlibs3
    CFLAGS_noboost=`echo $CFLAGS|sed -E "s|\-I[^ ]*boost||g"`    
    CFLAGS=$CFLAGS_noboost $SPECIAL_SH_TO_USE ./configure --prefix=$IRAP_DIR
    make clean
    make -j $J
    make -j $J check
    make install
    popd
    # wrappers
#     cat <<EOF > $IRAP_DIR/scripts/R3
# #!/bin/bash
# export PATH=\$IRAP_DIR/R3/bin:\$PATH
# export R_LIBS_USER=\$R3_LIBS_USER
# \$IRAP_DIR/R3/bin/R "\$@"
# EOF
#     chmod +x $IRAP_DIR/scripts/R3
#     cat <<EOF > $IRAP_DIR/scripts/Rscript3
# #!/bin/bash
# export PATH=\$IRAP_DIR/R3/bin:\$PATH
# export R_LIBS_USER=\$R3_LIBS_USER
# \$IRAP_DIR/R3/bin/Rscript "\$@"
# EOF
#     chmod +x $IRAP_DIR/scripts/Rscript3

    pinfo "Installing R3-x...done."
}

######################################################
# Yap 
function YAP_install {
    pinfo "Installing YAP..."
    if [ $USE_CACHE == "y" ]; then
	tar xzvf $SRC_DIR/download/yap.tgz
    else
	rm -rf mainline
	git clone https://github.com/vscosta/yap-6.3
	#http://gitorious.org/yap-git/mainline.git
	mv yap-6.3 mainline
	rm -rf mainline/.git
	tar czvf yap.tgz mainline
	cp yap.tgz $SRC_DIR/download/
    fi
    pushd mainline
    # cleanup
    rm -rf $IRAP_DIR/include/Yap $IRAP_DIR/share/Yap $IRAP_DIR/lib/Yap $IRAP_DIR/lib/libYap*
    CXXFLAGS= CFLAGS=  ./configure --prefix=$IRAP_DIR --disable-myddas --disable-horus
    CXXFLAGS=  CFLAGS=  make clean
    CXXFLAGS=  CFLAGS=  make
    CXXFLAGS=  CFLAGS=  make install
    popd
    pinfo "Installing YAP...done."
}

function deps_install {

    pinfo "Installing dependencies (make, perl, boost, gnuplot, R, samtools, ...)"
    make_install
    zlib_install
    perl_install
    #ruby_install
    if [ "$1-" != "minimal-" ] && [ "$BOOST_INSTALL" == "y" ]; then
	boost_install
    fi
    gnuplot_install
    #R_install
    if [ "$INSTALL_R3-" == "y-" ]; then
	R_install
    fi
    YAP_install
    # some mappers (e.g., tophat) still require samtools 0.x
    samtools_install
    samtools1_install
    bedtools_install
    vcftools_install
    python_packages_install
    ucsc_utils_install
    #picard_install
    pinfo "Installing dependencies...done."
}

######################################################
# Samtools
# Utilities for manipulating alignments in the SAM format
function samtools_install {
    pinfo "Downloading, compiling, and installing SamTools..."
    download_software SAMTOOLS
    tar xvjf $SAMTOOLS_FILE
    pushd samtools-${SAMTOOLS_VERSION}
    make -j $J
    make -j $J razip
    mkdir -p $BIN_DIR/samtools0.x/bin
    cp samtools razip bcftools/vcfutils.pl bcftools/bcftools $BIN_DIR/samtools0.x/bin
    mkdir -p $BIN_DIR/samtools0.x/lib
    mkdir -p $BIN_DIR/samtools0.x/include/bcftools
    cp *.h $BIN_DIR/samtools0.x/include
    cp bcftools/*.h $BIN_DIR/samtools0.x/include/bcftools
    cp bcftools/libbcf.a libbam.a $BIN_DIR/samtools0.x/lib
    popd
    pinfo "Downloading, compiling, and installing SAMTools...done."
}

function samtools1_install {
    pinfo "Downloading, compiling, and installing SamTools 1.x..."
    download_software SAMTOOLS1
    tar xvjf $SAMTOOLS1_FILE
    pushd samtools-${SAMTOOLS1_VERSION}
    make -j $J prefix=$IRAP_DIR
    make prefix=$IRAP_DIR install
    mkdir -p $INC_DIR/bam
    cp *.h $INC_DIR/bam
    mkdir  -p $INC_DIR/bam/htslib-1.3.1
    cp htslib-1.3.1/*.h $INC_DIR/bam/htslib-1.3.1
    cp libbam.a $INC_DIR/bam
    #
    pwd
    download_software BCFTOOLS
    tar xjvf $BCFTOOLS_FILE
    pushd bcftools-${BCFTOOLS_VERSION}
    sed -i -E "s|^prefix\s*=.*|prefix=$IRAP_DIR|"  Makefile
    make -j $J 
    make install
    popd
    popd
    pinfo "Downloading, compiling, and installing SAMTools 1.x...done."
}

function vcftools_install {
    pinfo "Downloading, compiling, and installing VCFTOOLS..."
    download_software VCFTOOLS
    tar xvzf $VCFTOOLS_FILE
    pushd vcftools-${VCFTOOLS_VERSION}
    ./configure prefix=$IRAP_DIR
    make -j $J prefix=$IRAP_DIR
    make prefix=$IRAP_DIR install
    popd
    pinfo "Downloading, compiling, and installing VCFTOOLS...done."       
}

######################################################
# zlib
function zlib_install {
    pinfo "Downloading, compiling, and installing zlib..."
    download_software ZLIB
    tar xvzf $ZLIB_FILE
    pushd zlib-${ZLIB_VERSION}
    ./configure --prefix $IRAP_DIR
    make 
    make install
    popd
    pinfo "Downloading, compiling, and installing zlib...done."
}

######################################################
# Bedtools
function bedtools_install {
    pinfo "Installing BEDTOOLS..."
    BEDTOOLS=bedtools-$BEDTOOLS_VERSION
    download_software BEDTOOLS
    tar xzvf $BEDTOOLS_FILE
    pushd bedtools2
    make
    cp bin/* $IRAP_DIR/bin
    popd
    pinfo "Installing BEDTOOLS...done."
}

######################################################
# Perl packages
# TODO: move from cpan to cpanm
function perl_cpan_install {
    if [ -e $IRAP_DIR/.cpan.irap.done ]; then
	pinfo "Skipping cpan init...already done"
    else
    pinfo "Initializing CPAN..."
    pinfo "! Internet access required !"
    if  [ ! -e ~/.cpan/ ]; then
	unset PERL5LIB
    else
	if [ "$INIT_CPAN-" != "no-" ]; then
	    cp -rab ~/.cpan  ~/.orig.cpan
	    rm -rf ~/.cpan
	fi
    fi
    # only reset cpan if the user wants it
    if [ "$INIT_CPAN-" != "no-" ]; then
	(echo y;) | perl -MCPAN -e shell
    else
	pinfo "Skipping CPAN initialization (you may enable it by doing export INIT_CPAN=yes)"
    fi
    # if myConfig.pm existed the previous command would not change it therefore force its creation
    ( echo mkmyconfig; ) | perl -MCPAN -e shell 
    ( echo o conf init urllist;echo y;echo o conf commit;) | perl -MCPAN -e shell 
#    (echo o conf init urllist;echo y;echo 3;echo 31; echo 1 2 3 4 5 6;echo o conf commit;) | perl -MCPAN -e shell 
    pinfo "Initializing CPAN...done."
    pinfo "Configuring CPAN..."
    # extra configuration
    #PREFIX=$IRAP_DIR/perl
    export INSTALL_BASE=$IRAP_DIR
    export PREFIX=$IRAP_DIR
    perl -MCPAN -e shell<<EOF
o conf makepl_arg  INSTALL_BASE=$IRAP_DIR PREFIX=$IRAP_DIR
o conf mbuildpl_arg "--install_base $IRAP_DIR/perl"
o conf prerequisites_policy follow
o conf mbuild_install_arg "--install_base $IRAP_DIR/perl"
o conf build_requires_install_policy yes
o conf commit
q
EOF
# 
    # upgrade cpan
    #cpan autobundle    
    set +e
    cpan  -f -i App::cpanminus
    cpanm -l $IRAP_DIR -n -f -i YAML   < /dev/null
    #cpanm -f -i Test::More@0.99 < /dev/null
    cpanm -l $IRAP_DIR  -n -i -f ExtUtils::MakeMaker  < /dev/null 
    # perhaps install the latest perl?
    cpan -f -u 
    # don't test
    #cpanm -n -i  Bundle::CPAN
    cpanm $IRAP_DIR  -f -n -i  CPAN < /dev/null
    set -e
    # set permissions 
    chmod +w $IRAP_DIR/bin/*
    pinfo "Configuring CPAN...done."
    touch $IRAP_DIR/.cpan.irap.done
  fi
}

function perl_bundle_install {

    perl_cpan_install
	
    pinfo "Installing perl bundle..."
    mkdir -p ~/.cpan/Bundle
    # cpan -a 
    mv $SRC_DIR/aux/irap.pm  ~/.cpan/Bundle
    cpan Bundle::irap < /dev/null
    pinfo "Installing perl bundle...done."
}

# short list of packages need by iRAP code
function perl_packages_install {
    mkdir -p $IRAP_DIR/perl
    #
    # ensure that make is in path
    export PATH=$IRAP_DIR/bin:$PATH
    
    # TODO: do it only once
    perl_cpan_install
    pinfo "Installing perl packages..."
    #  required by iRAP core: Bio::SeqIO
    cpanm -l $IRAP_DIR  --force -n http://www.cpan.org/authors/id/L/LE/LEONT/Module-Build-0.40_11.tar.gz
    cpanm -l $IRAP_DIR -f -n  Bio::SeqIO
    #cpan -fi CJFIELDS/BioPerl-1.6.1.tar.gz  
    # for jbrowse see jbrowse install
    #cpanm -fi CJFIELDS/BioPerl-1.6.924.tar.gz 
    #popd
    pinfo "Installing perl packages...done."
}

function perl_packages_jbrowse_install {
    mkdir -p $IRAP_DIR/perl
    #
    # ensure that make is in path
    export PATH=$IRAP_DIR/bin:$PATH
    
    # TODO: do it only once
    perl_cpan_install
    
    pinfo "Installing perl packages for jbrowse..."
    #    Test::Requisites
    PACKAGES="Algorithm::Munkres     Array::Compare    Math::Random    Sort::Naturally    Sub::Install Sub::Uplevel     Params::Util    List::MoreUtils    Math::Round    DB_File    Test     Test::Fatal    Test::Run    Test::NoWarnings  Test::Exception  Error    XML::Parser    XML::Simple    XML::SAX    XML::SAX::Writer    XML::Writer JSON Hash::Merge  Devel::Size  PerlIO::locale Compress::Raw::Zlib Locale::Maketext::Lexicon Build Module::CoreList ExtUtils::MakeMaker"    

    set +e
    for p in $PACKAGES; do
       pinfo "************ Package $p"
       cpanm -l $IRAP_DIR --notest --force -n $p < /dev/null
    done
    set -e
    # the tests fail...
    cpanm -l $IRAP_DIR -n  -i Heap::Simple::Perl
    # 
    cpanm -l $IRAP_DIR -n  -i L/LD/LDS/GD-2.50.tar.gz
    # SAMTOOLS needs to be recompiled :(
    mkdir -p $IRAP_DIR/tmp
    pushd $IRAP_DIR/tmp
    download_software SAMTOOLS
    tar xvjf $SAMTOOLS_FILE
    #svn export https://samtools.svn.sourceforge.net/svnroot/samtools/tags/samtools-0.1.7/;
    perl -i -pe 's/^CFLAGS=\s*/CFLAGS=-fPIC / unless /\b-fPIC\b/' samtools-${SAMTOOLS_VERSION}/Makefile;
    make -C samtools-${SAMTOOLS_VERSION} -j3 lib;
    export SAMTOOLS="$PWD/samtools-${SAMTOOLS_VERSION}";    
    #cpanm -fi CJFIELDS/BioPerl-1.6.924.tar.gz 

    popd
    pinfo "Installing perl packages for jbrowse...done."
}

##################################
# R packages
# Software environment for statistical computing and graphics
function R2_packages_install {
    export PATH=$IRAP_DIR/bin:$PATH
    pinfo "Installing R packages..."
    #export R_LIBS=
    CFLAGS_noboost=`echo $CFLAGS|sed -E "s|\-I[^ ]*boost||g"`

    CFLAGS=$CFLAGS_noboost R --no-save <<EOF
repo<-"$CRAN_REPO"

source("http://bioconductor.org/biocLite.R")
biocLite("DBI",ask=FALSE,suppressUpdates=TRUE)
# RSQLite_1 does not compile
download.file("http://cran.r-project.org/src/contrib/Archive/RSQLite/RSQLite_0.11.4.tar.gz","RSQLite_0.11.4.tar.gz")
install.packages("RSQLite_0.11.4.tar.gz",type="source",repos=NULL)
download.file("http://cran.r-project.org/src/contrib/Archive/gplots/gplots_2.11.0.tar.gz","gplots_2.11.0.tar.gz")

biocLite("gtools",ask=FALSE,suppressUpdates=TRUE)
biocLite("gdata",ask=FALSE,suppressUpdates=TRUE)
biocLite("caTools",ask=FALSE,suppressUpdates=TRUE)
install.packages("gplots_2.11.0.tar.gz",type="source",repos=NULL)
download.file("http://cran.r-project.org/src/contrib/Archive/xtable/xtable_1.7-4.tar.gz","xtable_1.7-4.tar.gz")
install.packages("xtable_1.7-4.tar.gz",type="source",repos=NULL)
#install.packages("gplots",repo=repo)

# agricolae deps
download.file("http://cran.r-project.org/src/contrib/Archive/sp/sp_1.0-11.tar.gz","sp_1.0-11.tar.gz")
install.packages("sp_1.0-11.tar.gz",type="source",repos=NULL)

biocLite("maptools",ask=FALSE,suppressUpdates=TRUE)
biocLite("deldir",ask=FALSE,suppressUpdates=TRUE)
biocLite("coda",ask=FALSE,suppressUpdates=TRUE)
download.file("http://cran.r-project.org/src/contrib/Archive/spdep/spdep_0.5-40.tar.gz","spdep_0.5-40.tar.gz")
install.packages("spdep_0.5-40.tar.gz",type="source",repos=NULL)



manually.installed.packages<-c("^gplots_","^xtable","^RSQLite","^spdep_","^sp_")

# not available in 2.15.2: "plyr", "reshape","lattice"
packages2install<-c("multicore","parallel","intervals","gclus","R2HTML","agricolae","optparse","brew","gtools","gdata","caTools","sfsmisc","GO.db","edgeR","DESeq","Rsamtools","DEXSeq","baySeq","limma","marray","org.Hs.eg.db","goseq","data.table")
for ( p in packages2install ) {
   cat("Installing ",p,":\n")
   biocLite(p,ask=FALSE,  suppressUpdates=TRUE)
}

species2db<-matrix(c('org.Ag.eg.db','Anopheles',
'org.At.tair.db','Arabidopsis',
'org.Bt.eg.db','Bovine',
'org.Ce.eg.db','Worm',
'org.Cf.eg.db','Canine',
'org.Dm.eg.db','Fly',
'org.Dr.eg.db','Zebrafish',
'org.EcK12.eg.db','E coli strain K12',
'org.Gg.eg.db','Chicken',
'org.Hs.eg.db','Human',
'org.Mm.eg.db','Mouse',
'org.Mmu.eg.db','Rhesus',
'org.Pf.plasmo.db','Malaria',
'org.Pt.eg.db','Chimp',
'org.Rn.eg.db','Rat',
'org.Sc.sgd.db','Yeast',
'org.Sco.eg.db','Streptomyces coelicolor',
'org.Ss.eg.db','Pig',
'org.Tgondii.eg.db','Toxoplasma gondii',
'org.Xl.eg.db','Xenopus'),byrow=T,ncol=2)
colnames(species2db)<-c("db","species")
for (p in species2db[,'db']) {
  biocLite(p,ask=FALSE,  suppressUpdates=TRUE)                     
}

q()
EOF
    # deprecated
    #pinfo "Installing EMBAM..."
    #embam_install $IRAP_DIR/bin/R
    #pinfo "installing EMBAM...done."
    pinfo "Installing R packages...done."
}
# 
# requires libcurl installed in the system

function R_packages_install {
    export PATH=$IRAP_DIR/bin:$PATH
    pinfo "Installing R-3.x packages..."

    CFLAGS_noboost=`echo $CFLAGS|sed -E "s|\-I[^ ]*boost||g"`
    # suppressUpdates should be TRUE otherwise it might try to update a package installed in the systems folder
    CFLAGS=$CFLAGS_noboost echo y | $SRC_DIR/scripts/install_R_packages.R
    pinfo "Installing R-3.x packages...done."
}
######################################################
# Stringtie
function stringtie_install {

# Short reads - Transcript assembly, abundance and differential expression estimations
    pinfo "Downloading and installing StringTie..."
    download_software stringtie
    tar xzvf $stringtie_FILE
    # file name conflict with cufflinks2
    mkdir -p $BIN_DIR/stringtie
    cp stringtie-${stringtie_VERSION}*/* $BIN_DIR/stringtie
    pinfo "Downloading and installing StringTie...done."
}

######################################################
# Cufflinks1
function cufflinks1_install {

# Short reads - Transcript assembly, abundance and differential expression estimations
    pinfo "Downloading and installing CuffLinks..."
    download_software cufflinks1
    tar xzvf $cufflinks1_FILE
    # file name conflict with cufflinks2
    mkdir -p $BIN_DIR/cufflinks1
    cp cufflinks-${cufflinks1_VERSION}*/* $BIN_DIR/cufflinks1
    pinfo "Downloading and installing conflinks1...done."
}

######################################################
# Cufflinks2
function cufflinks2_install {
# Short reads - Transcript assembly, abundance and differential expression estimations     
    pinfo "Downloading and installing CuffLinks 2..."
    download_software cufflinks2
    tar xzvf $cufflinks2_FILE
    # file name conflict with cufflinks1
    mkdir -p $BIN_DIR/cufflinks2
    cp `echo $cufflinks2_FILE|sed "s/.tar.gz//"`/* $BIN_DIR/cufflinks2
    pinfo "Downloading and installing cufflinks2...done."
}
#######################################################
#bitseq
#Bayesian Inference of Transcripts from Sequencing Data
function bitseq_install {
    pinfo "Installing BitSeq..."
    download_software BitSeq
    tar xzvf $BitSeq_FILE
    pushd BitSeq-$BitSeq_VERSION
    make
    #PROGS=`grep "^PROGRAMS =" Makefile| cut -f 2 -d=`     
    mkdir -p $IRAP_DIR/bin/bitseq/bin
    find .  -maxdepth 1 -executable -type f -exec  mv {}  $IRAP_DIR/bin/bitseq/bin \; ;
    popd
    pinfo "Installing BitSeq...done."
}

######################################################
# mmseq
# Haplotype and isoform specific expression estimation using multi-mapping RNA-seq reads
function mmseq_install {
    pinfo "Installing mmseq..."
    download_software MMSEQ
    unzip $MMSEQ_FILE
    pushd mmseq_${MMSEQ_VERSION}/
    rm -f *Darwin*
    cp *.rb *-x86_64 *.sh  $IRAP_DIR/bin
    # remove version from the name of the executables
    rename -- "-${MMSEQ_VERSION}*-Linux-x86_64" "" $IRAP_DIR/bin/*-$MMSEQ_VERSION*
    popd
    pinfo "Installing mmseq...done."
}

######################################################
# HTSeq
function htseq_install {
    pinfo "Installing HTSeq..."
    download_software htseq
    tar xzvf $htseq_FILE
    pushd `echo $htseq_FILE|sed "s/.tar.gz//"`
# python version needs to be equal or greater than  (2.6)
    #. ./build_it ;# not needed in 0.5.4p5
    # python setup.py install --user
    pip install --user .
    chmod +x scripts/*
    cp scripts/* $IRAP_DIR/bin
    popd
    pinfo "Installing HTSeq...done."
}

######################################################
# Flux capacitor
function flux_capacitor_install {
    pinfo "Installing Flux_capacitor..."
# stable version
# download http://sammeth.net/artifactory/barna/barna/barna.capacitor/1.2.2/flux-capacitor-$FLUX_VERSION.tgz
# devel version
    download_software FLUX_CAPACITOR
#
    tar xzvf $FLUX_CAPACITOR_FILE
    pushd `echo $FLUX_CAPACITOR_FILE|sed "s/.tgz//"`
    mv bin/* $IRAP_DIR/bin
    mv lib/* $IRAP_DIR/lib
    # note: $IRAP_DIR/lib/ must be in the LD_LIBRARY_PATH variable? maybe not...
    popd
    pinfo "Installing Flux_capacitor...done."
}
######################################################
# IsoEM
function isoem_install {
    pinfo "Installing IsoEM..."
    download_software IsoEM
    unzip $IsoEM_FILE
    pushd isoem-$IsoEM_VERSION
    # fix the PATH
    sed -i "s|isoEMDir=.*|isoEMDir=$IRAP_DIR|" bin/*
    cp bin/* $IRAP_DIR/bin
    cp lib/* $IRAP_DIR/lib
    popd
    pinfo "IsoEM installation complete."    
}

######################################################
# Sailfish
# requires boost
function sailfish_install {
    pinfo "Installing Sailfish..."
    download_software Sailfish
    tar xzvf $Sailfish_FILE
    pushd sailfish-$Sailfish_VERSION
    echo "Not working..."
    exit 2
    rm -rf $IRAP_DIR/bin/Sailfish/lib/*
    rm -rf $IRAP_DIR/bin/Sailfish/bin/*
    mkdir -p $IRAP_DIR/bin/Sailfish/lib
    mkdir -p $IRAP_DIR/bin/Sailfish/bin
    mv bin/* $IRAP_DIR/bin/Sailfish/bin
    mv lib/* $IRAP_DIR/bin/Sailfish/lib
    popd
    pinfo "Sailfish installation complete."    
}


######################################################
# rsem
function rsem_install {
    pinfo "Installing rsem..."
    download_software rsem
    tar xzvf $rsem_FILE
    pushd rsem-$rsem_VERSION
    make
    mkdir -p $IRAP_DIR/bin/rsem/bin
    cp rsem* extract-* convert-* $IRAP_DIR/bin/rsem/bin
    popd
    pinfo "rsem installation complete."    
}


######################################################
# kallisto
function kallisto_install {
    pinfo "Installing kallisto..."
    download_software kallisto
    tar xzvf $kallisto_FILE
    pushd kallisto_linux-v$kallisto_VERSION
    mkdir -p $IRAP_DIR/bin/kallisto/bin
    cp kallisto $IRAP_DIR/bin/kallisto/bin
    popd
    pinfo "kallisto installation complete."    
}

# salmon
function salmon_install {
    pinfo "Installing salmon..."
    download_software SALMON
    tar xzvf $SALMON_FILE
    pushd Salmon-${SALMON_VERSION}_linux_x86_64
    mkdir -p $IRAP_DIR/bin/salmon/bin $IRAP_DIR/bin/salmon/lib
    cp bin/* $IRAP_DIR/bin/salmon/bin
    cp lib/* $IRAP_DIR/bin/salmon/lib
    popd
    pinfo "salmon installation complete."    
}


######################################################
# Scripture
function scripture_install {
#download ftp://ftp.broadinstitute.org/pub/papers/lincRNA/scripture.jar
    pinfo "Installing scripture..."
    download_software SCRIPTURE
    mv scripture-$SCRIPTURE_VERSION.jar $IRAP_DIR/bin/scripture.jar
    cat <<EOF > $IRAP_DIR/scripts/scripture
#!$ENV_FP bash
java -Xmx8000m -jar \$IRAP_DIR/bin/scripture.jar \$*
EOF
    chmod +x $IRAP_DIR/scripts/scripture
    
# scripture requires IGVTools
# IGV 2.1 requires Java 6 or greater. 
    download_software IGV
    unzip $IGV_FILE
    cp IGV_$IGV_VERSION/* $IRAP_DIR/bin/
    chmod +x $IRAP_DIR/bin/igv*
    rm -rf IGV_$IGV_VERSION $IGV.zip
    
    cat <<EOF > $IRAP_DIR/bin/igv.sh
#!$ENV_FP bash
java -Dapple.laf.useScreenMenuBar=true -Xmx750m -jar $IRAP_DIR/bin/igv.jar $*
EOF
    
    IGVTOOLS=igvtools_nogenomes_$IGV_TOOLS_VERSION
    download_software IGV_TOOLS
    unzip $IGV_TOOLS_FILE
    cp IGVTools/* $IRAP_DIR/bin/
    echo "exit 0" >> $IRAP_DIR/bin/igvtools
    pinfo "Installing scripture...done"
}

function quant_install {
    cufflinks1_install
    cufflinks2_install
    htseq_install
    flux_capacitor_install
    scripture_install
    NURD_install
    stringtie_install
    rsem_install
    kallisto_install
    ##fusionmap_install
    salmon_install
    #isoem_install
    #sailfish_install
    #mmseq_install
    #ireckon_install
}

######################################################
# Fastq QC
function fastq_qc_install { 
    ###############
    # Fastx toolkit
    # Collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing.
    pinfo "Downloading and installing FASTX toolkit..."
    download_software FASTX
    tar xvjf $FASTX_FILE
    mv bin/* $BIN_DIR
    pinfo "Downloading and installing FASTX toolkit...done."
    
    ########
    # fastqc
    # Quality control tool for high throughput sequence data
    pinfo "Installing fastqc..."
    download_software FASTQC
    unzip $FASTQC_FILE
    rm -rf $IRAP_DIR/bin/FastQC  $IRAP_DIR/bin/fastqc
    mv FastQC $IRAP_DIR/bin
    pushd $IRAP_DIR/bin/FastQC
    chmod 755 fastqc
    sed "s.^#\!/usr/bin/perl.#\!$ENV_FP perl." -i fastqc
    ln -s `pwd`/fastqc $IRAP_DIR/bin
    popd
    pinfo "Installing fastqc...done."
}

######################################
function core_install {

    pinfo "*******************************************************"
    pinfo "IRAP_INSTALL_TOPLEVEL_DIRECTORY=$IRAP_DIR"
    pinfo "IRAP_SRC_TOPLEVEL_DIR=$SRC_DIR"
    pinfo "*******************************************************"

    gen_setup_irap $SETUP_FILE
    pinfo "Created setup file $SETUP_FILE"
    source $SETUP_FILE
    pinfo "Loaded $SETUP_FILE"

    #############
    #
    pinfo "Installing irap files..."
    if [ -h $IRAP_DIR/scripts ]; then
	pinfo "$IRAP_DIR/scripts is a symbolic link...skipping update."
    else
	if [ ! -e $IRAP_DIR/scripts ]; then
	    mkdir -p $IRAP_DIR/scripts
	fi
	cp -r $SRC_DIR/scripts/* $IRAP_DIR/scripts	
	# install should always be ran from the source 
	chmod -x $IRAP_DIR/scripts/irap_install.sh
	# update the env path
	if [ "$DEF_ENV" != "$ENV_FP" ]; then
	    sed -i "s|^#\!$DEF_ENV|#\!$ENV_FP|" $IRAP_DIR/scripts/*
	fi
    fi

    if [ -h $IRAP_DIR/aux ]; then
	pinfo "$IRAP_DIR/aux is a symbolic link...skipping update."
    else
	cp -r $SRC_DIR/aux $IRAP_DIR
    fi
    pinfo "Installing irap files...done."

    # examples
    cp -r $SRC_DIR/examples $IRAP_DIR

    #############
    # fastq utils
    # Fastq processing utilities
    pinfo "Compiling and installing fastq/bam processing programs..."
    pushd $SRC_DIR/src/fastq_utils
    make -j $J -B
    TARGETS=`grep "^TARGETS=" Makefile| cut -f 2 -d=`
    cp $TARGETS $BIN_DIR
    popd
    pushd $SRC_DIR/src/bamutils
    make -B
    cp bam_pe_insert bam_fix_NH bam_fix_se_flag bam_tophat2_fix bamRindex $BIN_DIR
    popd
    pinfo "Compiling and installing fastq/bam processing programs...done."

    # create the VERSION file
    echo $IRAP_VERSION > $IRAP_DIR/version
    pinfo "Core installation complete."
}
###############################################
# may be useful...but requires a web server
function jbrowse_install {
    pinfo "Installing jbrowse..."
    download_software JBROWSE
    rm -f  $IRAP_DIR/aux/jbrowse.zip
    mv $JBROWSE_FILE $IRAP_DIR/aux/jbrowse.zip
    unzip $IRAP_DIR/aux/jbrowse.zip
    # TODO: install deps
    pinfo "Uncompressing and installing jbrowse..."
    pushd JBrowse-*
    sed -i "s|-l extlib/| -l $IRAP_DIR|" setup.sh
    sed -i "s|bin/cpanm|cpanm|" setup.sh

    pinfo "Uncompressing and installing jbrowse...extra PERL packages"
    # TODO: do it only once
    perl_packages_jbrowse_install

    #cpan -f -i Module::CoreList  < /dev/null
    #cpan -f -i ExtUtils::MakeMaker < /dev/null
    #cpan -f -i Build < /dev/null
    #cpanm -v   GD
    #
    download_software SAMTOOLS
    tar xvjf $SAMTOOLS_FILE
    export SAMTOOLS="$PWD/samtools-${SAMTOOLS_VERSION}";    
    #svn export https://samtools.svn.sourceforge.net/svnroot/samtools/tags/samtools-0.1.7/;
    perl -i -pe 's/^CFLAGS=\s*/CFLAGS=-fPIC / unless /\b-fPIC\b/' samtools-${SAMTOOLS_VERSION}/Makefile;
    make -C samtools-${SAMTOOLS_VERSION} -j3 lib;
    ln -s samtools-${SAMTOOLS_VERSION} samtools
    # It is necessary to clean .cpan 
    # yes, weird!
    rm -rf $IRAP_DIR/.cpan.bak2
    if [ -e ~/.cpan ]; then
     	mv ~/.cpan $IRAP_DIR/.cpan.bak2
    fi
    # Files found in blib/arch: installing files in blib/lib into architecture dependent library tree
    #Installing /home/nf/perl5/lib/perl5/x86_64-linux/auto/DBI/DBI.so

    #unset INSTALL_BASE
    cpanm -l $IRAP_DIR --notest -f local::lib < /dev/null    
    set +e
    cpanm -v --notest -l $IRAP_DIR --installdeps . < /dev/null;
    cpanm -v --notest -l $IRAP_DIR --installdeps . < /dev/null;
    cpanm  -l $IRAP_DIR -f -v  --installdeps . < /dev/null;
    set -e
    if [ -e $IRAP_DIR/.cpan.bak2 ] ; then
	rm -rf ~/.cpan
	mv $IRAP_DIR/.cpan.bak2 ~/.cpan
    fi
    # cpanm -v --notest -l $IRAP_DIR --installdeps . < /dev/null;
    pinfo "Uncompressing and installing jbrowse...extra PERL packages (done)"
    pinfo "Uncompressing and installing jbrowse...compiling wig2png"
    pushd src/wig2png
    ./configure
    make    
    popd
    pinfo "Uncompressing and installing jbrowse...compiling wig2png (done)"
    ./setup.sh < /dev/null
    make 
    pinfo "Uncompressing and installing jbrowse...copying binaries and libs "
    cp -rf src $IRAP_DIR    
    #cp -rf extlib $IRAP_DIR
    # fix permissions
    chmod 755 bin/*
    cp bin/* $IRAP_DIR/bin 
    echo `pwd`   
    chmod 755 blib/script/*
    cp blib/script/* $IRAP_DIR/bin    
    #cp -rf blib $IRAP_DIR
    pinfo "Uncompressing and installing jbrowse...copying binaries and libs (done)"
    pinfo "Uncompressing and installing jbrowse... downloading and installing extra programs"
    # get some required 2utilities
    for f in $JBROWSE_EXTRA_UTILS; do
	pinfo "Uncompressing and installing jbrowse... downloading and installing extra programs ($f)"
	download $JBROWSE_EXTRA_UTILS_URL/$f
	chmod +x $f
	cp $f $IRAP_DIR/bin
	pinfo "Uncompressing and installing jbrowse... downloading and installing extra programs (done)"
    done
    pinfo "Uncompressing and installing jbrowse... downloading and installing extra programs (done)"    
    popd
    pinfo "jbrowse installation complete."
}
function new_jbrowse_install {
    pinfo "Installing new jbrowse..."
    download_software NEW_JBROWSE
    mv $NEW_JBROWSE_FILE $IRAP_DIR/aux/latest_jbrowse.zip
    unzip $IRAP_DIR/aux/latest_jbrowse.zip
    # 
    pinfo "Uncompressing and installing jbrowse..."
    pushd JBrowse-*
    sed -i "s|-l extlib/| -l $IRAP_DIR|" setup.sh
    sed -i "s|bin/cpanm|cpanm|" setup.sh

    pinfo "Uncompressing and installing jbrowse...extra PERL packages"
    perl_packages_jbrowse_install
    #
    download_software SAMTOOLS
    tar xvjf $SAMTOOLS_FILE
    export SAMTOOLS="$PWD/samtools-${SAMTOOLS_VERSION}";    
    #svn export https://samtools.svn.sourceforge.net/svnroot/samtools/tags/samtools-0.1.7/;
    perl -i -pe 's/^CFLAGS=\s*/CFLAGS=-fPIC / unless /\b-fPIC\b/' samtools-${SAMTOOLS_VERSION}/Makefile;
    make -C samtools-${SAMTOOLS_VERSION} -j3 lib;
    ln -s samtools-${SAMTOOLS_VERSION} samtools
    # It is necessary to clean .cpan 
    # yes, weird!
    rm -rf $IRAP_DIR/.cpan.bak2
    mv ~/.cpan $IRAP_DIR/.cpan.bak2
    #unset INSTALL_BASE
    cpanm -f -l $IRAP_DIR local::lib < /dev/null    
    set +e
    # cpanm -v --notest -l $IRAP_DIR --installdeps . < /dev/null;
    cpanm -f -v -l $IRAP_DIR --installdeps . < /dev/null;
    set -e
    rm -rf ~/.cpan
    mv $IRAP_DIR/.cpan.bak2 ~/.cpan 
    # cpanm -v --notest -l $IRAP_DIR --installdeps . < /dev/null;
    pinfo "Uncompressing and installing jbrowse...extra PERL packages (done)"
    pinfo "Uncompressing and installing jbrowse...compiling wig2png"
    pushd src/wig2png
    ./configure
    make    
    popd
    pinfo "Uncompressing and installing jbrowse...compiling wig2png (done)"
    ./setup.sh < /dev/null
    make 
    pinfo "Uncompressing and installing jbrowse...copying binaries and libs "
    cp -rf src $IRAP_DIR    
    #cp -rf extlib $IRAP_DIR
    # fix permissions
    chmod 755 bin/*
    cp bin/* $IRAP_DIR/bin 
    echo `pwd`   
    chmod 755 blib/script/*
    cp blib/script/* $IRAP_DIR/bin    
    #cp -rf blib $IRAP_DIR
    pinfo "Uncompressing and installing jbrowse...copying binaries and libs (done)"
    pinfo "Uncompressing and installing jbrowse... downloading and installing extra programs"
    # get some required utilities
    for f in $NEW_JBROWSE_EXTRA_UTILS; do
	pinfo "Uncompressing and installing jbrowse... downloading and installing extra programs ($f)"
	download $NEW_JBROWSE_EXTRA_UTILS_URL/$f
	chmod +x $f
	cp $f $IRAP_DIR/bin
	pinfo "Uncompressing and installing jbrowse... downloading and installing extra programs (done)"
    done
    pinfo "Uncompressing and installing jbrowse... downloading and installing extra programs (done)"    
    popd
    pinfo "jbrowse installation complete."
}

function data_install {
    pinfo "Creating data folder $IRAP_DIR/data..."
    mkdir -p $IRAP_DIR/data
    cp -r $SRC_DIR/data/contamination $IRAP_DIR/data
    mkdir -p $IRAP_DIR/data/reference
    mkdir -p $IRAP_DIR/data/raw_data
    pinfo "Creating data folder...done."
}

function ucsc_utils_install {

    pinfo "Installing extra programs from UCSC..."
    for f in $JBROWSE_EXTRA_UTILS; do
	pinfo "Downloading and installing $f..."
	download $JBROWSE_EXTRA_UTILS_URL/$f
	chmod +x $f
	cp $f $IRAP_DIR/bin
	pinfo "Downloading and installing $f...done"
    done
    pinfo "Installing extra programs from UCSC...done."
}
#############################################
# WIP
function ireckon_install {
    pinfo "Installing Ireckon...deprecated"
    exit
    download_software IReckon
    mkdir -p $IRAP_DIR/bin/ireckon
    cat <<EOF > $IRAP_DIR/bin/ireckon/ireckon
#!/bin/bash
if [ "$MEM-" = "-" ] ; then
MEM=15000M
fi
#Requirements:
# - iReckon works on linux. You need to have at least java-1.6 and the latest version of BWA installed and added to the PATH.
# - For large datasets and genomes, anticipate an important memory cost and running time (It is usually around 16G and 24 hours on 8 processors for human RNA-Seq with 60M read pairs).
#   The working directory (output directory) should have enough memory space (>100Go for the previous example). 
java -Xmx${MEM}M  -jar $IRAP_DIR/bin/IReckon-1.0.6.jar $* 
EOF
    chmod +x $IRAP_DIR/bin/ireckon
    mv IReckon-$IReckon_VERSION.jar $IRAP_DIR/bin/ireckon
    rm -f $IRAP_DIR/bin/ireckon/bwa
    ln -s  $IRAP_DIR/bin/bwa/bin/bwa  $IRAP_DIR/bin/ireckon/ 
    # deps
    download_software SAVANT
    chmod +x $SAVANT_FILE
    cat <<EOF > savant_responses
CreateDesktopShortcut: No
CreateQuickLaunchShortcut: Yes
InstallDir: $IRAP_DIR
InstallMode: Standard
InstallType: Typical
LaunchApplication: No
ProgramFolderName: Savant
SelectedComponents: Default Component
ViewReadme: No
EOF
    ./SAVANT_FILE -prefix $IRAP_DIR/bin --response-file savant_responses -mode silent
    #
    download http://genomesavant.com/savant/dist/v2_0_2/FormatTool.sh
    chmod +x FormatTool.sh
    mv FormatTool.sh $IRAP_DIR/bin/
    # requires BWA
    pinfo "Installing Ireckon...done."
}
######################################
#function soap_fusion_install {
#     MAPPER=soap_fusion
#      soap_fusion_VERSION=1.1
#      download ftp://public.genomics.org.cn/BGI/soap/Soapfusion/SOAPfusion_all_in_one_package.zip
#     unzip SOAPfusion_all_in_one_package.zip
#     pinfo "$MAPPER installation complete."    
#}

function NURD_install {
    pinfo "Installing NURD..."
    download_software NURD
    tar xzvf $NURD_FILE
#    pushd release_v$NURD_VERSION
    pushd NURD
    make
    cp NURD $IRAP_DIR/bin
    popd
    pinfo "NURD installation complete."    
}


######################################
function fusionmap_install {
# requires mono == osa

    pinfo "Installation of FusionMap..."
    download_software FUSIONMAP
    unzip $FUSIONMAP_FILE
    pushd FusionMap_$FUSIONMAP_VERSION
    cp bin/* $IRAP_DIR/bin
    popd
    pinfo "FusionMap installation complete."    
}
#FusionMap: detecting fusion genes from next-generation sequencing data at base-pair resolution Huanying Ge; Kejun Liu; Todd Juan; Fang Fang; Matthew Newman; Wolfgang Hoeck Bioinformatics (2011) 27 (14): 1922-1928. doi: 10.1093/bioinformatics/btr310

#      soap_fusion_VERSION=1.1
#      download ftp://public.genomics.org.cn/BGI/soap/Soapfusion/SOAPfusion_all_in_one_package.zip
#     unzip SOAPfusion_all_in_one_package.zip
#     pinfo "$MAPPER installation complete."    
#}

# python packages
function python_packages_install {
    pinfo "Installing python packages..."
    pinfo "Check python version... (2.6+ required)"
    min=$(python -c "import sys; print (sys.version_info[:])[1]")
    maj=$(python -c "import sys; print (sys.version_info[:])[0]")
    if [[ $maj == "2" ]] && [[ $min  -gt 5 ]] ; then
	pinfo "OK."
    else
	pinfo "You need Python2.6 or 2.7 to run this pipeline."
	exit 1
    fi
    ###########################################
    # only install pip if not already installed
    set +e
    PATH2PIP=`which pip 2> /dev/null`
    set -e
    CFLAGS_bak=$CFLAGS
    export CFLAGS=`echo $CFLAGS|sed -E "s|\-I[^ ]*boost||g"`
    if [ "$PATH2PIP-" == "-" ]; then
	wget --no-check-certificate https://bootstrap.pypa.io/get-pip.py
	python get-pip.py --user
	PATH2PIP=$IRAP_DIR/python/bin/pip
	pinfo "PIP installation complete."    
    else
	pinfo "pip already installed"
    fi
    $PATH2PIP install pysam==0.8.4 --user    
    export CFLAGS=$CFLAGS_bak
    pinfo "python packages installed"
}

function numpy_install {
    download http://sourceforge.net/projects/numpy/files/NumPy/1.8.1rc1/numpy-1.8.1rc1.tar.gz/download
    #download http://sourceforge.net/projects/numpy/files/latest/download
    tar xzvf numpy-1.8.1rc1.tar.gz
    pushd numpy-1.8.1rc1
    python setup.py build
    python setup.py install --user
    popd

}

function miso_install {
    pinfo "Installing miso..."
    # deps problems
    
    download http://sourceforge.net/projects/scipy/files/scipy/0.11.0/scipy-0.11.0.tar.gz
    tar xzvf scipy-0.11.0.tar.gz
    pushd scipy-0.11.0
    python setup.py build
    python setup.py install   --user
    popd

    download https://github.com/downloads/matplotlib/matplotlib/matplotlib-1.2.0.tar.gz
    tar xzvf matplotlib-1.2.0.tar.gz
    pushd matplotlib
    python setup.py build
    python setup.py install   --user
    popd

    download https://nodeload.github.com/yarden/MISO/legacy.zip/fastmiso
    unzip fastmiso  
    pushd yarden-MISO-*
    python setup.py build
    python setup.py install   --user
    pinfo "RSEM installation complete."    
}

PICARD_VERSION=1.119
PICARD_FILE=picard-tools-$PICARD_VERSION.zip
PICARD_URL=http://sourceforge.net/projects/picard/files/picard-tools/$PICARD_VERSION/$PICARD_FILE/download

function picard_install {
    pinfo "Installing Picard..."
    download_software PICARD
    unzip $PICARD_FILE
    mkdir -p $BIN_DIR/picard-tools
    mv picard-tools-$PICARD_VERSION/* $BIN_DIR/picard-tools
    pinfo "Picard installed"    
}

###############################
UPDATE_FILE_PERMS=n
INSTALL_JBROWSE=n
INSTALL_GCC=n
INSTALL_R3=n
INSTALL_BOOST=n
SPECIAL_SH_TO_USE=bash
OPTERR=0
while getopts "s:c:l:a:x:gmqpruhbdtfjvGKRB"  Option
do
    case $Option in
# update/reinstall
        a ) install=all;IRAP_DIR1=$OPTARG;;# send all output to a log file
	l ) install=minimal;IRAP_DIR1=$OPTARG;;# send all output to a log file
	b ) install=browser;IRAP_DIR1=$IRAP_DIR;;
	c ) install=core;IRAP_DIR1=$OPTARG;;  # run irap up to the given stage
	d ) USE_CACHE=n;install=download;IRAP_DIR1=$IRAP_DIR;; # download all the source packages
        f ) UPDATE_FILE_PERMS=y;;
	m ) install=mappers;IRAP_DIR1=$IRAP_DIR;;
	p ) install=p_pack;IRAP_DIR1=$IRAP_DIR;;
	q ) install=quant;IRAP_DIR1=$IRAP_DIR;;
	r ) install=r_pack;IRAP_DIR1=$IRAP_DIR;;
        s ) SRC_DIR=$OPTARG;;# send all output to a log file
	u ) install=core;IRAP_DIR1=$IRAP_DIR;; # update
	x ) install=software_install;IRAP_DIR1=$IRAP_DIR;SOFTWARE=$OPTARG;;
	t ) install=testing;IRAP_DIR1=$IRAP_DIR;;
	v ) install=collect_software_versions;IRAP_DIR1=$IRAP_DIR;;
	j ) INSTALL_JBROWSE=y;;
	G ) INSTALL_GCC=y;;
	R ) INSTALL_R3=y;;
	B ) INSTALL_BOOST=y;;
	K ) SPECIAL_SH_TO_USE=ksh;;
        h ) usage; exit;;
    esac
done


if [ "$IRAP_DIR1-" != "-" ]; then
    export IRAP_DIR=$IRAP_DIR1
fi


if [ "$IRAP_DIR-" = "-" ]; then
    echo ERROR: IRAP directory not defined. > /dev/stderr
    usage
    exit 1    
fi

if [ "$SRC_DIR-" = "-" ]; then
    usage
    exit 1    
fi

if [ "`uname 2>/dev/null`-" = "Linux-" ]; then
    OS=linux
else
    # dual world :)
    OS=mac
    pinfo " WARNING: This script will install binaries for Linux."
fi

# 
pinfo "iRAP $IRAP_VERSION"
# Check if env is available
DEF_ENV="/usr/bin/env"
ENV_FP=$DEF_ENV
if [ -x $ENV_FP ]; then
    pinfo "env found in $ENV_FP"
else
    ENV_FP="/bin/env"
    if [ -x $ENV_FP ]; then
	pinfo "env found in $ENV_FP"
    else
	echo "ERROR: env command not found - please ensure that it is in the PATH"
	exit 1
    fi
fi
# print system info
uname -a
# Check dependencies
check_dependencies
# Full path
pinfo "Checking paths..."
#IRAP_DIR=$(realpath -f "$IRAP_DIR")
IRAP_DIR=$(readlink -f "$IRAP_DIR")
SRC_DIR=$(readlink -f "$SRC_DIR")
pinfo "Checking paths...done."
#############################################
# print a few variables to help troubleshooting
pinfo "IRAP_DIR=$IRAP_DIR"
pinfo "SRC_DIR=$SRC_DIR"
pinfo "PATH=$PATH"
pinfo "CFLAGS=$CFLAGS"
pinfo "CXXFLAGS=$CXXFLAGS"
pinfo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
pinfo "R_LIBS=$R_LIBS"
pinfo "R_LIBS_USER=$R_LIBS_USER"

#############################################
#
BIN_DIR=$IRAP_DIR/bin
TMP_DIR=$IRAP_DIR/tmp
LIB_DIR=$IRAP_DIR/lib
INC_DIR=$IRAP_DIR/include

SETUP_FILE=$IRAP_DIR/irap_setup.sh

if [ "$CRAN_REPO-" == "-" ]; then
    #CRAN_REPO=http://cran.ma.imperial.ac.uk/
    CRAN_REPO=http://www.stats.bris.ac.uk/R/
fi
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIB_DIR
J=8
set -e

#################
# Setup directories
mkdir -p $TMP_DIR
mkdir -p $BIN_DIR
mkdir -p $LIB_DIR
mkdir -p $INC_DIR
mkdir -p $IRAP_DIR/scripts
# Python
export PYTHONUSERBASE=$IRAP_DIR/python
mkdir -p $PYTHONUSERBASE
mkdir -p $PYTHONUSERBASE/lib/python2.7/site-packages $PYTHONUSERBASE/lib/python2.6/site-packages
######################################
cd $TMP_DIR
# clean up before proceeding
pinfo "Cleaning up $TMP_DIR..."
rm -rf *
pinfo "Cleaning up $TMP_DIR...done."

# Keep a log file
mkdir -p $IRAP_DIR/install_logs
logfile=$IRAP_DIR/install_logs/`date "+%d%m%y%H%M"`.log
mkfifo ${logfile}.pipe
tee < ${logfile}.pipe $logfile &
exec &> ${logfile}.pipe
rm ${logfile}.pipe

# ensure that all files can be modified
if [ $UPDATE_FILE_PERMS == "y" ]; then
    pinfo "Fixing permissions..."
    chmod -R +w $IRAP_DIR
    pinfo "Fixing permissions...done."
fi

#if [  "$install" == "all"   &&  ! -e $IRAP_DIR  ]; then
#    echo "Directory $IRAP_DIR already exists. Please delete it before proceeding with the installation"
#    exit 1
#fi
if [ "$install" == "download" ]; then
    download2cache
    pinfo "Log saved to $logfile"
    exit 0
fi

if [ "$install" == "collect_software_versions" ]; then
    check_for_irap_env
    collect_software_versions
    pinfo "Log saved to $logfile"
    exit 0
fi

if [ "$install" == "software_install" ]; then
    # TODO check if software exists
    ${SOFTWARE}_install
    pinfo "Log saved to $logfile"
    exit 0
fi

if [ "$install" == "core" ]; then
    # TODO: check if IRAP is installed (otherwise it will fail)
    core_install
    pinfo "Log saved to $logfile"
    exit 0
fi

if [ "$install" == "mappers" ]; then
    check_for_irap_env
    mappers_install
    pinfo "Log saved to $logfile"
    exit 0
fi

if [ "$install" == "r_pack" ]; then
    check_for_irap_env
    R_packages_install
    pinfo "Log saved to $logfile"
    exit 0
fi

if [ "$install" == "p_pack" ]; then
    check_for_irap_env
    perl_packages_install
    pinfo "Log saved to $logfile"
    exit 0
fi

if [ "$install" == "quant" ]; then
    check_for_irap_env
    quant_install
    pinfo "Log saved to $logfile"
    exit 0
fi

if [ "$install" == "browser" ]; then
    check_for_irap_env
    jbrowse_install
    pinfo "Log saved to $logfile"
    exit 0
fi

if [ "$install" == "data" ]; then
    check_for_irap_env
    data_install
    pinfo "Log saved to $logfile"
    exit 0
fi

if [ "$install" == "testing" ]; then
    set -e
    check_for_irap_env
    isoem_install
    #sailfish_install
    #NURD_install
    rsem_install
    pinfo "Log saved to $logfile"
    exit 0
fi

#############
# all
deps_install $install
core_install
pinfo "Loading environment $SETUP_FILE..."
source $SETUP_FILE
# force cpan installation/reconfiguration
rm -f $IRAP_DIR/.cpan.irap.done
pinfo "PATH=$PATH"
pinfo "IRAP_DIR=$IRAP_DIR"
env |  grep IRAP_DIR
pinfo "Loading environment $SETUP_FILE...done."
#check_for_irap_env
#R_packages_install
R_packages_install

fastq_qc_install
perl_packages_install

if [ "$install" == "minimal" ]; then
    bowtie2_install
    bowtie1_install
    tophat2_install
    star_install
    htseq_install
    cufflinks2_install
   
   pinfo "WARNING: You chose to install the minimal installation of iRAP. Only the following tools will be available: bowtie1, bowtie2, tophat2, star, cufflinks2 "

else
    mappers_install
    quant_install
    # report
    if [ $INSTALL_JBROWSE == "y" ] ; then
	jbrowse_install
    fi
fi

collect_software_versions
# data directory
data_install
pinfo "Installation complete."
pinfo "Log saved to $logfile"
exit 0

