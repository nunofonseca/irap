#!/bin/bash
# update iRAP's version
VERSION=$1
FILES_CHANGED=$2
if [ "$VERSION-" = "-" ]; then
    echo "ERROR: new_version.sh VERSION OFILE"
    exit 1
fi
set -e
FILESCHANGED=
F1="`ack -l VERSION src/*/*.c`"
for f in $F1; do
    echo "Updating $f"
    sed -i "s/#define VERSION .*/#define VERSION \"$VERSION\"/" $f
done

F2="docker/install_tests/Makefile docker/Makefile"
for f in $F2; do
    echo "Updating $f"
    sed -i "s/#define VERSION .*/#define VERSION \"$VERSION\"/" $f
done

F2a="`find docker -name "*.docker"` `find docker -name "Dockerfile"`"
for f in $F2a; do
    echo "Updating $f"
    sed -i "s/ENV IRAP_VERSION=.*/ENV IRAP_VERSION=$VERSION/" $f
done


#F3="scripts/irap_single_lib2report  scripts/irap_fastq_info aux/mk/irap_core.mk"
F3=" scripts/irap_fastq_info"
for f in $F3; do
    echo "Updating $f"
    sed -i "s/^version=.*/version=$VERSION/" $f
done
F4=
#aux/R/irap_utils.R
#sed -i "s/^irap_version<-.*/irap_version<-\"$VERSION\"/" aux/R/irap_utils.R

F5=
#F5=scripts/irap_install.sh
#sed -i "s/^IRAP_VERSION=.*/IRAP_VERSION=$VERSION/" $F5

F6=
#scripts/irap_single_lib
#sed -i "s/^VERSION=.*/VERSION=$VERSION/" $F6

echo $VERSION > version

echo $F1 $F2 $F2a $F3 $F4 $F5 $F6 version > $FILES_CHANGED

exit 0
