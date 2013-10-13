#!/bin/bash
set -e
# Collect and save the versions of the software installed
install_script=../scripts/irap_install.sh
DEST_FILE=$IRAP_DIR/aux/mk/irap_versions.mk

echo  "# Generated `date` " > $DEST_FILE
grep "VERSION=" $IRAP_DIR/scripts/irap_install.sh |sort -u | grep -v "^#"| sed -e "s/^ *//" >> $DEST_FILE
# R packages
irap_R_package_version.R DESeq edgeR  Rsamtools DEXSeq baySeq  >> $DEST_FILE

sed -i "s/VERSION=/version=/" $DEST_FILE
exit 0


