#!/bin/bash
set -e
# Collect and save the versions of the software installed
install_script=../scripts/irap_install.sh
DEST_FILE=$IRAP_DIR/aux/mk/irap_versions.mk

echo  "# Generated `date` " > $DEST_FILE
grep "VERSION=" $IRAP_DIR/scripts/irap_install.sh |sort -u | grep -v "^#"| sed -e "s/^ *//" >> $DEST_FILE
# R packages
if [ "`which R 2>/dev/null`-" = "-" ]; then
    echo "R not found" > /dev/stderr
else
    irap_R_package_version.R edgeR Rsamtools limma DEXSeq baySeq >> $DEST_FILE
fi

if [ -x $IRAP_DIR/scripts/R3 ]; then
    irap_R3_package_version.R DESeq DESeq2 piano  >> $DEST_FILE
else 
    echo "R3 not found" > /dev/stderr
fi
sed -i "s/VERSION=/version=/" $DEST_FILE
exit 0


