#!/usr/bin/env bash

if [ "$1-" = "-" ] ; then
    echo "atlas_clean_html.sh HTML file" > /dev/stderr
    exit 1
fi

if [ ! -e "$1" ] ; then
    echo "ERROR: file not found $1" > /dev/stderr
    exit 1
fi

# remove header
LINE=`grep -n -i -F "<body" $1 | cut -f 1 -d:`
if [ $? != 0 ] || [ "$LINE-" = "-" ]; then
    echo "ERROR: Unable to find body in $1"
    exit 1
fi
HLINE=`expr $LINE + 1`
#<div id="footer">
FLINE=`grep -n -i -F "<div id=\"footer" $1 | cut -f 1 -d:`
if [ $? != 0 ]|| [ "$FLINE-" = "-" ]; then
    FLINE=`grep -n -i -F "</body>" $1 | cut -f 1 -d:`
    if [ $? != 0 ]|| [ "$FLINE-" = "-" ]; then
	echo "ERROR: Unable to find footer div or /body in $1"
	exit 1
    fi
fi
# footer
FLINE=`expr $FLINE - 1`
set -e
head -n +$FLINE $1 | tail -n +$HLINE | grep -v -F "<nav>  <ul>   <li><a href='info.html' >Info</a>" > $1.tmp
mv $1.tmp $1
exit 0
