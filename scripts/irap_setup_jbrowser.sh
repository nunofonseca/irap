#!/bin/sh
#
# setup a browser with the results
set -e 

DIR=$1
#NAME=$2

if [ "$DIR-" = "-" ]; then
    echo "irap_setup_jbrowser.sh dir" >&2
    exit 1
fi

############################
#
mkdir -p $DIR
pushd $DIR > /dev/null

rm -rf JBrowse-*
unzip -q  -u -o $IRAP_DIR/aux/jbrowse.zip 
#
JBDIRS="bin docs img src"
rm -rf $JBDIRS
JBROWSE_DIR=`ls -d JBrowse-*` 
if [ ! -e $JBROWSE_DIR ]; then
    echo "Error: $JBROWSE_DIR not found"  >&2
    exit 1
fi
pushd $JBROWSE_DIR > /dev/null
# clean up before moving the new folders
rm -rf ../plugins sample_data
mv * ..
popd >/dev/null
# cleanup
rm -rf sample_data docs setup.sh Makefile.PL INSTALL 
#rm -rf bin
mkdir -p data

# 
mkdir -p data/raw/bam
mkdir -p data/raw/wig

#############
# custumize the main page
sed -i "s|<title>JBrowse</title>|<title>IRAP report - browser</title>|" index.html

# add link to main report
menu=menu.html
TTOPLEVEL=../

if [ ! -e "$TTOPLEVEL$menu" ] ; then
    echo "Error: file $TTOPLEVEL$menu not found"  >&2
    exit 1
fi

inject_menu=`cat $TTOPLEVEL$menu| tr "\n" " "`
CSS_FILE='<link rel="stylesheet" href="TTOPLEVELmenu.css"  type=text/css>'
sed  -E -i "s|(<body[^>]*)>.*|\1> $inject_menu|;s|TTOPLEVEL|$TTOPLEVEL|g;" index.html
sed  -E -i "s|<head>.*|<head>$CSS_FILE|" index.html


#sed -i "s|<body>.*|<body><div><H1 style='text-align:center'>IRAP - $2</H1><h2><a href='..'>Return to main report</a></h2></div>|" index.html

# 
cat <<EOF > jbrowse_conf.json
// JBrowse JSON-format configuration file
{

    // // uncomment and edit the example below to configure a faceted track selector
    trackSelector: {
        type: 'Faceted',
	displayColumns: ['key','Stage','Observation']
    },
    trackMetadata: {
       sources: [
           { type: 'csv', url: 'data/trackMetadata.csv' }
    ]
    },

    // the variable below does nothing
    placeholder: 1
}
EOF
popd > /dev/null
exit 0
