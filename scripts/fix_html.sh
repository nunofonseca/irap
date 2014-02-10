#!/bin/bash

# adds the top level toolbar to all html pages
# makes some extra fixes on the pages
menu=menu.html
# include /
toplevel_path=
inject_menu=`cat $menu| tr "\n" " "`

function fix_html_files {
    path=$1

    CSS_FILE='<link rel="stylesheet" href="TTOPLEVELmenu.css">'

    for f in `find . -maxdepth 1 -name  "*.html"`; do
	if [ `basename $f` != $menu ]; then
	    echo "Fixing $f..."
	    # add CSS FILE
	    #echo adding menu.css
	    sed  -E -i "s|<head>.*|<head>$CSS_FILE|" $f
	    # add menu
	    # replace toplevel path
	    # remove r2html
	    # add
	    sed  -E -i "s|(<body[^>]*)>.*|\1> $inject_menu|;s|TTOPLEVEL|$path|g;s|R2HTML||;s|aehts.css|irap.css|;s|href=\"*.*irap.css\"|href=${path}irap.css|" $f
#	    sed  -E -i "s|(<body[^>]*)>.*|\1> $inject_menu|;s|TTOPLEVEL|$path|g;s|R2HTML||;s|aehts.css|irap.css|;s|href=../../ |href=${path}irap.css |" $f
	    # remove the line (temporary)
	    sed -i "s|<p class='character'>Project: <a href='../index.html'>go to main page</a></p>||" $f
	    # keep this temporarily
	    sed -i "s|<object id=\"mathplayer\" .*||" $f
	    echo "'Fixing' $f...done."
	fi    
    done
    #echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>$$ $path<<<<<<<<<<<<<<<<<<<<<<<<<"
    #bash does not have a stack when invoking a function recursively
    if [ "$path-" != "../../../../../-" ]; then
	echo $path > .path
	for d in `find -L .  -maxdepth 1  -type d -name  "*"`; do	
	    if [ "$d" != "." ] ; then
		pushd $d > /dev/null
		if [ "$d" != "./data" ]; then
		    echo "Going to dir $d"
		    fix_html_files $path../
		fi
		popd > /dev/null
		path=`cat .path`
		#echo ">>>>>>>>>>>>>>>>>>>>>>>>>>$path $PWD"
	    fi
	done
    fi
}

#
fix_html_files
exit 0




