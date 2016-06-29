#!/bin/sh
# wrapper to irap_single_lib.sh
set -e
if which $0 >/dev/null 2>/dev/null; then
 source `dirname $0`/../irap_setup.sh
fi
irap_single_lib $@
exit 0
