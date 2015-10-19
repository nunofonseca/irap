#!/bin/sh
export MERGE_EXTRA_OPTIONS="--basename_header --exclude_irap_suffix"
exec irap_single_lib2report $@
