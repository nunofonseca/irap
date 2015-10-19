#!/bin/sh
export MERGE_EXTRA_OPTIONS=--basename_header
exec irap_single_lib2report $@
