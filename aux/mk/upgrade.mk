#; -*- mode: Makefile;-*-
# =========================================================
# Copyright 2012-2016,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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
#    $Id: 0.1.3 Nuno Fonseca Fri Dec 21 11:56:56 2012$
# =========================================================
# Rules for upgrading irap (software and files)



#
define main_version=
$(shell echo $(1)|sed 's/[dp].*$$//')
endef
# clone
irap_clone_dir=$(clone_dir)
irap_clone_version=$(call main_version,$(shell cat $(clone_dir)/version))
# current version
version=$(shell cat $(irap_path)/../version)
irap_cur_version=$(call main_version,$(version))


$(info $(irap_cur_version))
$(info $(irap_clone_version))
irap_clone_version=0.5.1
#
ifeq ($(irap_cur_version),$(irap_clone_version))
# just apply any iRAP patches
upgrade_irap: 
	$(irap_clone_dir)/scripts/irap_install.sh -s $(irap_clone_dir) -u
	$(irap_clone_dir)/scripts/irap_install.sh -s $(irap_clone_dir) -v
else
upgrade_irap: upgrade_irap_$(irap_cur_version)_to_$(irap_clone_version)
endif

upgrade_irap_0.5.0_to_0.5.1: 
	$(irap_clone_dir)/scripts/irap_install.sh -s $(irap_clone_dir) -x R3_packages
	$(irap_clone_dir)/scripts/irap_install.sh -s $(irap_clone_dir) -u
	$(irap_clone_dir)/scripts/irap_install.sh -s $(irap_clone_dir) -v


