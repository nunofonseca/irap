#; -*- mode: Makefile;-*-
# =========================================================
# Copyright 2012-2018,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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
# =========================================================
# Rules for producing a simple status of the execution of the pipeline

#1-files
define stage1_status=
	irap_stage_stats.pl stage1 $(name) $(qual_filtering) "x" "x" "x" $(1)
endef
define stage2_status=
	irap_stage_stats.pl stage2 $(name) $(qual_filtering) $(mapper) "x" "x" $(1)
endef
define stage3_status=
	irap_stage_stats.pl stage3 $(name) $(qual_filtering) $(mapper) $(quant_method) "x"  $(1)
endef
define stage4_status=
	irap_stage_stats.pl stage4 $(name) $(qual_filtering) $(mapper) $(quant_method) $(de_method) $(1)
endef
define stage5_status=
	irap_stage_stats.pl stage5 $(name) $(qual_filtering) $(mapper) $(quant_method) $(de_method) $(1)
endef

##################################################################

$(report_toplevel_folder)/status.html: $(name)/$(name).status.tsv
	irap_report_status.R $(conf) $< $@.tmp && mv $@.tmp.html $@

$(name)/$(name).status.tsv: $(qc_toplevel_folder)/stage1.status $(mapper_toplevel_folder)/stage2.status $(quant_toplevel_folder)/stage3.status $(de_toplevel_folder)/stage4.status
	@cat $^ > $@

$(qc_toplevel_folder)/stage1.status: FORCE
	@$(call stage1_status,$(STAGE1_OUT_FILES)) > $@.tmp && mv $@.tmp $@

$(mapper_toplevel_folder)/stage2.status: FORCE
	@mkdir -p $(@D)
	@$(call stage2_status,$(STAGE2_OUT_FILES)) > $@.tmp && mv $@.tmp $@

$(quant_toplevel_folder)/stage3.status:  FORCE
	@mkdir -p $(@D)/
	$(info $(STAGE3_OUT_FILES))
	$(call stage3_status,$(STAGE3_OUT_FILES)) > $@.tmp && mv $@.tmp $@

#
ifeq (,$(STAGE4_OUT_FILES))
$(de_toplevel_folder)/stage4.status:
	touch $@
else
$(de_toplevel_folder)/stage4.status: FORCE
	@mkdir -p $(@D)
	@$(call stage4_status,$(STAGE4_OUT_FILES)) > $@.tmp && mv $@.tmp $@
endif

ifeq (,$(STAGE5_OUT_FILES))
$(de_toplevel_folder)/stage5.status:
	touch $@
else
$(de_toplevel_folder)/stage5.status: FORCE
	@mkdir -p $(@D)
	@$(call stage5_status,$(STAGE5_OUT_FILES)) > $@.tmp && mv $@.tmp $@
endif


FORCE:
