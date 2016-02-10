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
#    $Id$
# =========================================================

# Fusion search

#
#fusion_valid_methods=tophat2 defuse fusionmap soapfusion 
fusion_valid_methods=fusionmap
fusion_valid_modes=single_method
# ensemble
#fusion_valid_modes=single_method ensemble

fusion_mode=single_method
fusion_method=none
#tophat2 defuse fusionmap 

fusion_support_reads=3
fusion_support_pairs=3

ifeq ($(fusion_method),none)
fusion_search=n
else
fusion_search=y
endif
FUSION_TARGETS=

ifeq ($(fusion_search),y)

# to
ifeq (tophat2,$(filter tophat2,$(fusion_method)))

# tophat2 - enable fusion search
# include non human
ishuman=n
ifneq (,$(findstring homo,$(species)))
ishuman=y
endif
ifneq (,$(findstring human,$(species)))
ishuman=y
endif

tophat2_map_options+= --fusion-search --fusion-min-dist 1000000 -num-fusion-reads $(fusion_support_reads) --num-fusion-pairs $(fusion_support_pairs)
ifeq ($(ishuman),y)
tophat2_map_options+= --non-human
endif


# rule to get the post processed tophat-fusion output
#tophat-fusion-post (to filter out fusion candidates)
#    result.html. A list of fusion candidates is given in HTML format. A sample list is found here.
#    result.txt. A text version of result.html.
# $(data)/$(mapper)/fusion/$(fusion_method)/lib_path/
# $(data)/$(mapper)/fusion/summary.tsv

endif


# PE only
ifeq (defuse,$(filter defuse,$(fusion_method)))

endif

########################################################
# FusionMap
# RL > 75bp
ifeq (fusionmap,$(filter fusionmap,$(fusion_method)))

FUSIONMAP_BASEDIR=$(name)/data/fusionmap

# Setup
SETUP_DATA_FILES+=$(FUSIONMAP_BASEDIR)/fusionmap.index $(FUSIONMAP_BASEDIR)/fusionmap.gm

# 
$(name)/data/fusionmap/fusionmap.index: $(reference_abspath)
	mkdir -p $(FUSIONMAP_BASEDIR)
	mono $(IRAP_DIR)/bin/FusionMap.exe --buildref $(FUSIONMAP_BASEDIR) $< fusionmap_index  && touch $@

$(name)/data/fusionmap/fusionmap.gm: $(gtf_file_abspath)
	mkdir -p $(FUSIONMAP_BASEDIR)
	mono $(IRAP_DIR)/bin/FusionMap.exe --buildgm $(FUSIONMAP_BASEDIR) $< fusionmap_index  fusionmap_refgene && touch $@

# Run FusionMap
# lib
# bamfile
# se|pe
# outfile
# index and gm need to be absolute paths?
define run_fusionmap=
	mkdir -p $(dir $(4))
	irap_gen_fusionmap_conf.sh $(abspath $(2)) $(max_threads) $(3) $(dir $(4))/$(notdir $(1)) > $(dir $(4))/$(notdir $(1)).conf && \
	mono $(IRAP_DIR)/bin/FusionMap.exe --semap $(abspath $(FUSIONMAP_BASEDIR))  fusionmap_index fusionmap_refgene  $(dir $(4))/$(notdir $(1)).conf && \
	cp $(dir $(4))/$(notdir $(1)).FusionReport.txt $(4)
endef

#
# Use the BAMs
$(name)/$(mapper)/fusionmap/%.fusion.tsv: $(name)/$(mapper)/%.se.hits.bam $(name)/data/fusionmap/fusionmap.index $(name)/data/fusionmap/fusionmap.gm
	$(call run_fusionmap,$*,$<,se,$@.tmp) && mv $@.tmp $@

$(name)/$(mapper)/fusionmap/%.fusion.tsv: $(name)/$(mapper)/%.pe.hits.bam $(name)/data/fusionmap/fusionmap.index $(name)/data/fusionmap/fusionmap.gm
	$(call run_fusionmap,$*,$<,pe,$@.tmp)  && mv $@.tmp $@

$(name)/$(mapper)/fusionmap/%.fusion.sum.tsv: $(name)/$(mapper)/fusionmap/%.fusion.tsv
	irap_Fusion_fm2descr --tsv "$^" --gtf $(gtf_file_abspath) -c $(max_threads)  -o $@.tmp && mv $@.tmp $@

############################
# Add to stage3 output files
FUSION_LIB_TARGETS=$(foreach p,$(pe),$(call lib2fusion_folder,$(p))$(p).fusion.tsv) $(foreach s,$(se),$(call lib2fusion_folder,$(s))$(s).fusion.tsv)


STAGE3_S_TARGETS+=$(FUSION_LIB_TARGETS) 

ifdef sop
ifeq ($(sop),pawg3_th2_mapping)
STAGE3_S_TARGETS+=$(foreach p,$(pe),$(call lib2fusion_folder,$(p))$(p).fusion.sum.tsv) $(foreach s,$(se),$(call lib2fusion_folder,$(s))$(s).fusion.sum.tsv)
FUSION_TARGETS+=$(foreach p,$(pe),$(call lib2fusion_folder,$(p))$(p).fusion.sum.tsv) $(foreach s,$(se),$(call lib2fusion_folder,$(s))$(s).fusion.sum.tsv)
endif
endif

STAGE3_OUT_FILES+=$(name)/$(mapper)/fusionmap/fusionmap_readcounts.tsv $(name)/$(mapper)/fusionmap/fusionmap_fusions.tsv
FUSION_TARGETS+=$(FUSION_LIB_TARGETS)

# counts file
$(name)/$(mapper)/fusionmap/fusionmap_readcounts.tsv:  $(FUSION_LIB_TARGETS)
	$(call pass_args_stdin,irap_Fusion_fm2tsv,$@.tmp, --tsv "$^" -o $@.tmp) && mv $@.tmp $@

$(name)/$(mapper)/fusionmap/fusionmap_fusions.tsv:  $(FUSION_LIB_TARGETS)
	$(call pass_args_stdin,irap_Fusion_fm2descr,$@.tmp, --tsv "$^"  -c $(max_threads) --gtf $(gtf_file_abspath) -o $@.tmp) && mv $@.tmp $@


endif

ifeq (soapfusion,$(filter soapfusion,$(fusion_method)))

endif

else


endif

do_fusion: $(FUSION_TARGETS)
