#####################################################################
# default options for each protocol
# obviously, may be overriden by the users
#####################################################################

# Drop-seq
ifeq ($(sc_protocol),drop-seq)
mapper?=kallisto
quant_method?=umi_count
def_umi_read=index1
def_umi_offset=12
def_umi_size=8
def_cell_read=index1
def_cell_offset=0
def_cell_size=12
endif

###############################
# 10X
#
#bam_umi_count_params?=--min_reads 0 --multi_mapped
# --cb_cutoff Number of counts to filter cellular barcodes.
#umis_params?=--cb_cutoff auto

define set_10x_params=
mapper?=kallisto
umi_quant?=umi_count
sc_non_zero_rows?=1
cell_filt_min_features?=0.3
cell_filt_max_ERCC?=0.75
cell_filt_outliers?=y
qc?=report
endef


# theoretical locations of the barcodes
ifeq ($(sc_protocol),10x_v1)
# index1=read2/I7, index2=I5
$(eval $(call set_10x_params))
def_umi_read=index1
def_umi_offset=14
def_umi_size=10
def_cell_read=index1
def_cell_offset=0
def_cell_size=14
def_sample_read=index2
def_sample_offset=0
def_sample_size=8
## not 100% sure that the file is correct
def_known_umi_file=$(IRAP_DIR)/data/10x/737K-august-2016.txt.gz
endif

# files obtained from the 10x website
ifeq ($(sc_protocol),10x_v1p)
# index1=*I1*
# index2=*I2*
# index3=read_RA*_2
$(info set 10x_v1p)
$(eval $(call set_10x_params))
def_umi_read=index3
def_umi_offset=0
def_umi_size=10

def_cell_read=index1
def_cell_offset=0
def_cell_size=14

## sample
def_sample_read=index2
def_sample_offset=0
def_sample_size=8

## not 100% sure that the file is the correct one for 10xv1
def_known_umi_file=
def_known_cells_file=$(IRAP_DIR)/data/10x/737K-april-2016.txt.gz

$(info set 10x_v1p)
endif


ifeq ($(sc_protocol),10x_v2)
$(eval $(call set_10x_params))

# default values?
# index1=R1 file index2=I1 index1=I5 reads=R2
def_umi_read=index1
def_umi_offset=16
def_umi_size=10

def_cell_read=index1
def_cell_offset=0
def_cell_size=16

## sample
def_sample_read=index2
def_sample_offset=0
def_sample_size=8
endif



# smart-seq
ifeq ($(sc_protocol),smart-seq)
mapper?=none
quant_method?=kallisto
#trim_poly_at?=y
#trim_poly_at_len=8
qc?=on
endif
