


## Column Nr	Name	Value	Description
## 1	Locus	chrom:start-end[W|C]	identifier of the transcriptional locus, given by the chromosome (chrom), start respectively end position, and the strand (Watson or Crick).
## 2	Transcript_ID	String	transcript identifier from the reference annotation.
## 3	Coding	[CDS|NC]	specifies whether the transcript has an annotated coding sequence (CDS) or not (NC)
## 4	Length	Integer	the mature length of the transcript after splicing out introns, disregarding the poly-A tail, as annotated in the reference annotation
## 5	Expressed Fraction	Float	fraction of RNA molecules that represent transcripts that are qualitatively equal to this RNA form
## 6	Expressed Number	Integer	absolute number of expressed RNA molecules
## 7	Library Fraction	Float	fraction of cDNA molecules in the final library that have been produced from this transcript
##  	Library Number	Integer	absolute number of cDNA fragments generated from this transcript
## 9	Sequenced Fraction	Float	fraction of total reads that have been sequenced from this transcript
## 10	Sequenced Number	Integer	absolute number of reads sequenced from this transcript
## 11	Covered Fraction	Float	fraction of the transcript that is covered by reads
## 12	Chi Square	Integer	chi-square goodness of fit measurement of coverage uniformity
## 13	Coefficient of Variation	Float	coefficient of variation for transcript coverage

# data in simulation
read.flux.transcript.pro <- function(filename) {
  # transcript profile
  transcript.profile <- read.csv(filename,sep="\t",header=F)
  colnames(transcript.profile) <- c("locus","tid","coding","length","expressed_frac","expressed_number","library_fraction",
                                    "lib_number",
                                    "sequenced_fraction",
                                    "sequenced_number",
                                    "covered_fraction",
                                    "chi_square",
                                    "coeff_var")
  return(transcript.profile)
}



merge.prof.gff3 <- function(transcript.profile,mrna.gff3) {

  #merge the annot with profile info
  m<-merge(transcript.profile,mrna.gff3,by.x="tid",by.y="ID",all.x=TRUE,sort=F)
  m
}

IRAP.DIR <- Sys.getenv(c("IRAP_DIR"))
source(paste(IRAP.DIR,"aux/R","irap_utils.R",sep="/"))
source(paste(IRAP.DIR,"aux/R","irap_misc.R",sep="/"))


mygff3 <- "/home/nf/Research/Projects/WIP/IRAP/irap/tests/test_files/simulation/Homo_sapiens.GRCh37.66.gff3"
trans.pro.file <- "/home/nf/Research/Projects/WIP/IRAP/irap/tests/test_files/simulation/test_se_100_5_3.1.pro"
counts.file <- "/home/nf/Research/Projects/WIP/IRAP/irap/tests/test_files/simulation/genes.raw.htseq2.tsv"

# gene id in parent
gff3 <- load.gff3(mygff3,"mRNA")


prof <- read.flux.transcript.pro(trans.pro.file)

m<-merge.prof.gff3(prof,gff3)

# now agregate the counts per gene
prof.by.gene <- aggregate(m[,c("sequenced_number","expressed_number")],by=list(as.character(m$Parent)),FUN=sum)

#

# Compare
counts <- read.table(counts.file,sep="\t",header=T)

all <- merge(prof.by.gene,counts,by.x="Group.1",by.y="Gene",all.x=TRUE)



# compare
v <- all[,c("wvL1","sequenced_number")]
rownames(v)<-all$Group.1
colnames(v) <- c(colnames(v)[1],"Ref")

# barplot
gen.plot(paste(counts.file,".barplot",sep=""),
         dir="",
         to.plot=function() {
           barplot(as.matrix(v),ylab="Reads")
         }
         )

# scatterplot
gen.plot(paste(counts.file,".scatter",sep=""),
         dir="",
         to.plot=function() {
           cor.scatterplot(v,colnames(v)[1],"Ref","Reads per gene")
         }
         )

# ma_plot?
min_count <- 1
row.sel <- apply(v,1,max)>min_count 
cds <- newCountDataSet(v[row.sel,],c('A','B'))
cds <- estimateSizeFactors(cds)
## estimates a dispersion value for each gene, then fits a curve
result <- try(cds <- estimateDispersions(cds));
if(class(result) == "try-error") {
  print("Parametric dispersion failed with default settings. Trying with local fit.");
  result <- try(cds <- estimateDispersions(cds,fitType="local"));
  if(class(result) == "try-error") {
    print("Parametric dispersion failed with default settings. Trying with method=blind.");
                                        # blind=estimate across conditions
    result <- try(cds <- estimateDispersions(cds,method="blind",fitType="local",sharingMode="fit-only"));
    if(class(result) == "try-error") {
                                        # pooled-CR: crossed factors
      print("Ooops, it seems that you need to manually tune DEseq  to estimate the dispersion.")
      q("no")
    }
  }
  print("OK")
}
print("Dispersion estimation complete.")

t<-counts(cds,normalized=TRUE)

gen.plot(paste(counts.file,".dispersion.png",sep=""),
         dir="",
         to.plot=function() {
           plotDispEsts(cds)
         }
         )

#v<-cbind(c(10,30,20),c(5,15,20))
#rownames(v) <- c('a','b','c')
# rank
o1<-order(v[,1])
v<-v[o1,]

o2<-order(v[,2])

order(v[,1])-o2
wilcox.test(v[,1],v[,2],paired=T)




# boxplot
# how many reads were aligned?
# how many reads were used in the quant?
