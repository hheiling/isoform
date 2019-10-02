
library(isoform)

# -------------------------------------------------------------------------
# read genotype data and info
# -------------------------------------------------------------------------

setwd("~/research/eQTL_seq/CEU/data/")

cts = read.table("counts_per_cluster.txt", sep="\t",
header=TRUE, as.is=TRUE)
dim(cts)
cts[1:2,]

pcts  = cts[,-(1:2)]
rdall = colSums(pcts)
rd75  = apply(pcts, 2, quantile, probs=0.75)

cor(rdall, rd75)

rdth = data.frame(sampleID = names(pcts), rd = rd75,
  row.names = NULL, stringsAsFactors=FALSE)
dim(rdth)
head(rdth)

ginfo = read.table("CEU_genotype_1000bp+genebody_info.txt", sep="\t", 
header=TRUE, as.is=TRUE)
dim(ginfo)
ginfo[1:2,]

# number of SNPs pre cluster
table(table(ginfo$names))

gdata = read.table("CEU_genotype_1000bp+genebody_data.txt", sep="\t", 
header=TRUE, as.is=TRUE)
dim(gdata)
gdata[1:2,]

# ----------------------------------------------------------------------------
# prepare isoDu, gene models, only keep the first 10 samples
# ----------------------------------------------------------------------------

g2tests = c("chr1_555", "chr19_490")

setwd("~/research/eQTL_seq/CEU/tophat/")

Routs = list.files(pattern="_geneModel_knownIsoforms.RData")
length(Routs)

tags  = gsub("_geneModel_knownIsoforms.RData", "", Routs)

all(tags == rdth$sampleID)

gL  = list()

for(i in 1:10){
  load(Routs[i])
  gL[[tags[i]]] = geneMod[g2tests]
}

# ----------------------------------------------------------------------------
# prepare isoDu, fragment size distribution, genotype data
# ----------------------------------------------------------------------------

fragSizeFiles = paste(tags, "_insertSize_dist.txt", sep="")
pdDistL = readFragmentSizes(fragSizeFiles, lmax=300)

g2test = "chr19_490"

xData = t(gdata[which(ginfo$names == g2test),,drop=FALSE])
dim(xData)
colnames(xData) = paste("X", 1:ncol(xData), sep="")

# ---------------------------------------------------------
# run isoDu
# ---------------------------------------------------------

isoDuEx2 = list(readDepth=rdth[1:10,], xData=xData[1:10,],
  pdDistL=pdDistL[1:10], gL=gL)

save(isoDuEx2, file="../tophat_example/isoDuEx2.RData")


