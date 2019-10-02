
library(isoform)

data(isoDuEx2)

names(isoDuEx2$gL[[1]])

g2test     = "chr19_490"
tags       = names(isoDuEx2$gL)
xData      = isoDuEx2$xData[,1:2]
readDepth  = isoDuEx2$readDepth$rd
pdDistL    = isoDuEx2$pdDistL
outputFile = sprintf("eQTL_%s_known_iso_duOnly.txt", g2test)

all(isoDuEx2$readDepth$sampleID == tags)

isoDu(tags, isoDuEx2$gL, xData, readDepth, outputFile, pdDistL,
g2test, readLen=37, method="permutation", lmax=300,
nResample=10, maxTime=1000, duOnly=TRUE)


