isoDetector <-
function(countFile, count_col = 1, bedFile, fragSizeFile, readLen, output, lmax=500, 
        minCount=5, minObs=2, knownIsoforms=NULL, maxBreaks=5, 
        pvalBreaks=0.05, pvalExpress=0.01, foldExpress=1/5,
        eLenMin=1.0, muMin=0.01, pMaxRel=10, pMaxAbs=2000, 
        verbose=1)
{
  
  # --------------------------------------------------------- 
  # check intput files
  # --------------------------------------------------------- 
  
  if(file.access(countFile, mode=4) != 0){
    stop(sprintf("cannot read file %s.\n", countFile))
  }
  
  if(file.access(bedFile, mode=4) != 0){
    stop(sprintf("cannot read file %s.\n", bedFile))
  }
  
  if(file.access(fragSizeFile, mode=4) != 0){
    stop(sprintf("cannot read file %s.\n", fragSizeFile))
  }
  
  if(!is.null(knownIsoforms)){
    
    if(! grepl("RData$", knownIsoforms, perl=TRUE)){
      stop("knownIsoforms must be NULL or an .RData file\n")
    }
    
    if(file.access(knownIsoforms, mode=4) != 0){
      stop(sprintf("cannot read file %s.\n", knownIsoforms))
    }
  }
  
  if(!(count_col %in% c(1,2))){
    stop("count_col must be either 1 or 2 \n")
  }
  
  # --------------------------------------------------------- 
  # utility function to check wether it worth to dissect
  # RNA isoforms for a transcript cluster
  # --------------------------------------------------------- 

  checkIt <- function(gene, minCount, minObs){
    kpit  = 1
    w2use = which(gene$count$count >= minCount)
    if(length(w2use)<=minObs){ kpit = 0 }
    kpit
  }
  
  # --------------------------------------------------------- 
  # read in fragment size data 
  # --------------------------------------------------------- 
  
  md = read.table(fragSizeFile)
  
  if(ncol(md) != 2){
    stop(fragSizeFile, " should have 2 columns for Freq and Len\n")
  }
  
  names(md) = c("Freq", "Len")
  
  pd = rep(0, lmax)
  w2 = which(md$Len <= lmax)
  pd[md$Len[w2]] = md$Freq[w2]
  pdDist = pd/sum(pd)
  
  # --------------------------------------------------------- 
  # load counts data and summarize them gene by gene
  # --------------------------------------------------------- 
  
  geneMod = loadData(countFile, count_col, bedFile, pdDist, readLen, lmax)

  if(verbose >= 1){
    message("finish loading data\n")
  }
  
  # --------------------------------------------------------- 
  # select those genes that are non-trival for isoform study
  # --------------------------------------------------------- 
    
  keepIt  = sapply(geneMod, checkIt, minCount=minCount, minObs=minObs)
  w2kp    = which(keepIt==1)
  
  if(verbose >= 1){
    nn1   = length(geneMod)
    nn2   = length(w2kp)
    str1  = "non-trival for isoform study"
    str2  = sprintf("among %d genes, %d are %s\n", nn1, nn2, str1)
    message(str2)
  }
  
  geneMod = geneMod[w2kp]
  
  # --------------------------------------------------------- 
  # load information of known isoforms
  # --------------------------------------------------------- 
  
  if(!is.null(knownIsoforms)){
    load(knownIsoforms)
    if(! exists("isoAll")){
      stop("ah.. I expect an list of isoAll for knownIsoforms...\n")
    }
  }
  
  # --------------------------------------------------------
  # next, gene by gene, dissect isoform
  # --------------------------------------------------------
  
  nms = names(geneMod)
  sep = "\n--------------------------------------------------------\n"
  
  if(verbose >= 1){
    str1 = sprintf("isoform selection: 1, ..., %d", length(nms))
    message(sep, str1, sep)
  }
  
  for(i in 1:length(nms)){
    
    if(verbose >= 1){
      if(i %% 100 == 0){
        message(sep, i, "  ", date(), sep)
      }
    }
    
    nm1  = nms[i]
    ge1  = geneMod[[nm1]]
    
    if(!is.null(knownIsoforms)){
      isoforms = isoAll[[nm1]]
    }else{
      isoforms = NULL
    }
    
    gm1  = geneModel(ge1, d=readLen, pdDist, isoforms, lmax, 
                     maxBreaks, pvalBreaks, pvalExpress, foldExpress, 
                     isoformOnly=FALSE, runIsoSelect=TRUE, eLenMin,  
                     muMin, pMaxRel, pMaxAbs, verbose)
    geneMod[[nm1]] = gm1
  }
  
  save(geneMod, file=output)
  
}

