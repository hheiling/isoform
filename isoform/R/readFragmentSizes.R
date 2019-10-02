readFragmentSizes <-
function(fragSizeFiles, lmax=500){
  pdDistL = list()

  for(i in 1:length(fragSizeFiles)){
    
    fragSizeFile = fragSizeFiles[i]
    md = read.table(fragSizeFile)
    
    if(ncol(md) != 2){
      stop(fragSizeFile, " should have 2 columns for Freq and Len\n")
    }
    
    names(md) = c("Freq", "Len")
    
    pd = rep(0, lmax)
    w2 = which(md$Len <= lmax)
    pd[md$Len[w2]] = md$Freq[w2]
    pdDistL[[i]]   = pd/sum(pd)
  }
  
  pdDistL
}
