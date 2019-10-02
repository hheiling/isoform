knownIsoforms <-
function(nonOverlapExonFile){
  
  inf = read.table(nonOverlapExonFile, as.is=TRUE, sep="\t")
  
  dim(inf)
  inf[1:2,]
  
  if(ncol(inf) != 9){
    warning("I expect 9 columns.. hm. I will treat the last column as annotation\n")
  }else{
    names(inf) = c("chr", "source", "feature", "start", "end", 
                   "score", "strand", "frame", "anno")
  }
  
  # --------------------------------------------------------- 
  # obtain gene_id
  # --------------------------------------------------------- 

  reg1   = regexpr('gene_id\\s(\\S+);', inf$anno, perl=TRUE)
  
  if(any(reg1 == -1)){
    stop("cannot find 'gene_id' for certain records\n")
  }
  
  len1   = attributes(reg1)[[1]]
  nadd   = length(unlist(strsplit("gene_id", split=""))) + 1
  geneId = substr(inf$anno, reg1+nadd, reg1+len1-2)

  # --------------------------------------------------------- 
  # obtain transcript_id
  # --------------------------------------------------------- 

  reg1   = regexpr('transcript_id\\s(\\S+);', inf$anno, perl=TRUE)
  
  if(any(reg1 == -1)){
    stop("cannot find 'transcript_id' for certain records\n")
  }
  
  len1   = attributes(reg1)[[1]]
  nadd   = length(unlist(strsplit("transcript_id", split=""))) + 1
  tranId = substr(inf$anno, reg1+nadd, reg1+len1-2)

  # --------------------------------------------------------- 
  # obtain exon_id
  # --------------------------------------------------------- 

  reg1   = regexpr('exon_id\\s(\\S+);', inf$anno, perl=TRUE)
  
  if(any(reg1 == -1)){
    stop("cannot find 'exon_id' for certain records\n")
  }
  
  len1   = attributes(reg1)[[1]]
  nadd   = length(unlist(strsplit("exon_id", split=""))) + 1
  exonId = substr(inf$anno, reg1+nadd, reg1+len1-2)

  # --------------------------------------------------------- 
  # obtain clustId
  # --------------------------------------------------------- 
  
  reg1   = regexpr('clustId\\s(\\S+);', inf$anno, perl=TRUE)
  
  if(any(reg1 == -1)){
    stop("cannot find 'clustId' for certain records\n")
  }
  
  len1   = attributes(reg1)[[1]]
  nadd   = length(unlist(strsplit("clustId", split=""))) + 1
  clustId = substr(inf$anno, reg1+nadd, reg1+len1-2)
  
  # --------------------------------------------------------- 
  # for each clustId collect the isoform information
  # --------------------------------------------------------- 
  
  info = data.frame(geneId=geneId, tranId=tranId, exonId=exonId, 
                    clustId=clustId, stringsAsFactors=FALSE)
  dim(info)
  info[1:2,]
    
  infoG = split(info[,1:3], info$clustId)
  length(infoG)
  infoG[[1]]
  
  isoforms = list()
  nms      = names(infoG)
  
  for(i in 1:length(infoG)){
    if(i %% 1000 == 0){ cat(i, date(), "\n") }
    
    u1    = nms[i]
    info1 = infoG[[u1]]
    tId1  = paste(info1$tranId, collapse=":")
    tId1  = unique(unlist(strsplit(tId1, split=":")))
    
    if(any(info1$exonId != 1:nrow(info1))){
      stop("hm.. some exon Ids do not match\n")
    }
    
    iso1  = matrix(0, nrow=nrow(info1), ncol=length(tId1))
    
    for(j in 1:length(tId1)){
      t1       = tId1[j]
      mapj     = grep(t1, info1$tranId)
      iso1[mapj,j] = 1
    }
    
    colnames(iso1) = tId1
    
    iso1 = unique(iso1, MARGIN=2)
    isoforms[[u1]] = iso1
  }
  
  
  isoforms
}



