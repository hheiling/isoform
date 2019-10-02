compareAnno.isoDetector <-
function(geneMod, gtf){

  # --------------------------------------------------------- 
  # read in gene annotation
  # --------------------------------------------------------- 
  
  info = read.table(gtf, sep="\t", as.is=TRUE)
  dim(info)
  info[1:2,]
  
  names(info) = c("chr", "source", "feature", "start", "end", 
                  "score", "strand", "frame", "anno")
  
  table(info$chr)
  table(info$source)
  table(info$feature)
  
  info   = info[info$feature == "exon",]
  dim(info)
  
  # --------------------------------------------------------- 
  # obtain gene_id
  # --------------------------------------------------------- 
  
  reg1   = regexpr('gene_id\\s(\\S+);', info$anno, perl=TRUE)
  len1   = attributes(reg1)[[1]]
  nadd   = length(unlist(strsplit("gene_id", split=""))) + 1
  info$geneId = substr(info$anno, reg1+nadd, reg1+len1-2)
  
  # --------------------------------------------------------- 
  # obtain transcript_id
  # --------------------------------------------------------- 
  
  reg1   = regexpr('transcript_id\\s(\\S+);', info$anno, perl=TRUE)
  len1   = attributes(reg1)[[1]]
  nadd   = length(unlist(strsplit("transcript_id", split=""))) + 1
  info$tranId = substr(info$anno, reg1+nadd, reg1+len1-2)

  # --------------------------------------------------------- 
  # look at geneMod
  # --------------------------------------------------------- 
    
  gnameNB   = NULL
  isoId     = NULL
  tnameNB   = NULL
  tabunNB   = NULL
  distNB    = NULL
  lenNB     = NULL
  status    = NULL

  clustNames = names(geneMod)
  length(clustNames)

  for(kk in 1:length(clustNames)){
    
    if(kk %% 100 == 0){
      cat(kk, date(), "\n")
    }
    
    cname1 = clustNames[kk]
    
    NB1    = geneMod[[cname1]]
    isoNB  = list()
    
    if(is.null(NB1$w2kp)){ next }
    
    # ------------------------------------------------------
    # extract NB model results of this gene
    # ------------------------------------------------------
    
    if(nrow(NB1$candiIsoform) != nrow(NB1$info)){
      stop("dimension of 'isoform' and 'info' do not match\n")
    }
    
    for(i in 1:length(NB1$w2kp)){
      is1        = NB1$candiIsoform[,NB1$w2kp[i]]
      isoNB[[i]] = NB1$info[which(is1==1),c("start", "end")]
    }
    
    gnameNB = c(gnameNB, rep(cname1, length(isoNB)))
    tabunNB = c(tabunNB, NB1$abundance)
    status  = c(status, rep(NB1$status, length(isoNB)))
    isoId   = c(isoId, 1:length(isoNB))
    
    # ------------------------------------------------------
    # comparison between NB model and annoation
    # ------------------------------------------------------
    
    lenV   = rep(NA, length(isoNB))
    distV  = rep(NA, length(isoNB))
    tnameV = rep(NA, length(isoNB))
    
    inChr  = info[info$chr==NB1$info$chr[1],]
    
    # ------------------------------------------------------
    # for each isoform in isoNB, find the transcript 
    # from annotation with the smallest distance
    # ------------------------------------------------------
    
    for(i in 1:length(isoNB)){
      iDis    = Inf
      iso1    = isoNB[[i]]
      pos1    = min(iso1$start)
      pos2    = max(iso1$end)
      lenV[i] = sum(iso1$end - iso1$start + 1)

      # ----------------------------------------------------
      # extract annotation overlapping this transcript
      # ----------------------------------------------------
      
      wExons  = which(inChr$start<=pos2 & inChr$end>=pos1)
      
      if(length(wExons) == 0){ next }
      
      info1   = inChr[wExons,]
      uTIds   = unique(info1$tranId)
      isoAnno = list()
      
      for(k in 1:length(uTIds)){
        u1            = uTIds[k]
        isoAnno[[u1]] = info1[info1$tranId==u1,c("start", "end")]
      }
      
      for(j in 1:length(isoAnno)){
        
        ijDis = distSet(iso1, isoAnno[[j]])
        
        if(ijDis < iDis){ 
          iDis       = ijDis 
          tnameNBi   = names(isoAnno)[j]
        }
        
      }
      
      distV[i]  = iDis
      tnameV[i] = tnameNBi
    }
        
    tnameNB = c(tnameNB, tnameV)
    distNB  = c(distNB,  distV)
    lenNB   = c(lenNB,   lenV)
  }
  
  r1 = data.frame(clustId=gnameNB, isoId=isoId, len=lenNB, abundance=tabunNB, 
                  status=status, nearestTranscript=tnameNB, dist=distNB, 
                  stringsAsFactors=FALSE)
  
  r1
}

