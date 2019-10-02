compareAnno.cufflinks <-
function(cufflinksGtf, isoTracking, gtf, knownIso){
  
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
  # read cufflinks results
  # --------------------------------------------------------- 
    
  cuff = read.table(cufflinksGtf, sep="\t", as.is=TRUE)

  dim(cuff)
  cuff[1:2,]

  names(cuff) = c("chr", "source", "feature", "start", "end", 
  "score", "strand", "frame", "anno")

  table(cuff$chr)
  table(cuff$feature)
  
  # --------------------------------------------------------- 
  # obtain transcript_id
  # --------------------------------------------------------- 

  reg1   = regexpr('transcript_id\\s(\\S+);', cuff$anno, perl=TRUE)
  len1   = attributes(reg1)[[1]]
  nadd   = length(unlist(strsplit("transcript_id", split=""))) + 1
  cuff$tranId = substr(cuff$anno, reg1+nadd, reg1+len1-2)

  # --------------------------------------------------------- 
  # split the resutls for transcript and exon
  # --------------------------------------------------------- 

  table(cuff$feature)

  cuffT = cuff[which(cuff$feature=="transcript"),]
  dim(cuffT)

  cuffE = cuff[which(cuff$feature=="exon"),]
  dim(cuffE)

  length(unique(cuffT$tranId))
  
  if(any(sort(unique(cuffE$tranId)) != sort(cuffT$tranId))){
    stop("mismatched transcript Id\n")
  }

  # --------------------------------------------------------- 
  # combine with tracking information
  # --------------------------------------------------------- 
    
  track = read.table(isoTracking, sep="\t", as.is=TRUE, header=TRUE)
  
  dim(track)
  track[1:2,]
  
  if(any(sort(track$tracking_id) != sort(cuffT$tranId))){
    stop("mismatched transcript Id\n")
  }
  
  # --------------------------------------------------------- 
  # comparison between annotation and cufflinks estimates
  # ---------------------------------------------------------
  
  if(!knownIso){
    distCuff  = tnameCuff = gnameCuff = rep(NA, nrow(cuffT))
    
    for(kk in 1:nrow(track)){
      
      if(kk %% 1000 == 0){
        cat(kk, date(), "\n")
      }
      
      tname1 = track$tracking_id[kk]
      
      # ------------------------------------------------------
      # extract cufflink results of this gene
      # ------------------------------------------------------
      
      cuffw   = which(cuffE$tranId == tname1)
      if(length(cuffw) == 0){
        stop("man... cannot find this:", tname1, "\n")
      }
      
      cuff1   = cuffE[cuffw,]
      
      # ------------------------------------------------------
      # comparison between cufflinks model and annotation, 
      # find the name of the closest transcript name and  
      # the corresponding distance
      # ------------------------------------------------------
      
      inChr   = info[info$chr==cuff1$chr[1],]
      cuffIso = cuff1[,c("start", "end")]
      
      iDis    = Inf      
      pos1    = min(cuff1$start)
      pos2    = max(cuff1$end)
      
      # ----------------------------------------------------
      # extract annotation of this gene
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
        
        ijDis = distSet(cuffIso, isoAnno[[j]])
        
        if(ijDis < iDis){ 
          iDis       = ijDis 
          tnameCuffi = names(isoAnno)[j]
          gnameCuffi = unique(info1$geneId[info1$tranId==tnameCuffi])
          if(length(gnameCuffi) > 1){
            stop("more than one gene ID for one transcript ID\n")
          }
        }
        
      }
      
      distCuff[kk]  = iDis
      tnameCuff[kk] = tnameCuffi
      gnameCuff[kk] = gnameCuffi
    }
    
  }else{
    tnameCuff = track$tracking_id
    gnameCuff = track$gene_id
    distCuff  = rep(0, nrow(track))
  }

  
  r1 = data.frame(cuffTranId=track$tracking_id, cuffGeneId=track$gene_id,
                  tranId=tnameCuff, geneId=gnameCuff, dist=distCuff, 
                  len=track$length, FPKM=track$FPKM, confLo=track$FPKM_conf_lo,
                  confHi=track$FPKM_conf_hi, status=track$FPKM_status,
                  stringsAsFactors=FALSE)
  
  r1
}

