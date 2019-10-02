compareAnno.cuffcompare <-
function(cuffcompareGtf, gtf){
  
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
    
  cuff = read.table(cuffcompareGtf, sep="\t", as.is=TRUE)

  dim(cuff)
  cuff[1:2,]

  names(cuff) = c("chr", "source", "feature", "start", "end", 
  "score", "strand", "frame", "anno")

  table(cuff$chr)
  table(cuff$feature)
  
  cuff = cuff[which(cuff$feature == "exon"),]
  
  # --------------------------------------------------------- 
  # obtain gene_id
  # --------------------------------------------------------- 

  reg1   = regexpr('gene_id\\s(\\S+);', cuff$anno, perl=TRUE)
  len1   = attributes(reg1)[[1]]
  nadd   = length(unlist(strsplit("gene_id", split=""))) + 1
  cuff$geneId = substr(cuff$anno, reg1+nadd, reg1+len1-2)

  # --------------------------------------------------------- 
  # comparison between annotation and cufflinks estimates
  # ---------------------------------------------------------
  
  geneIds = unique(info$geneId)
  length(geneIds)
  
  distCuff  = gnameCuff = rep(NA, length(geneIds))
  
  for(kk in 1:length(geneIds)){
    
    if(kk %% 1000 == 0){
      cat(kk, date(), "\n")
    }
    
    gIdk = geneIds[kk]
    
    # ------------------------------------------------------
    # extract exon annotation of this gene
    # ------------------------------------------------------
    
    infow   = which(info$geneId == gIdk)
    if(length(infow) == 0){
      stop("man... cannot find this:", gIdk, "\n")
    }
    
    info1   = info[infow,]
    
    # ------------------------------------------------------
    # comparison between gene model and cuffcompare 
    # annotation, find the name of the closest gene nand the 
    # corresponding distance
    # ------------------------------------------------------

    iDis    = Inf      

    infoIso = info1[,c("start", "end")]
    pos1    = min(info1$start)
    pos2    = max(info1$end)
    
    cuffChr = cuff[cuff$chr==info1$chr[1],]
    wExons  = which(cuffChr$start<=pos2 & cuffChr$end>=pos1)
    
    if(length(wExons) == 0){ next }
    
    cuff1   = cuffChr[wExons,]
    ugIds   = unique(cuff1$geneId)
    
    isoCuff = list()
    
    for(k in 1:length(ugIds)){
      u1            = ugIds[k]
      isoCuff[[u1]] = cuff1[cuff1$geneId==u1,c("start", "end")]
    }
    
    for(j in 1:length(isoCuff)){
      
      ijDis = distSet(infoIso, isoCuff[[j]])
      
      if(ijDis < iDis){ 
        iDis       = ijDis 
        gnameCuffi = names(isoCuff)[j]
        if(length(gnameCuffi) > 1){
          stop("more than one gene ID for one transcript ID\n")
        }
      }
      
    }
    
    distCuff[kk]  = iDis
    gnameCuff[kk] = gnameCuffi
  }

  r1 = data.frame(geneId=geneIds, cuffId=gnameCuff, dist=distCuff)
  
  r1
}

