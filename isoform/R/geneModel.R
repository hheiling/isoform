geneModel <-
function(gene, d, pdDist, isoforms=NULL, lmax=length(pdDist), 
         maxBreaks=5, pvalBreaks=0.05, pvalExpress=0.01, foldExpress=1/5,  
         isoformOnly=FALSE, runIsoSelect=TRUE, eLenMin=1.0, muMin=0.01, 
         pMaxRel=10, pMaxAbs=2000, verbose=1)
{
  status = "OK"

  # --------------------------------------------------------- 
  # calculate the effective length of each exon
  # --------------------------------------------------------- 
    
  effLen3 <- function(rj, d, pdDist, lmax){
    lRange  = d:lmax
    sum(pdDist[lRange]*(rj - 1 + lRange)) 
  }

  # --------------------------------------------------------- 
  # construct isoforms for one gene
  # --------------------------------------------------------- 
  
  iso1 <- function(exonsJ, infoJ, countJ, effLnJ, d, pdDist, lmax, 
                   maxBreaks, pvalBreaks, pvalExpress, foldExpress,
                   pMaxRel, pMaxAbs){
    
    statusJ = "OK" 

    # ------------------------------------------------------- 
    # first, identify all the isofroms with consecutive exons
    # -------------------------------------------------------
    
    nn    = length(exonsJ)
    nExnJ = nrow(infoJ)
    
    # if there is only one exon set, return the only isoform
    if(nn==1){
      isoJ = matrix(1,nrow=1,ncol=1)
      return(list(isoform=isoJ, status=statusJ))
    }
    
    isoJ  = NULL

    # --------------------------------------------------------- 
    # for each exon, choose the connected exons, and count
    # the number of reads that overlap with each exon
    # --------------------------------------------------------- 
    
    ancester = offspring = list()
    eCounts  = rep(0, nExnJ)
    exn2RM   = numeric(0)
    
    for(ex1 in 1:nExnJ){
      wmat = which(sapply(exonsJ, function(x){ ex1 %in% x }))
      
      if(length(wmat)>0){
        eCounts[ex1] = sum(countJ[wmat])
        tmp = unique(unlist(exonsJ[wmat]))
        w1  = which(tmp<ex1)
        w2  = which(tmp>ex1)
        ancester[[ex1]]  = tmp[tmp<ex1]
        offspring[[ex1]] = tmp[tmp>ex1]
      }else{
        ancester[[ex1]]  = -1
        offspring[[ex1]] = -1
        exn2RM = c(exn2RM, ex1)
      }
    }
    
    starts = which(sapply(ancester, length)==0)
    ends   = which(sapply(offspring, length)==0)
        
    if(length(starts)==0){
      stop("no starting exon is found\n")  
    }
    
    if(length(ends)==0){
      stop("no ending exon is found\n")  
    }
    
    pLens = infoJ$end - infoJ$start + 1
    eLens = rep(NA, length(pLens))
    
    for(i in 1:nExnJ){
      eLens[i] = effLen3(pLens[i], d, pdDist, lmax)
    }
    
    # --------------------------------------------------------- 
    # calculate the p-value of differential isoform expression
    # --------------------------------------------------------- 
    
    pBreaks = rep(1, nExnJ)
    
    for(i in 2:nExnJ){
      x = eCounts[(i-1):i]
      if(any(x==0)){ next }
      p = eLens[(i-1):i]
      p = p/sum(p)
      pBreaks[i] = chisq.test(x, p=p)$p.value
    }
    
    where2Break = which(pBreaks < pvalBreaks)
    
    if(length(where2Break)>0){
      pvalsB = pBreaks[where2Break]
      odrB   = order(pvalsB)
      where2Break = where2Break[odrB]
      
      if(length(where2Break) > maxBreaks){
        where2Break = where2Break[1:maxBreaks]
      }
      
      starts = c(starts, where2Break)
      ends   = c(ends,   where2Break-1)
    }
    
    # ----------------------------------------------------- 
    # consecutive isoforms with alternative starts and ends
    # ----------------------------------------------------- 

    for(i in starts){
      for(j in ends){
        
        if(i >= j){ next }
        
        isoij = rep(0, nExnJ)
        isoij[i:j] = 1
        isoJ = cbind(isoJ, isoij, deparse.level = 0)
        
      }
    }
    
    # ----------------------------------------------------- 
    # remove those exons with no reads at all
    # ----------------------------------------------------- 
    
    if(length(exn2RM) > 0){
      isoJ[exn2RM,] = 0
    }
    
    # -------------------------------------------------------
    # find those exon-skipping events that have enough 
    # read-depth 
    # ------------------------------------------------------- 
    
    wSkip  = which(sapply(exonsJ, function(v){ any(diff(v) > 1) }))
    
    # adjust exon length
    effLnJ[effLnJ < 0.1] = 0.1
    effLnJ = effLnJ + eLenMin
    
    # average density
    totalCnt = sum(countJ)
    totalLen = sum(effLnJ)
    
    pvalExps = rep(NA, nn)
    
    for(kk in 1:length(exonsJ)){
      pvalExps[kk] = pbinom(countJ[kk], totalCnt, effLnJ[kk]/totalLen)
    }
    
    ratioExps  = foldExpress*totalCnt/totalLen
    wExpressed = which(pvalExps > pvalExpress | countJ/effLnJ > ratioExps)
    
    length(wExpressed)
    length(wSkip)

    length(intersect(wExpressed, wSkip))

    # -------------------------------------------------------
    # second, construct all the isoforms with exon-skipping
    # ------------------------------------------------------- 
    
    wSkip  = intersect(wSkip, wExpressed)
    pSkip  = pvalExps[wSkip]
    odSkip = order(pSkip, decreasing=TRUE)
    
    wSkip  = wSkip[odSkip]
    pSkip  = pSkip[odSkip]

    nSkip  = length(wSkip)
    
    if(nSkip > 0){
      
      for(w1 in wSkip){
        
        ## wExons are the set of exons not to drop
        wExons = exonsJ[[w1]]
        nEx2kp = length(wExons)
        checkW = colSums(isoJ[wExons,,drop=FALSE])
        
        ## wIsoms are the index of isoforms where we may skip some exons
        wIsoms = which(checkW == nEx2kp)
        if(length(wIsoms) == 0){ next }
        
        ## iso2 are the new set of isoforms after exon skipping
        iso2   = isoJ[,wIsoms,drop=FALSE]
        w2skip = setdiff((min(wExons):max(wExons)), wExons)
        iso2[w2skip,] = 0
        
        isoJ1 = isoJ
        isoJ  = cbind(isoJ, iso2)
        isoJ  = unique(isoJ, MARGIN=2)
        
        if(ncol(isoJ)/nn > pMaxRel || ncol(isoJ) > pMaxAbs){
          isoJ = isoJ1
          
          statusJ = "LOWDATA" 

          if(verbose > 1){
            message(sprintf("n=%d, p=%d, ", nn, ncol(isoJ)))
            message("skip other exon junctions due to large number of candidate isoJ\n")
          }
          
          break
        }
        
        
      }
    }
    
    list(isoform=isoJ, status=statusJ)
  }
  
  # --------------------------------------------------------- 
  # log likelihood of negative Binomial distribution
  # --------------------------------------------------------- 
    
  loglikNB <- function(phi, mu, y, muMin){
    mu[mu < muMin] = muMin
    th = 1/phi
    lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(th) +
         y * log(mu + (y == 0)) - (th + y) * log(th + mu)
  }
  
  # --------------------------------------------------------- 
  # log likelihood of Poisson distribution
  # --------------------------------------------------------- 
  
  loglikPoisson <- function(mu, y, muMin){
    mu[mu < muMin] = muMin
    y*log(mu) - mu - lgamma(y + 1)
  }
  
  # --------------------------------------------------------- 
  # check gene
  # --------------------------------------------------------- 
  
  if(!is.list(gene)){
    stop("gene must be a list\n")
  }
  
  if(any(! c("info", "count") %in% names(gene))){
    stop("gene must be a list with component 'exon' and 'junc'\n")
  }
  
  info      = gene$info
  count     = gene$count
  nn        = nrow(count)
  info$exon = as.numeric(info$exon)
  
  colnms = c("chr", "start", "end", "exon")
  if(any(! colnms %in% names(info))){
    stop("gene$info must be a list with components: ", colnms, "\n")
  }
  
  exonLens  = info$end - info$start

  colnms = c("count", "exons")
  if(any(! colnms %in% names(count))){
    stop("count must be a list with components: ", colnms, "\n")
  }
  
  # --------------------------------------------------------- 
  # other variables
  # --------------------------------------------------------- 
  
  if(length(d) != 1 || (!is.numeric(d)) || d < 0){
    stop("d must a positive number\n")  
  }
  
  if(length(pdDist) != lmax || (!is.numeric(pdDist))){
    stop("pdDist must be numeric vector of length lmax\n")  
  }
  
  if(abs(sum(pdDist) - 1) > 1e-10 ){
    stop("summation of pdDist must be 1\n")  
  }
  
  if(any(pdDist < 0)){
    stop("pdDist must be non-negative\n")  
  }
  
  # --------------------------------------------------------- 
  # choose possible starting exons and ending exons, if
  # isoforms are given, add some additional informative exons
  # --------------------------------------------------------- 
  
  nExons = length(info$exon)
  if(! all(info$exon == 1:nExons)){
    stop("exon IDs are not consecutive ordered numbers\n")
  }
  
  exons = count$exons
  
  if(!is.null(isoforms)){
    
    for(ki in 1:ncol(isoforms)){
      isoKi = isoforms[,ki]
      wKi1  = as.character(which(isoKi == 1))
      
      if(length(wKi1)==1){
        wKi2 = character(0)
      }else{
        wKi2  = paste(wKi1[-length(wKi1)], wKi1[-1], sep=";")
      }
      
      exons = union(exons, c(wKi1, wKi2))
    }
    
  }
  
  exons  = lapply(strsplit(exons, ";"), function(v){sort(as.numeric(v))})
  
  exn1st = sapply(exons, function(v){ v[1] })
  exn2nd = sapply(exons, function(v){ v[2] })
  exn3rd = sapply(exons, function(v){ v[3] })
  exnsum = sapply(exons, sum)

  exnodr = order(exn1st, exn2nd, exn3rd, exnsum, na.last=FALSE)
  exons  = exons[exnodr]
  
  # --------------------------------------------------------- 
  # obtain isoforms for gene cluster
  # --------------------------------------------------------- 
    
  if(is.null(isoforms)){
    count  = count[exnodr,,drop=FALSE]

    # unique gene names
    uGnms = unique(unlist(strsplit(info$gene, split=":")))
    
    ## calculate effective length
    eL = rep(0, length(exons))
    for(iE in 1:length(exons)){
      eSeti = exons[[iE]]
      nExi  = length(eSeti)
      eL[iE] = effLen(1:nExi, exonLens[eSeti], d, pdDist, lmax)
    }
    
    count$eL = eL
    
    if(length(uGnms) == 1){
      isoList = iso1(exons, info, count$count, count$eL, d, pdDist, lmax,  
                      maxBreaks, pvalBreaks, pvalExpress, foldExpress, 
                      pMaxRel, pMaxAbs)

      isoforms = isoList$isoform
      status   = isoList$status
      
      if(is.null(isoforms)){
        y  = count$count
        gm = list(info=info, count=count, y=y, candiIsoform=NULL, X=NULL)
        gm$design     = NULL
        gm$abundance  = NULL
        gm$isoform    = NULL
        gm$likelihood = NULL
        gm$status     = "FAIL"
        return(gm)
      }
      
    }else{
      
      isoforms = NULL
      
      for(j in 1:length(uGnms)){
        
        # construct isoforms for each gene

        uGnmj  = uGnms[j]
        wEx    = grep(uGnmj, info$gene)
        infoJ  = info[wEx,]
        infoJ$exon = 1:nrow(infoJ)

        exonsJ = list()
        countJ = rep(NA, nrow(count))
        effLnJ = rep(NA, nrow(count))

        kk     = 0
        
        for(k in 1:length(exons)){
          matk = match(exons[[k]], wEx)
          
          if(! any(is.na(matk))){
            kk = kk + 1
            exonsJ[[kk]] = matk
          }
          
          countJ[kk] = count$count[k]
          effLnJ[kk] = count$eL[k]
        }
        
        countJ = countJ[1:kk]
        effLnJ = effLnJ[1:kk]
        
        if(kk > 0){
          isoList = iso1(exonsJ, infoJ, countJ, effLnJ, d, pdDist, lmax, 
                         maxBreaks, pvalBreaks, pvalExpress, foldExpress, 
                         pMaxRel, pMaxAbs)
          
          isoformJ = isoList$isoform
          
          if(isoList$status == "LOWDATA"){
            status = "LOWDATA"
          }
          
          if(is.null(isoformJ)){
            isoforms = NULL
            break
          }
          
          isoformsjj = matrix(0, nrow=nExons, ncol=ncol(isoformJ))
          isoformsjj[wEx,] = isoformJ
          
          isoforms = cbind(isoforms, isoformsjj)
        }        
        
      }
      
      if(is.null(isoforms)){
        y  = count$count
        gm = list(info=info, count=count, y=y, candiIsoform=NULL, X=NULL)
        gm$design     = NULL
        gm$abundance  = NULL
        gm$isoform    = NULL
        gm$likelihood = NULL
        gm$status     = "FAIL"

        return(gm)
      }
      
    }
  }
  
  # -------------------------------------------------------
  # based on the isoforms, construct the covariates
  # each row of X is an observation and each column
  # is an isoform
  # ------------------------------------------------------- 
  
  if(verbose > 1){
    message(sprintf("# of isoforms=%d\n", ncol(isoforms)))
  }
  
  if(isoformOnly){
    return(isoforms)
  }
  
  X = matrix(0, nrow=length(exons), ncol=ncol(isoforms))
  
  ## iterate for each feature
  for(i in 1:length(exons)){
    eSeti = exons[[i]]
    nExi  = length(eSeti)
    
    if(all(diff(eSeti)<=1)){
      ## if no exon is skipped in this feature
      wmat      = which(colSums(isoforms[eSeti,,drop=FALSE]) == nExi)
      X[i,wmat] = effLen(1:nExi, exonLens[eSeti], d, pdDist, lmax)
      # X[i,wmat] = count$eL[i]
    }else{
      ## if some exon(s) are skipped in this exon set
      ## there are two cases:
      ##
      ## (1) the isoform matches this exon set, i.e., the exons 
      ## skipped in the exon set is also skipped in the isoform
      ##
      ## (2) one consecutive set of the exons that are skipped in the  
      ## exon set is kept in the isoform. 
      ##
      
      eSkipped = setdiff(min(eSeti):max(eSeti), eSeti)
      sum1     = colSums(isoforms[eSeti,,drop=FALSE]) 
      sum0     = colSums(isoforms[eSkipped,,drop=FALSE]) 
      
      ## case (1)
      wCase1   = which(sum1==length(eSeti) & sum0==0)
      
      if(length(wCase1) > 0){
        X[i,wCase1] = effLen(1:nExi, exonLens[eSeti], d, pdDist, lmax)
      }
      
      ## case (2)
      ## for single end read, situation (2) is impossible.
      ## for paired-end read, situation (2) is possible, 
      ## but it is impossible that two or more non-consecutive  
      ## sets of exons are skipped in the exon set, but kept in  
      ## the isoform. The function effLen will return 0 if there 
      ## are two or more break points in the exon set

      wCase2   = which(sum1==length(eSeti) & sum0>0)

      if(length(wCase2) > 0){
        # unqIso count the unique patterns among the skipped exons. 
        # for example, the exon set is 1-2-4, and
        # many of the isoforms have exons 1-2-3-4
        unqIso = unique(isoforms[eSkipped,wCase2,drop=FALSE], MARGIN=2)
        
        for(jj in 1:ncol(unqIso)){
          
          isojj = isoforms[eSkipped,wCase2,drop=FALSE]-unqIso[,jj]
          w3    = wCase2[which(colSums(abs(isojj))==0)]
          
          ## w3 index all those isoforms having the same
          ## pattern as unqIso[,jj] among wCase2.
          ## eSeti should have the same effective length  
          ## for any isoforms indexed by w3
          w2    = w3[1]
          w2Exn = which(isoforms[,w2]>0)
          ids   = match(eSeti, w2Exn)
          wbrk  = which(diff(ids) > 1)
          if(length(wbrk) <= 1){
            rjs     = exonLens[w2Exn]
            X[i,w3] = effLen(ids, rjs, d, pdDist, lmax)
          }
          
        }
        
      }
      
    }
  }
  
  w2rm     = which(apply(X, 2, function(v){all(v <= eLenMin)}))
  if(length(w2rm) > 0){
    isoforms = isoforms[,-w2rm, drop=FALSE]
    X        = X[,-w2rm, drop=FALSE]
  }
  
  X = X + eLenMin
  # X = unique(X, MARGIN=2)

  if(verbose > 1){
    message(sprintf("size of design matrix X = (%d,%d)\n", nrow(X), ncol(X)))
  }
  
  # --------------------------------------------------------- 
  # now, penalized regression
  # --------------------------------------------------------- 
  
  y  = rep(0, nrow(X))
  exnsStr = sapply(exons, paste, collapse=";")
  mm      = match(count$exon, exnsStr)
  
  y[mm]   = count$count
  countN  = data.frame(exons=exnsStr, count=y)
  
  gm = list(info=info, count=countN, y=y, candiIsoform=isoforms, X=X)

  if(runIsoSelect){
    if(ncol(X)==0){
      gm$w2kp       = NULL
      gm$status     = "FAIL"
      gm$abundance  = NULL
      gm$fitted     = NULL
      gm$likelihood = NULL
    }else if(nrow(X)==1){
      gm$w2kp       = 1
      gm$status     = "OK"
      gm$abundance  = y*(X[1,]/sum(X[1,]))
      gm$fitted     = NULL
      gm$likelihood = NULL
      gm$family       = "Poisson"
      gm$phi          = 0
    }else{
      
      g1 = isoSelect(y, X, muMin=muMin)
      
      if(ncol(X)/nrow(X) > pMaxRel || ncol(X) > pMaxAbs){
        status = "LOWDATA"
      }
      
      if(is.null(g1$w2use)){
        gm$w2kp       = NULL
        gm$status     = "FAIL"
        gm$abundance  = NULL
        gm$fitted     = NULL
        gm$likelihood = NULL
      }else{
        designX       = X[,g1$w2use,drop=FALSE]
        gm$w2kp       = g1$w2use
        gm$status     = status
        gm$abundance  = g1$b2use*colSums(designX)
        names(gm$abundance) = colnames(gm$isoform)

        mu            = designX %*% g1$b2use
        gm$fitted     = mu
        
        if(g1$family=="Poisson"){
          gm$likelihood = loglikPoisson(mu, y, muMin)
        }else if(g1$family=="Negative Binomial"){
          gm$likelihood = loglikNB(g1$phi, mu, y, muMin)
        }else{
          stop("invalid family: ", g1$family, "\n")
        }
      }
      
      gm$family = g1$family
      gm$phi    = g1$phi
    }
    
  }
  
  gm
}

