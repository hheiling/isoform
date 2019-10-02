
isoDu <-
function(tags, gL, xData, readDepth, outputFileName, pdDistL, g2test,
  readLen, method=c("bootstrap", "permutation"), nlambda=10, ntau=3,
  lmax=500, duOnly=FALSE, nResample=10000, pUpper=0.05,
  pConf=0.99, misAligns=rep(0.001, length(tags)), muMin = 0.01, 
  maxTime=Inf, saveRData="", simulationPhi=NULL, simuNULL=FALSE)
{

  # --------------------------------------------------------- 
  # given isoform, calculate design matrix 
  # --------------------------------------------------------- 
  
  getX <- function(exons, isoforms, exonLens, readLen, pdDist, misAlign){
    
    X = matrix(0, nrow=length(exons), ncol=ncol(isoforms))
    physicalLen = rep(0, nrow(X))
    
    ## iterate for each feature
    for(i in 1:nrow(X)){
      eSeti = exons[[i]]
      physicalLen[i] = sum(exonLens[eSeti])
      
      if(all(diff(eSeti)<=1)){
        ## if no exon is skipped in this feature
        wmat      = which(colSums(isoforms[eSeti,,drop=FALSE]) == length(eSeti))
        ids       = 1:length(eSeti)
        X[i,wmat] = effLen(ids, exonLens[eSeti], readLen, pdDist, lmax)
      }else{
        ## if some exon(s) are skipped in this exon set
        ## there are two cases:
        ##
        ## (1) the isoform matches this exon set, i.e., the exons 
        ## skipped in the exon set is also skipped in the isoform
        ##
        ## (2) a consecutive set of the exons are skipped in the
        ## exon set but kept in the isoform.
        ##
        
        eSkipped = setdiff(min(eSeti):max(eSeti), eSeti)
        sum1     = colSums(isoforms[eSeti,,drop=FALSE]) 
        sum0     = colSums(isoforms[eSkipped,,drop=FALSE]) 
        
        ## case (1)
        wCase1   = which(sum1==length(eSeti) & sum0==0)
        
        if(length(wCase1) > 0){
          ids         = 1:length(eSeti)
          X[i,wCase1] = effLen(ids, exonLens[eSeti], readLen, pdDist, lmax)
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
          # unqIso is the unique patterns among the skipped exons. 
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
              X[i,w3] = effLen(ids, exonLens[w2Exn], readLen, pdDist, lmax)
            }
            
          }
          
        }
        
      }
    }
    
    X = X + misAlign*physicalLen
    X
  }
  
  # --------------------------------------------------------- 
  # check parameters
  # --------------------------------------------------------- 
  
  if( (! is.character(tags)) ){
     stop("tags should be character vectors\n")
  }
  
  if( (! is.numeric(readDepth)) ){
    stop("readDepth should be numeric vectors\n")
  }

  if(length(readDepth) != length(tags)){
    stop("tags and readDepth should have the same length\n")
  }
  
  readDepth = readDepth/median(readDepth)
  
  ## nn is the number of individuals
  nn = length(tags)
  
  method = method[1]
  
  if(nn <= 4 && method=="permutation"){
    warning("sample is small, will do boostrap rather than permutation.\n")
    method = "bootstrap"
  }
  
  # ---------------------------------------------------------
  # log sum exp
  # --------------------------------------------------------- 
  
  logsumexp <- function(v){
    wv  = which.max(v)
    mv  = v[wv]
    res = sum(exp(v[-wv] - mv))
    lse = mv + log(1+res)
    lse
  }
  
  # -----------------------------------------------------------------
  # log likelihood of (zero inflated) Negative Binomial distribution
  # -----------------------------------------------------------------

  loglikNB <- function(phi, mu, y, muMin){
    mu[mu < muMin] = muMin
    th = 1/phi
    lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(th) +
      y * log(mu + (y == 0)) - (th + y) * log(th + mu)
  }

  
  loglikZINB <- function(pZI, phi, mu, y, muMin){
    mu[mu < muMin] = muMin
    th   = 1/phi
    lk1  = lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(th) +
           y * log(mu + (y == 0)) - (th + y) * log(th + mu)
    wZI1 = which(pZI==1)
    wZIp = which(pZI>0 & pZI<1)
    
    if(length(wZI1) > 0){
      lk1[wZI1]  = 0
    }
    
    if(length(wZIp) > 0){
      for(w1 in wZIp){
        lk1[w1] = logsumexp(c(log(pZI[w1]), log(1 - pZI[w1]) + lk1[w1]))
      }
    }
    lk1
  }
  
  # -----------------------------------------------------------------
  # log likelihood of (zero inflated) Poisson distribution
  # -----------------------------------------------------------------

  loglikPoisson <- function(mu, y, muMin){
    mu[mu < muMin] = muMin
    y*log(mu) - mu - lgamma(y + 1)
  }
  
  loglikZIP <- function(pZI, mu, y, muMin){
    mu[mu < muMin] = muMin
    lk1  = y*log(mu) - mu - lgamma(y + 1)
    wZI1 = which(pZI==1)
    wZIp = which(pZI>0 & pZI<1)
    
    if(length(wZI1) > 0){
      lk1[wZI1]  = 0
    }
    
    if(length(wZIp) > 0){
      for(w1 in wZIp){
        lk1[w1] = logsumexp(c(log(pZI[w1]), log(1 - pZI[w1]) + lk1[w1]))
      }
    }
    lk1
  }
  
  # --------------------------------------------------------- 
  # checking data
  # --------------------------------------------------------- 
  
  ng  = length(g2test)

  sep = "\n------------------------------------------------------------------\n"
  ann = "testing differential isoform usage for"
  ann = paste(ann, length(g2test), "trancription clusters")
  message(sep, ann, sep)
  
  # ---------------------------------------------------------
  # checking xData
  # ---------------------------------------------------------
  
  if(is.vector(xData)){
    if(length(xData) != nn){
      stop("length of xData does not match sample size ", nn, "\n")
    }
    
    xData = xData - min(xData)
    xData = xData/max(xData)
    xData = matrix(xData, ncol=1)
  }else if(is.matrix(xData)){
    if(nrow(xData) != nn){
      stop("number of rows of xData does not match sample size ", nn, "\n")
    }
    
    for(xjj in 1:ncol(xData)){
      xData[,xjj] = xData[,xjj] - min(xData[,xjj])
      xData[,xjj] = xData[,xjj]/max(xData[,xjj])
    }
    
  }else{
    stop("xData must be a numeric vector or a data matrix\n")
  }

  # ---------------------------------------------------------
  # Start to analyze each transcript cluster one by one
  # ---------------------------------------------------------

  maxLRS = pvalBoots = nB = rep(NA, length(g2test))
  ns  = ps = k0s = k1s = k0p = k1p = rep(NA, length(g2test))

  for(i in 1:length(g2test)){
    
    sep = "\n----------------------------------------------------------\n"
    message(sep, "trancription cluster ", i, " @ ", date(), sep, appendLF=FALSE)
    
    gi       = g2test[i]
    exonLens = gL[[1]][[gi]]$info$end - gL[[1]][[gi]]$info$start
    ## kk is the number of exons
    kk       = length(exonLens)
    
    # -----------------------------------------------------------
    # obtain the union of exon sets acros samples, then generate
    # the fragment count matrix for each exon set and each sample
    # -----------------------------------------------------------
    
    exns = gL[[1]][[gi]]$count$exons
    for(ii in 2:nn){
      exns = union(exns, gL[[ii]][[gi]]$count$exons)
    }
    
    ctsI = matrix(0, nrow=length(exns), ncol=nn)
    
    for(ii in 1:nn){
      mat1 = match(gL[[ii]][[gi]]$count$exons, exns)
      ctsI[mat1,ii] = gL[[ii]][[gi]]$count$count
    }
    
    sumCtsI = colSums(ctsI)
    sumCtsI = sumCtsI/median(sumCtsI)

    # -------------------------------------------------------
    # obtain the union of candidate isoforms acros samples,
    # then recalculate effective length for certain exon sets
    # -------------------------------------------------------  
        
    candIso  = gL[[1]][[gi]]$candiIsoform
    for(ii in 2:nn){
      candIsoii = gL[[ii]][[gi]]$candiIsoform
      
      if(nrow(candIsoii) != kk){
        stop("Number of exon sets mismatched: ", nrow(candIsoii), "!=", kk, "\n")
      }
      
      candIso = unique(cbind(candIso, candIsoii), MARGIN=2)
    }
    dim(candIso)
    
    exons = lapply(strsplit(exns, ";"), function(v){sort(as.numeric(v))})
    
    ## mm is the number of exon sets, pp is the number of candiate isoforms
    mm    = length(exons)
    pp    = ncol(candIso)
    
    str1  = sprintf("\ncalculating effective length for %d isoforms", pp)
    str1  = sprintf("%s from %d exon sets", str1, mm)
    message(str1)
    XiiL  = list()
    
    for(ii in 1:nn){
      message(".", appendLF=FALSE)
      if(ii %% 70 == 0){ message("\n", appendLF=FALSE) }

      XiiL[[ii]] = getX(exons, candIso, exonLens, readLen, pdDistL[[ii]],
                        misAligns[ii])
    }
    message("\n", appendLF=FALSE)
    
    # -------------------------------------------------------
    # only keep those exon sets whose effect length > 0 for
    # at least one candidate isoform
    # ---------------------------------------------------------
    
    maxELen = sapply(XiiL, function(mm){apply(mm, 1, max)} )
    maxELen = apply(maxELen, 1, max)
    
    physicalLen = rep(0, mm)
    for(jj in 1:mm){
      eSetj = exons[[jj]]
      physicalLen[jj] = sum(exonLens[eSetj])
    }
    
    wEet2use = which(maxELen > physicalLen*min(misAligns))
    mm       = length(wEet2use)
    
    X0 = matrix(NA, nrow=mm*nn, ncol=pp)
    
    if(duOnly){
      for(ii in 1:nn){
        X0[((ii-1)*mm+1):(ii*mm),] = XiiL[[ii]][wEet2use,]*sumCtsI[ii]
      }
    }else{
      for(ii in 1:nn){
        X0[((ii-1)*mm+1):(ii*mm),] = XiiL[[ii]][wEet2use,]*readDepth[ii]
      }
    }

    ctsI = ctsI[wEet2use,]

    # -------------------------------------------------------
    # run isoSelect under H0
    # -------------------------------------------------------
    
    ns[i]   = nn*mm
    ps[i]   = ncol(X0)

    y   = as.numeric(ctsI)
    
    ## if this is a simulation
    if(!(is.null(simulationPhi))){
      beta0 = runif(ncol(X0), 0, 1)
      beta0[beta0 < 0.5] = 0
      
      gjj   = xData[,1]
      gjj   = rep(gjj, each=mm)
      
      X1    = cbind(X0*(1-gjj), X0*gjj)
      
      beta1 = runif(ncol(X1), 0, 1)
      beta1[beta1 < 0.5] = 0
      
      if(simuNULL){
        mu0 = X0 %*% beta0
      }else{
        mu0 = X1 %*% beta1
      }
      
      y = rnegbin(length(mu0), mu0, simulationPhi)
    }
    
    if(nrow(X0) > 2*ncol(X0)){
      useBIC = "BIC"  
    }else{
      useBIC = "extended BIC"  
    }
    
    g0 = isoSelectZI(y, X0, muMin, useBIC, nlambda=nlambda, ntau=ntau)
    
    if(is.null(g0)){ next }
        
    mu0  = g0$mu
    phi0 = g0$phi
    pZI0 = g0$pZI
    
    if(g0$family == "Negative Binomial"){
      if(max(g0$pZI) < 1e-16){
        lR0  = loglikNB(g0$phi, g0$mu, y, muMin)
      }else{
        lR0  = loglikZINB(g0$pZI, g0$phi, g0$mu, y, muMin)
      }
    }else{
      if(max(g0$pZI) < 1e-16){
        lR0  = loglikPoisson(g0$mu, y, muMin)
      }else{
        lR0  = loglikZIP(g0$pZI, g0$mu, y, muMin)
      }
    }
    
    k0s[i]  = g0$n2use
    
    # -------------------------------------------------------  
    # run isoSelect under alternative
    # -------------------------------------------------------  
    
    maxLRSi = -Inf
    
    for(xjj in 1:ncol(xData)){
      
      gjj = xData[,xjj]
      gjj = rep(gjj, each=mm)
      
      X1 = cbind(X0*(1-gjj), X0*gjj)
      g1 = isoSelectZI(y, X1, muMin, useBIC, nlambda=nlambda, ntau=ntau)
      
      if(is.null(g1)){ next }
            
      if(g1$family == "Negative Binomial"){
        if(max(g1$pZI) < 1e-16){
          lR1  = loglikNB(g1$phi, g1$mu, y, muMin)
        }else{
          lR1  = loglikZINB(g1$pZI, g1$phi, g1$mu, y, muMin)
        }
      }else{
        if(max(g1$pZI) < 1e-16){
          lR1  = loglikPoisson(g1$mu, y, muMin)
        }else{
          lR1  = loglikZIP(g1$pZI, g1$mu, y, muMin)
        }
      }
      
      lRSjj = sum(lR1 - lR0)
      
      if(lRSjj > maxLRSi){
        maxLRSi  = lRSjj
        k1s[i]   = g1$n2use
        g1Max    = g1
        g1Max$Xdata = xData[,xjj]
        g1Max$Xindx = xjj
      }
    }
    
    if(is.infinite(maxLRSi)){
      next
    }else{
      maxLRS[i]   = maxLRSi
      
      if(saveRData!=""){
        designX            = X1[,g1Max$w2use,drop=FALSE]
        g1Max$abundance    = g1Max$b2use*colSums(designX)
        g1Max$candiIsoform = candIso
        
        save(g1Max, file=sprintf("%s_%s.RData", saveRData, gi))
      }
    }
    
    # -------------------------------------------------------  
    # Bootstrapping...
    # -------------------------------------------------------  
    
    if(nResample > 0){
      # -----------------------------------------------------
      # estimate null distribution by parametric bootstrap
      # or permutation
      # -----------------------------------------------------
      
      LRSBoots = rep(NA, nResample)
      time0    = proc.time()[3]
      ks       = matrix(NA, nrow=nResample, ncol=2)
      
      for(j in 1:nResample){
        # ---------------------------------------------------
        # simulate read counts
        # ---------------------------------------------------
        
        if(method=="bootstrap"){
          
          if(g0$family=="Negative Binomial"){
            yb = rnegbin(rep(1,nn*mm), mu0, 1/phi0)
          }else{
            yb = rpois(nn*mm, mu0)
          }
          
          ru1  = runif(nn*mm)
          ww0  = which(ru1 <= pZI0)

          if(length(ww0)>0){
            yb[ww0] = 0
          }
          
          # -------------------------------------------------
          # liklihood under null
          # -------------------------------------------------
          
          g0b = isoSelectZI(yb, X0, muMin, useBIC, nlambda=nlambda, ntau=ntau)
          
          if(is.null(g0b)){ next }
          
          k0  = g0b$n2use
          
          if(g0b$family == "Negative Binomial"){
            if(max(g0b$pZI) < 1e-16){
              lR0  = loglikNB(g0b$phi, g0b$mu, yb, muMin)
            }else{
              lR0  = loglikZINB(g0b$pZI, g0b$phi, g0b$mu, yb, muMin)
            }
          }else{
            if(max(g0b$pZI) < 1e-16){
              lR0  = loglikPoisson(g0b$mu, yb, muMin)
            }else{
              lR0  = loglikZIP(g0b$pZI, g0b$mu, yb, muMin)
            }
          }
          
          xDataResample = xData
          
        }else if(method=="permutation"){
          
          p1 = sample(nn)
          xDataResample = xData[p1,,drop=FALSE]
          yb = y
          k0 = g0$n2use
          
        }
        
        # -----------------------------------------------------
        # run isoSelect under alternative
        # -----------------------------------------------------
        
        maxLRSi = -Inf

        for(xjj in 1:ncol(xDataResample)){
          
          gjj = xDataResample[,xjj]
          gjj = rep(gjj, each=mm)
          
          X1  = cbind(X0*(1-gjj), X0*gjj)
          g1b = isoSelectZI(yb, X1, muMin, useBIC, nlambda=nlambda, ntau=ntau)
          
          if(is.null(g1b)){ next }
          
          if(g1b$family == "Negative Binomial"){
            if(max(g1b$pZI) < 1e-16){
              lR1  = loglikNB(g1b$phi, g1b$mu, yb, muMin)
            }else{
              lR1  = loglikZINB(g1b$pZI, g1b$phi, g1b$mu, yb, muMin)
            }
          }else{
            if(max(g1b$pZI) < 1e-16){
              lR1  = loglikPoisson(g1b$mu, yb, muMin)
            }else{
              lR1  = loglikZIP(g1b$pZI, g1b$mu, yb, muMin)
            }
          }
          
          lRSjj  = sum(lR1 - lR0)
          
          if(lRSjj > maxLRSi){
            maxLRSi = lRSjj
            k1      = g1b$n2use
          }
        }

        if(is.finite(maxLRSi)){
          LRSBoots[j] = maxLRSi
          ks[j,] = c(k0, k1)
        }
        
        time1 = proc.time()[3]
        
        if(time1 - time0 > maxTime){ break }
        
        if(j %% 100 == 0){
          message(" ", j, " ", date())
          ntest = length(na.omit(LRSBoots))
          kpass = length(which(LRSBoots >= maxLRS[i]))
          ptest = pbinom(kpass, ntest, pUpper, lower.tail=FALSE)
          
          if(ptest < 1-pConf){ break }
        }
      }
      
      nB[i]   = length(na.omit(LRSBoots))
      pvalBoots[i] = (length(which(LRSBoots >= maxLRS[i]))+1)/(nB[i] + 1)
      
      k0p[i]  = median(ks[,1], na.rm=TRUE)
      k1p[i]  = median(ks[,2], na.rm=TRUE)
    }
    
  }
  
  
  maxLRS    = signif(maxLRS, 5)
  pvalBoots = signif(pvalBoots, 5)
  
  if(nResample > 0){
    db = data.frame(cluster=g2test, maxLRS=maxLRS, n=ns, 
                    p=ps, k0=k0s, k1=k1s, k0p=k0p, k1p=k1p, 
                    pResample=pvalBoots, nResample=nB)
  }else{
    db = data.frame(cluster=g2test, maxLRS=maxLRS, n=ns,  
                    p=ps, k0=k0s, k1=k1s, k0p=k0p, k1p=k1p)
  }
  
  write.table(db, file = outputFileName, append = FALSE, quote = FALSE, 
              sep = "\t", row.names = FALSE, col.names = TRUE)
  
}

