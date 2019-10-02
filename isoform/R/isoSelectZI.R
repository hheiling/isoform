isoSelectZI <-
function(y, X, muMin, useBIC, nlambda, ntau, maxIt=5,
        min_pi_ZI=0.001, trace=0){
    
    g0  = isoSelect(y, X, modelSelection=useBIC, nlambda=nlambda, ntau=ntau)
    
    if(is.null(g0$w2use)){ return(NULL) }
    
    design = X[, g0$w2use, drop = FALSE]
    mu0    = design %*% g0$b2use
    mu0[which(mu0 < muMin)] = muMin
    phi0   = g0$phi
    
    if(g0$family=="Negative Binomial"){
      db1  = dnbinom(0, size=1/phi0, mu=mu0)
    }else{
      db1  = dpois(0, lambda=mu0)
    }
    
    wy0    = which(y==0)
    ## pi1 : the probability of zero infation
    ## pZI0: posterior probability of ZI component
    pi1    = max(min_pi_ZI, (length(wy0) - sum(db1))/length(y))
    pZI0   = as.numeric(pi1*(y==0)/(pi1 + (1-pi1) *db1))
    maxDi  = 1.0
    kk     = 0
    
    while(pi1 > 0 && maxDi > 1e-5 && kk < maxIt){
      kk = kk + 1
      # -------------------------------------------------------  
      # update model fit by incoporating zero inflated component
      # -------------------------------------------------------  
      
      g0  = isoSelect(y, X, priorWeights=1-pZI0, modelSelection=useBIC, 
                      nlambda=nlambda, ntau=ntau)
      
      if(is.null(g0$w2use)){ return(NULL) }
      
      design = X[, g0$w2use, drop = FALSE]
      mu0    = design %*% g0$b2use
      mu0[which(mu0 < muMin)] = muMin
      phi0   = g0$phi
      
      if(g0$family=="Negative Binomial"){
        db1  = dnbinom(0, size=1/phi0, mu=mu0)
      }else{
        db1  = dpois(0, lambda=mu0)
      }
      wy0    = which(y==0)
      ## pi1 = the probability of zero infation
      ## prw = posterior probability of none ZI component
      pi1    = max(min_pi_ZI, (length(wy0) - sum(db1))/length(y))
      pZI    = as.numeric(pi1*(y==0)/(pi1 + (1-pi1) *db1))
      maxDi  = max(abs(pZI0 - pZI))
      pZI0   = pZI
      
      if(trace){
        message(sprintf("iteration %d: pi1=%e, maxDi=%e", kk, pi1, maxDi))
      }
      
    }
    
    g0[["pZI"]] = pZI
    g0[["mu"]]  = mu0
    
    g0
  }
