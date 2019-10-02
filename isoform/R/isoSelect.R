isoSelect <-
function(y, X, lambda=NULL, tau=NULL, nlambda=10, ntau=3,
         pMax=min(round(0.5*length(y)),100), offset=NULL, 
         naPercent=0.4, maxit=20, maxitIAL=20, conv=1e-5, 
         scoreTestP=0.05, trace=1, muMin=0.01, 
         priorWeights=NULL, modelSelection=NULL)
{

  if(!is.numeric(y)){
    stop("y must be a numeric vector\n")
  }
  
  if(!is.matrix(X)){
    stop("X must be a matrix\n")
  }
  
  M = ncol(X)
  N = length(y)
  
  if(N <= 2){
    warning("sample size is too small\n")
    return(NULL)
  }
  
  if(nrow(X) != N){
    stop("the dimension of X and y do not match\n")
  }
  
  if(!is.null(offset)){
    useOffset = 1
    if((! is.numeric(offset)) || length(offset) != N){
      stop("offset must be a numeric vector of the same length as y\n")
    }
  }else{
    useOffset = 0
    offset    = 0.0 
  }

  isNA = apply(X, 1, function(v){ any(is.na(v)) })
  isNA = is.na(y) | isNA
    
  if(length(which(isNA))/N > naPercent){
    cat(y)
    print(X)
    stop("percent of missing data is too high\n")
  }
  
  w2kp  = which(!isNA)
  
  if(length(w2kp) == 0){
    warning("isoSelect: there is no non-missing values\n")
    return(NULL)
  }
  
  yR    = y[w2kp]
  XR    = X[w2kp,,drop=FALSE]
  XNorm = sqrt(apply(XR, 2, function(v){sum(v*v)})/N)
  XR    = t(t(XR)/XNorm)
  
  # calculate sqrt weights assuming overdispersion is 1
  ws    = 1/(0.1 + yR)
  if(!is.null(priorWeights)){
    ws = ws*priorWeights
  }
  
  XN    = XR*ws
  yN    = yR*ws
  
  if(!is.null(lambda)){
    nlambda = length(lambda)
  }
  
  if(!is.null(tau)){
    ntau = length(lambda)
  }
  
  if(is.null(tau)){
    tau = exp(seq(log(0.1/sqrt(N)), log(0.1), length.out=ntau))
  }
  tau = rep(tau, each=nlambda)

  if(is.null(lambda)){
    ratioMax = max(100*abs(t(XN) %*% yN)/apply(XN, 2, function(v){sum(v*v)}))
    ratioMin = min(0.1/N, ratioMax/1e5)
    ratios   = exp(seq(log(ratioMax), log(ratioMin),  length.out=nlambda))
  
    lambda = tau*ratios
  }
  
  if(!is.null(lambda)){
    lambda = c(lambda, 0)
  }
  
  if(!is.null(tau)){
    tau    = c(tau, 0.1)
  }
  
  if(useOffset){
    offset = offset[w2kp]
  }

  N    = length(yR)
  init = 0

  dims     = numeric(11)
  dims[1]  = N
  dims[2]  = M
  dims[3]  = maxit
  dims[4]  = maxitIAL
  dims[5]  = init
  dims[6]  = useOffset
  dims[7]  = length(lambda)
  dims[8]  = pMax
  
  if(is.null(modelSelection)){
    dims[9]  = 0
  }else if(modelSelection == "BIC"){
    dims[9]  = 1
  }else if(modelSelection == "extended BIC"){
    dims[9]  = 2
  }else{
    stop("modelSelection must be NULL, BIC, or extended BIC\n")
  }
  
  if(is.null(priorWeights)){
    dims[10]  = 0
    priorWeights = rep(1, N)
  }else if(is.numeric(priorWeights) && length(priorWeights)==N){
    dims[10]  = 1
  }else{
    stop(sprintf("priorWeights must be a vector of length %d.\n", N))
  }
  
  nIter    = 0
  phi      = 0.0
  
  w2use    = rep(-9, pMax)
  b2use    = numeric(pMax)
  score    = rep(0, nlambda*ntau+1)
  n2kp     = rep(0, nlambda*ntau+1)
  convg    = 0
  
  n2use      = 0
  score2use  = 0.0
  lambda2use = 0.0
  tau2use    = 0.0
  family     = 0
  likelihood = rep(0, maxit)

  Z = .C("isoSelect", as.integer(dims), as.double(yR), as.double(muMin), 
         as.double(offset), as.double(XR), as.double(conv), 
         as.double(scoreTestP), as.integer(trace), as.double(lambda), 
         as.double(tau), nIter=as.integer(nIter), phi=as.double(phi),  
         as.double(priorWeights), n2use=as.integer(n2use), 
         w2use=as.integer(w2use), n2kp=as.integer(n2kp), 
         b2use=as.double(b2use), score=as.double(score), 
         score2use=as.double(score2use), lambda2use=as.double(lambda2use), 
         tau2use=as.double(tau2use), likelihood=as.double(likelihood),
         family=as.integer(family), PACKAGE="isoform")
  
  n2use=Z$n2use
  
  if(n2use > 0){
    w2use = Z$w2use[1:n2use] + 1
    b2use = Z$b2use[1:n2use]
    b2use = b2use/XNorm[w2use]
  }else{
    w2use = b2use = NULL
  }
  
  if (Z$nIter >= maxit){
    warning(sprintf("reach mamximum %d iterations\n", maxit))  
  }
  
  score  = data.frame(lambda=lambda, tau=tau, score=Z$score, n2kp=Z$n2kp)
  likelihood = Z$likelihood[1:Z$nIter]
  
  if(Z$family == 2){
    family = "Poisson"
  }else if(Z$family == 5){
    family = "Negative Binomial"
  }else if(Z$family == 0){
    str = "considering run isoSelect with extended tunning parameters\n"
    warning("no isoform is selected, ", str)
  }else{
    stop("invalid family: ", Z$family, "\n")
  }
  
  result = list(n2use=n2use, w2use=w2use, b2use=b2use, 
                score2use=Z$score2use, lambda2use=Z$lambda2use, 
                tau2use=Z$tau2use, score=score, nIter=Z$nIter, 
                family=family, phi=Z$phi, plogLikelihood=likelihood)
  
  class(result)="iSelect"
  
  result
}

