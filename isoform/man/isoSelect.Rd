\name{isoSelect}
\alias{isoSelect}
\title{
  select the most likely isoform combination for a gene
}
\description{
  select isoform subset by penalized generalized regression
}
\usage{
isoSelect(y, X, lambda=NULL, tau=NULL, nlambda=10, ntau=3,  
         pMax=min(round(0.5*length(y)),100), offset=NULL, 
         naPercent=0.4, maxit=20, maxitIAL=20, conv=1e-5, 
         scoreTestP=0.05, trace=1, muMin=0.01, 
         priorWeights=NULL, modelSelection=NULL)
}
\arguments{
	
  \item{y}{
	a vector of the response variable. }
	
  \item{X}{
	a matrix of covariate data, each column is one covariate and 
	each row is one observation. Intercept should not be included in 
	this covariate matrix. }

  \item{lambda}{
	one of the two penalization parameters for log penalty }
	
  \item{tau}{
	one of the two penalization parameters for log penalty }

  \item{nlambda}{
	number of values of tuning parameter lambda. Ignored if lambda is specified. }
	
  \item{ntau}{
	number of values of tuning parameter tau. Ignored if tau is specified. }

  \item{pMax}{
	the maximum number of covariates to be kept in the model }
	
  \item{offset}{
	a vector specifying offset in the linear predictor scale, eta = Xbeta + offset.}
	
  \item{naPercent}{
	the maximum percent of missing that are allowed in the data. }
	
  \item{maxit}{
	the maximum number of iteration of Iterated Re-weighted Least Square }
	 
  \item{maxitIAL}{
	the maximum number of iteration of for IAL }

  \item{conv}{
	convergence tolerance. Let wss be the weighted sum of Pearson residual square, 
	the algorithm converges if |wss-wss_last|/wss_last < conv }
		
  \item{scoreTestP}{ p-value cutoff for score test of the overdispersion in 
  Poisson model. For each association regression, we first fit a Poisson model, 
  and move to negative binomial if the score test p-value smaller than scoreTestP}

  \item{trace}{ print out tracking information if trace > 1.}
  
  \item{muMin}{
  if mu=0, set it to muMin. This is used during GLM fitting for stability of the 
  algorithm.}
  
  \item{priorWeights}{
  the prior weights of the observations. If it is NULL, all observations have 
  the same prior weights.}
  
  \item{modelSelection}{
  the information criterion to be used for choosing tuning parameter. Possible values
  include NULL, "BIC", and "extended BIC". If it is NULL, BIC will use used if sample 
  size is larger than the number of parameters. Otherwise extended BIC will be used. }

}

\value{
a list of the following components:
  \item{n2use }{number of covariates selected}

  \item{w2use }{the selected covariates}

  \item{b2use }{the coefficients of the selected covariates}

  \item{score2use }{the smallest (BIC) score}
  
  \item{lambda2use }{lambda value corresponding to the smallest score}
  
  \item{score }{(BIC) scores across all the tuning parameters}
  
  \item{nIter }{number of iterations}
  
  \item{family }{Poisson or Negative Binomial}

  \item{phi }{overdispersion parameter}

  \item{likelihood }{likelihoods across iterations}

}
\author{
Wei Sun weisun@email.unc.edu
}
\references{
  Sun et al. IsoDOT Detects Differential RNA-isoform expression with respect to 
  a Categorical or Continuous Covariate with High Sensitivity and Specificity.
}

\examples{
data(simuE5I6)

y  = simuE5I6$y
X  = data.matrix(simuE5I6[,-1])
is = isoSelect(y, X)
is

# mu = as.numeric(X[,is$w2use] %*% is$b2use)
# plot(mu, y)
# abline(0,1)

}
\keyword{ methods }
