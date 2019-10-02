\name{isoSelectZI}
\alias{isoSelectZI}
\title{
  isoSelect plus zero-inflated component
}
\description{
  select isoform subset by penalized generalized regression, while extending 
  isoSelect by employing zero-inflated negative binomial distribution rather than 
  negative binomial distribution used by isoSelect.  
}
\usage{
isoSelectZI(y, X, muMin, useBIC, nlambda, ntau, maxIt = 5, 
	min_pi_ZI = 0.001, trace = 0)
}

\arguments{
  \item{y}{
	a vector of the response variable. 
	}
	
  \item{X}{
	a matrix of covariate data, each column is one covariate and 
	each row is one observation. Intercept should not be included in 
	this covariate matrix. 
	}
	
  \item{muMin}{
  	mu is the expected value of read count for an exon set, which is modelled 
  	by a negative binomial distribution. if mu=0, we set mu=muMin. This is used 
  	during GLM fitting to improve stability of the algorithm. It should be set to 
  	be the same value as is used by function geneModel.
	}
	
  \item{useBIC}{
  	useBIC may be specified as "BIC" or "extended BIC", which are two criterion for 
  	tuning parameter selection. extended BIC should be used if the number of 
  	candidate isoform is close to or larger than sample size (i.e., the number 
  	of exon sets. )
	}
	
  \item{nlambda}{
	number of values of tuning parameter lambda. }
	
  \item{ntau}{
	number of values of tuning parameter tau. }

  \item{maxIt}{
    the maximum number of iterations to estimate zero-inflated component
	}
	
  \item{min_pi_ZI}{
    the minimum value of the mixture proportion for ZI component in the mixture model. 
	}
	
  \item{trace}{
    to print out tracing report or not
    }
    
}

\value{
  a list, with the same elements as the output of isoSelect and two extra elements: 
  \item{pZI }{posterior probability for each observation arising from the 
  	zero-inflated component}

  \item{mu }{expected value for each observation}
}

\author{
Wei Sun weisun@email.unc.edu
}
\references{
  Sun et al. IsoDOT Detects Differential RNA-isoform expression with respect to 
  a Categorical or Continuous Covariate with High Sensitivity and Specificity.
}

\seealso{
\code{\link{isoSelect}}, \code{\link{isoDu}}
}
\keyword{ methods }
