\name{effLen1}
\alias{effLen1}
\title{
 effective length of a single exon
}
\description{
 calculate the effective length of a single exon
}
\usage{
effLen1(rj, d, pdDist, lmax)
}
\arguments{
  \item{rj}{
  the length (in bps) of the exon of interest. 
}
  \item{d}{
  size of a read at one end of a paired-end fragment
}
  \item{pdDist}{
  the distribution of the lengths of the paired-end fragments
}
  \item{lmax}{
  maximum RNA-seq fragment length that includes two reads at both ends
}
}
\details{
Intuitively, effective length is the usable length for the data 
generation mechanism, i.e., the number of positions where a randomly 
selected RNA-seq fragment can be sampled.
}
\value{
a numeric value of effective length.
}
\references{
  Sun et al. IsoDOT Detects Differential RNA-isoform expression with respect to 
  a Categorical or Continuous Covariate with High Sensitivity and Specificity.
}
\author{ Wei Sun weisun@email.unc.edu }
\seealso{
 \code{\link{effLen}}, \code{\link{effLen2}}
}
\keyword{ utilities }
