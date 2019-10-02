\name{effLen2}
\alias{effLen2}
\title{
 effective length of two non-adjancent exons
}
\description{
 calcualte the effective length of two non-adjancent exons
}
\usage{
effLen2(rj, rh, rk, d, pdDist, lmax)
}
\arguments{
  \item{rj}{
  the length (in bps) of the first exon of interest.
}
  \item{rh}{
  the length (in bps) of the second exon of interest (the skipped exon).
}
  \item{rk}{
  the length (in bps) of the third exon of interest.
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
selected RNA-seq fragment can be sampled. This function calculate 
the effective length of the reads that overlap both first and third 
exons (of length rj and rk bps), but skip the second exon of rh bps. 
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