\name{effLen}
\alias{effLen}
\title{
	effective length of a series of exons
}
\description{
  calculate the effective length of a series of exons, allowing 
  exon-skipping and given the genomic lengths of all exons
}
\usage{
effLen(ids, rjs, d, pdDist, lmax)
}
\arguments{
  \item{ids}{
  the ids of the exons to be used. For example, ids=c(1,2,4) means 
  the exon set that includes the 1st, 2nd, and 4th exon. 
}
  \item{rjs}{
  the lengths of all the exons. 
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
  
For example, if there are 4 exons, with lengths of 100, 200, 300, and 
400 bps. Then rjs=c(100,200,300,400). If ids=c(1,2,4), we calculate the 
effective length of those reads that cover exons 1, 2, and 4, but does not 
cover exon 3. 
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
 \code{\link{effLen1}}, \code{\link{effLen2}}
}
\keyword{ utilities }
