\name{readFragmentSizes}
\alias{readFragmentSizes}
\title{
read in the distribution of RNA-seq fragment sizes.
}
\description{
The distribution of RNA-seq fragment sizes per sample is saved in a text file 
and this function read in such information for a number of samples. 
}
\usage{
readFragmentSizes(fragSizeFiles, lmax=500)
}
\arguments{
  \item{fragSizeFiles}{
  names of the text files recording the distribution of the fragment lengths. 
  The file names should be one-to-one correspondence with the samples 
  specified in tags. Each file should have 2 columns: 
  "Frequency" and "Length", without header. 
  For example:
  \verb{
	20546 75
	40465 76
	37486 77
	27533 78
	25344 79
	}
}
  \item{lmax}{
  maximum fragment length.
}
}
\value{
  a list where each element includes the fragment length distribution for a sample
}
\references{
  Sun et al. IsoDOT Detects Differential RNA-isoform expression with respect to a 
  Categorical or Continuous Covariate with High Sensitivity and Specificity.
}
\author{ Wei Sun weisun@email.unc.edu }
\seealso{
 \code{\link{isoDu}}
}
\keyword{ utilities }
