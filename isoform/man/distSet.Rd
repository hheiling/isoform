\name{distSet}
\alias{distSet}
\title{
  calculate the distance of two exon sets
}
\description{
each exon set can be an isoform. Thus this function calculates the distance of 
two isoforms
}
\usage{
distSet(set1, set2)
}
\arguments{
  \item{set1}{
  a data frame of k rows and 2 columns, where k is the number of exons and 
  the two columns have names "start", "end". 
}
  \item{set2}{
  a data frame of k rows and 2 columns, where k is the number of exons and 
  the two columns have names "start", "end". 
}
}
\value{
a ratio of the number of base pairs belonging only one of set1 and set2 over the 
number of base pairs belonging to either set1 or set2
}
\references{
  Sun et al. IsoDOT Detects Differential RNA-isoform expression with respect to 
  a Categorical or Continuous Covariate with High Sensitivity and Specificity.
}
\author{ Wei Sun weisun@email.unc.edu }
\keyword{ utilities }
