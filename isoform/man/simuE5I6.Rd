\name{simuE5I6}
\alias{simuE5I6}
\docType{data}
\title{
a simulated dataset of a gene with 5 isoform and 6 possible isoforms.
}
\description{
the data are simulated such as only the first three isoforms are expressed. 
}
\usage{data(simuE5I6)}
\format{
  A data frame with 15 observations on the following 7 variables.
  \describe{
    \item{\code{y}}{ the number of RNA-seq reads per exon set. }
    \item{\code{x1,x2,x3,x4,x5,x6}}{ the effective length of each exon set within the 6 isoforms. }
  }
}
\details{
  the data are simulated such as only the first three isoforms are expressed.
}

\examples{
data(simuE5I6)

y  = simuE5I6$y
X  = data.matrix(simuE5I6[,-1])
is = isoSelect(y, X)
is

}
\keyword{datasets}
