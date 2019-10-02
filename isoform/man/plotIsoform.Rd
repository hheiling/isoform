\name{plotIsoform}
\alias{plotIsoform}
\title{
  plot isoform information
}
\description{
  plot exon and exon junction expression levels and inferred
  isoforms
}
\usage{
  plotIsoform(gene, iso, cexLabel=0.9)
}
\arguments{
  \item{gene}{
  a list with two components, exon and junction
}
  \item{iso}{
  a matrix where the columns corresponding to isoforms and rows
  corresponding to exon and exon junctions specified in parameter
  \code{gene}
}
  \item{cexLabel}{
  font size of the exon or exon-junction labels
}

}
\references{
  Sun et al. IsoDOT Detects Differential RNA-isoform expression with respect to 
  a Categorical or Continuous Covariate with High Sensitivity and Specificity.
}
\author{ Wei Sun weisun@email.unc.edu }
\keyword{ methods }
