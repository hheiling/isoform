\name{knownIsoforms}
\alias{knownIsoforms}
\title{
extract known isoform information into a R 'list' object.
}
\description{
extract known isoform information from a gtf file into a R 'list' object.
The gtf file should have an annotation of non-overlapping exons and the 
last column of this file should have information like:
gene_id "ENSMUSG00000090025"; transcript_id "ENSMUST00000160944"; 
gene_name "Gm16088"; transcript_name "Gm16088-001"; exon_id "1";
}
\usage{
knownIsoforms(nonOverlapExonFile)
}
\arguments{
  \item{nonOverlapExonFile}{
input file name
}
}
\value{
a list which each element for one gene.
}
\references{
  Sun et al. IsoDOT Detects Differential RNA-isoform expression with respect to 
  a Categorical or Continuous Covariate with High Sensitivity and Specificity.
}
\author{ Wei Sun weisun@email.unc.edu }
\keyword{ methods }
