\name{loadData}
\alias{loadData}
\title{
  load RNA-seq counts data
}
\description{
  load count data, which includes the number of RNA-seq fragments per exon 
  or exon set. 
}
\usage{
  loadData(countFile, bedFile, pdDist, readLen, lmax=500)
}
\arguments{
  \item{countFile}{
  Name of a text file recording the number of RNA-seq fragments per exon set, 
  which should  have 2 columns "count" and "exons", without header.
  For example:
  \verb{
 		110 ENSMUSG00000000001:1;
  		 16 ENSMUSG00000000001:1;ENSMUSG00000000001:2;
  }
}
  \item{bedFile}{
  Name of a bed file recording information of non-overlapping exons, which   
  have 6 columns: "chr", "start", "end", "exon", "score", and "strand", 
  without header. For example:
  \verb{
		chr1    3044314 3044814 ENSMUSG00000090025:1    666     +
		chr1    3092097 3092206 ENSMUSG00000064842:1    666     +
  }

}
  \item{pdDist}{
  the distribution of the lengths of the paired-end fragments. The i-th 
  component of this vector indicates the probability observing a 
  fragment of i bps long. 
}
  \item{readLen}{
  the length of a read, e.g., 75.
}
  \item{lmax}{
  maximum fragment length.
}
}
\value{
a list for all the genes. Each gene is a list of exon sets and the corresponding 
number of fragments and effective length. 
}
\references{
  Sun et al. IsoDOT Detects Differential RNA-isoform expression with respect to 
  a Categorical or Continuous Covariate with High Sensitivity and Specificity.
}
\author{ Wei Sun weisun@email.unc.edu }
\keyword{ methods }
