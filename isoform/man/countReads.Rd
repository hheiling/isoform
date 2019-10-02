\name{countReads}
\alias{countReads}
\title{
  count the number of RNA-seq reads
}
\description{
  count the number of RNA-seq reads from a bam file that are mapped to features
  described in a bed file. 
}
\usage{
countReads(bamFile, bedFile, outputFile, overlapFraction = 0.9, 
	forceStrand = FALSE)
}
\arguments{
  \item{bamFile}{
  a bam file that stores the mapped RNA seq reads
}
  \item{bedFile}{
  a bed file that stores the annotations of the genomic features, for example, 
  the information of each exon. 
}

  \item{outputFile}{
  the name of the output file. 
}

  \item{overlapFraction}{
  minimum overlap required as a fraction of a RNA-seq read. 
  The default is 0.9. 
}
  \item{forceStrand}{
  force strandedness.  That is, only report hits in bed file that overlaps 
  with a read in bam file on the same strand.
}
}
\value{
  the results will be written into the output file. 
}
\references{
  Sun et al. IsoDOT Detects Differential RNA-isoform expression with respect to 
  a Categorical or Continuous Covariate with High Sensitivity and Specificity.
}
\author{ Wei Sun weisun@email.unc.edu }
\keyword{ methods }

