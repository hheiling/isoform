\name{isoDuEx2}
\alias{isoDuEx2}
\docType{data}
\title{
  a testing dataset for isoDu (example 2)
}
\description{
 this testing dataset is build using function geneModel, 
 using RNAseq data from 10 HapMap CEU samples
}
\usage{data(isoDuEx2)}
\format{
	a list with four elements
  \describe{
    \item{\code{readDepth}}{ read depth across the 10 samples. }
    \item{\code{xData}}{ genotype data of 8 SNPs, which are around a transcription cluster }
    \item{\code{pdDistL}}{ distributino of RNA-seq fragment length. }
    \item{\code{gL}}{ gene model list }
  }
  }
\details{
 only two transcription clusters are included to reduce file size
}
\source{
  Montgomery et al. Transcriptome genetics using second generation sequencing in a Caucasian population. Nature. 2010 Apr 1;464(7289):773-7. }
\keyword{datasets}
