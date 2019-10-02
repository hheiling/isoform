\name{isoDetector}
\alias{isoDetector}
\title{
  RNA isoform detector
}
\description{
 RNA isform selection and abundance estimation
}
\usage{
isoDetector(countFile, bedFile, fragSizeFile, readLen, output, lmax = 500, 
			minCount = 5, minObs = 2, knownIsoforms = NULL, maxBreaks = 5, 
			pvalBreaks = 0.05, pvalExpress = 0.01, foldExpress = 1/5, 
			eLenMin = 1, muMin = 0.01, pMaxRel = 10, pMaxAbs = 2000, 
			verbose = 1)
}
\arguments{
  \item{countFile}{
  Name of a text file recording the number of RNA-seq fragments per exon set, 
  which should  have 2 colums "count" and "exons", without header.
  For example:
  \verb{
 		110 ENSMUSG00000000001:1;
  		 16 ENSMUSG00000000001:1;ENSMUSG00000000001:2;
  }
}
  \item{bedFile}{
  Name of a bed file recording information of non-overlapping exons, which   
  have 6 colums: "chr", "start", "end", "exon", "score", and "strand", 
  without header. For example:
  \verb{
		chr1    3044314 3044814 ENSMUSG00000090025:1    666     +
		chr1    3092097 3092206 ENSMUSG00000064842:1    666     +
  }

}
  \item{fragSizeFile}{
  Name of a text file recording the distribution of the fragment lengths, 
  which should have 2 columns: "Frequency" and "Length", without header
  For example:
  \verb{
	20546 75
	40465 76
	37486 77
	27533 78
	25344 79
	}
}
  \item{readLen}{
  the length of a read, e.g., 75.
}
  \item{output}{
  the name of an output file, which is an RData file. 
}
  \item{lmax}{
  maximum fragment length.
}
  \item{minCount}{
  we say a transcription cluster is non-trivial for the study of RNA isoforms
  if it has at least \code{minCount} exon sets with \code{minObs} or 
  more fragments. 
}
  \item{minObs}{
  see the description for minCount
}
  \item{knownIsoforms}{
  If knownIsoforms is NULL, candidate isoforms will be selected from all 
  possible isoforms; otherwise, knownIsoforms is the name of an .RData file
  which has a list named `isoAll'. Each component of isoAll is a binary matrix 
  that specifies a set of possible isoforms (e.g., isoforms from 
  annotations). Specifically, it is a binary matrix 
  of k rows and m columns, where k is the number of 
  non-overlapping exons and m is the number of isoforms. isoforms[i,j]=1 
  indicates that the i-th exon belongs to the j-th isoform. For example, 
  the following matrix indicates the three isoforms for one gene 
  ENSMUSG00000000003:
  \verb{
       ENSMUST00000000003 ENSMUST00000166366 ENSMUST00000114041
[1,]                  1                  1                  1
[2,]                  1                  1                  1
[3,]                  1                  1                  1
[4,]                  1                  1                  0
[5,]                  1                  1                  1
[6,]                  1                  1                  1
[7,]                  1                  1                  1
[8,]                  1                  0                  0
  }
}
  \item{maxBreaks}{
  the maximum number of break points. We detect possible boundaries of 
  RNA isoforms by checking whether the read counts of adjacent exons 
  are consistent. maxBreaks indicates the maximum number of such boundaries. 
}
  \item{pvalBreaks}{
  We detect possible boundaries of RNA isoforms by checking whether the 
  read counts of adjacent exons are consistent. A Chis-square p-value is 
  calculated to assess the consistency. We only use the breaks with p-value
  smaller than pvalBreaks. If there are more than maxBreaks p-values smaller 
  than pvalBreaks, we choose the smallest maxBreaks break points. 
}
  \item{pvalExpress}{
  We test whether one exon set has less reads than expected by chance. 
  Smaller p-value implies the corresponding exon set may not be part of 
  the transcriptome or it belongs to a lowly expressed transcript. We add a 
  particular exon set with at least one skipped exon into the RNA 
  isoforms if the corresponding p-value is larger than pvalExpress. 
}
  \item{foldExpress}{
  We also add a particular exon set with at least one skipped exon 
  into the RNA isoforms if the ratio of count vs. effective length is larger 
  than foldExpress multiplies the average raio across the gene. 
}
  \item{eLenMin}{
  if effective length of an exon or exon junction is smaller than eLenMin,
  i.e., if this exon is not included in the corresponding isoform, 
  set it to eLenMin. This is to account for possible sequencing error or
  mapping errors.
}
  \item{muMin}{
  if mu=0, set it to muMin. This is used during GLM fitting for stability of the 
  algorithm.
}
  \item{pMaxRel}{
  quit the computation if the number of possible isoforms over the number of 
  exon sets is larger than pMaxRel.
}
  \item{pMaxAbs}{
  quit the computation if the total number of possible isoforms is larger 
  than pMaxRel.
}
  \item{verbose}{
  if verbose = 1, print out basic running information. If verbose>1, 
  print out more details e.g., the dimension of the design matrix 
  for each transcription cluster
}
}
\value{
  A list with each component for a transcription cluster. Then each 
  component is a list with following elements:
  
    \item{info}{ a matrix of the location information of each exon.}
    
    \item{count}{ a matrix of the effective length, and sequence counts of
    each exon set. Here effective length is calculated assuming all 
    exons of the exon set are consecutive. }
    
    \item{y}{ a vector of sequence counts of each exon set. }

    \item{candiIsoform}{ a matrix of all the candidate isoforms. Each row
    corresponds to an exon and each column corresponds to an isoform. }

    \item{X}{ The design matrix of all the candidate isoforms. Each row
    corresponds to an exon set (matching row in matrix count), and 
    each column match to the isoform specified in candiIsoform. }

    \item{w2kp}{ indices of the isoforms selected by penalized regression. }
    
	\item{status}{ "FAIL", "LOWDATA", or "OK". }

    \item{abundance}{ estimated number of fragments of each selected isoform. }
    
    \item{fitted}{ a vector of fitted value for each exon set 
    based on the penalized regression. }
    
    \item{likelihood}{ a vector likelihood for each exon set 
    based on the penalized regression. }

    \item{family}{ either `Poisson' or `Negative Binomial'. }

    \item{phi}{ the overdispersion parameter for `Negative Binomial' 
    distribution, and it is 0 if the distribution is `Poisson'. }

}
\references{
  Sun et al. IsoDOT Detects Differential RNA-isoform expression with respect to 
  a Categorical or Continuous Covariate with High Sensitivity and Specificity.
}
\author{ Wei Sun weisun@email.unc.edu }
\keyword{ methods }

