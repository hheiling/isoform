\name{geneModel}
\alias{geneModel}
\title{
  Build a gene model
}
\description{
  construct RNA isoforms, and estimate their abundances
}
\usage{
geneModel(gene, d, pdDist, isoforms=NULL, lmax = length(pdDist), 
		 maxBreaks=5, pvalBreaks=0.05, pvalExpress=0.01, 
		 foldExpress=1/5, isoformOnly=FALSE, runIsoSelect = TRUE,  
		 eLenMin=1.0, muMin=0.01, pMaxRel=10, pMaxAbs=2000, verbose=1)
}
\arguments{
  \item{gene}{
  a list with two components, info and count
}
  \item{d}{
  size of a read at one end of a paired-end fragment
}
  \item{pdDist}{
  the distribution of the lengths of the paired-end fragments. The i-th 
  component of this vector indicates the probability observing a 
  fragment of i bps long. 
}
  \item{isoforms}{
  If isoforms is NULL, candidate isoforms will be selected from all 
  possible isoforms; otherwise, isoforms should be a binary matrix 
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
  \item{lmax}{
  maximum RNA-seq fragment length that includes two reads at both ends
  (or insert size)
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
  \item{isoformOnly}{
  if isoformOnly is TRUE, return the selected isoform, do not calculate effective
  length and do not apply penalized regression to further select a subset of 
  isoforms.
}
  \item{runIsoSelect}{
  if runIsoSelect is TRUE, run function isoSelect with default setting to select
  isoforms. Otherwise geneModel return a list with y and X and user can run 
  isoSelect in a customized way.
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

  If isoformOnly is TRUE, only return the isoform matrix. 

  Otherwise if runIsoSelect is FALSE, return a list with the following components:
  
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

  If runIsoSelect is TRUE, return a list with all the above components, 
  plus the following components:

    \item{w2kp}{ indices of the isoforms selected by penalized regression. }
    
	\item{status}{ "FAIL", "LOWDATA", or "OK". }

    \item{abundance}{ estimated number of fragments of each selected isoform. }
    
    \item{fitted}{ a vector of fitted value for each exon set 
    based on the penalized regression. }
    
    \item{likelihood}{ a vector likelihood for each exon set 
    based on the penalized regression. }

    \item{family}{ either `Poisson' or `Negative Binomial'. }

    \item{phi}{ the overdisperson parameter for `Negative Binomial' 
    distribution, and it is 0 if the distribtion is `Poisson'. }

}
\references{
  Sun et al. IsoDOT Detects Differential RNA-isoform expression with respect to 
  a Categorical or Continuous Covariate with High Sensitivity and Specificity.
}
\author{ Wei Sun weisun@email.unc.edu }
\keyword{ methods }
