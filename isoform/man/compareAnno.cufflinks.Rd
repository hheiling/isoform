\name{compareAnno.cufflinks}
\alias{compareAnno.cufflinks}
\title{
  compare the results of cufflinks and annotation at transcript level
}
\description{
  compare the results of cufflinks and annotation. For each isoform identified by 
  cufflinks, identify the closest isoform from the annotation. 
}
\usage{
compareAnno.cufflinks(cufflinksGtf, isoTracking, gtf, knownIso)
}
\arguments{
  \item{cufflinksGtf}{
  the file name of cufflinks output of gtf format. By cuffllinks's default, 
  it is "transcripts.gtf".
}
  \item{isoTracking}{
  the file name of cufflinks output of isoform tracking. By cuffllinks's default, 
  it is "isoforms.fpkm_tracking".
}

  \item{gtf}{
  gtf annotation file name, which should have the columns "chr", "source", 
  "feature", "start", "end", "score", "strand", "frame", and "anno".
}
  \item{knownIso}{
  if knowIso=TRUE, cufflinks was used to estimate expresion on a know set of 
  isoforms instead of detecting de novo isoforms. 
}
}
\value{
  A data.frame of 10 columns:
  \item{cuffTranId }{cufflinks transcript ID}
  \item{cuffGeneId }{cufflinks gene ID}
  \item{tranId }{transcript ID of the transcript (from annotation) that is closest 
  	to the cufflinks transcript}
  \item{geneId }{the corresponding gene ID of the transcript "tranId"}
  \item{dist }{distance between the cufflinks transcript and the closest 
  	transcript from annotation}
  \item{len }{length of the cufflinks transcript}
  \item{FPKM }{FPKM estimated from cufflinks}
  \item{confHi }{Upper bound of the confidence of FPKM}
  \item{confLo }{Lower bound of the confidence of FPKM}
  \item{status }{status of the isoform estimation from cufflinks}
}
\author{ Wei Sun weisun@email.unc.edu }
\keyword{ utilities }
