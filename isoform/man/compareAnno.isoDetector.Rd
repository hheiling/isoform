\name{compareAnno.isoDetector}
\alias{compareAnno.isoDetector}
\title{
  compare the results of isoDetector and annotation
}
\description{
  Compare the results of isoDetector and annotation gene by gene. For each gene, 
  	check each of its isoforms identified by isoDetector to identify the 
  	closest isoform from the annotation. 
}
\usage{
compareAnno.isoDetector(geneMod, gtf)
}
\arguments{
  \item{geneMod}{
	a list, the result of function geneModel
}
  \item{gtf}{
  gtf annotation file name, which should have the columns "chr", "source", 
  "feature", "start", "end", "score", "strand", "frame", and "anno".
}
}

\value{
  A data.frame of 5 columns:
  \item{tranId }{transcript ID, obtained from gtf annotation}
  \item{clustId }{cluster ID, obtained from nTE}
  \item{dist }{distance between the isoDetector transcript and the closest transcript from annotation}
  \item{len }{length of this isoform identified by isoDetector}
  \item{abundance }{abundance (number of fragments) of the isoform identified by isoDetector}
  \item{status}{status from isoDetector, could be "OK" or "LOWDATA". }
}
\author{ Wei Sun weisun@email.unc.edu }
\keyword{ utilities }
