\name{compareAnno.cuffcompare}
\alias{compareAnno.cuffcompare}
\title{
  compare the results of compare and annotation at gene level
}
\description{
  for each gene from annotation, identify the closest gene from 
  cuffcompare output
}
\usage{
compareAnno.cuffcompare(cuffcompareGtf, gtf)
}
\arguments{
  \item{cuffcompareGtf}{
  the file name of cuffcompare output of gtf format. By cuffcompare's default, it 
  is "*.combined.gtf".
}
  \item{gtf}{
  gtf annotation file name, which should have the columns "chr", "source", 
  "feature", "start", "end", "score", "strand", "frame", and "anno".
}
}
\details{
 one gene is a set of exons (that may overlap). The distance of two genes are the 
 ratio of the number of bases that belong to only one gene divided by the number 
 of bases that belong to either gene. 
}
\value{
  A data.frame of 3 columns:
  \item{geneId }{gene ID from annotation}
  \item{cuffId }{ gene ID of the cuffcompare gene that is closest to the annotation gene `geneId'}
  \item{dist }{distance between the annotation gene `geneId' and the cuffcompare gene `cuffId'.}
}
\author{ Wei Sun weisun@email.unc.edu }
\keyword{ utilities }
