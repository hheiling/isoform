countReads <-
function(bamFile, bedFile, outputFile, overlapFraction=0.9, forceStrand=FALSE)
{
  if(file.access(bamFile, mode = 0) != 0){
    stop(sprintf("bamFile %s does not exist\n", bamFile))
  }
  
  if(file.access(bedFile, mode = 0) != 0){
    stop(sprintf("bedFile %s does not exist\n", bedFile))
  }
  
  tmpFile  = sprintf("%s_mapped_2Exon.txt", bamFile)

  Z = .C("countReads", as.character(bamFile), as.character(bedFile), 
         as.character(tmpFile), as.double(overlapFraction), 
         as.integer(forceStrand), PACKAGE="isoform")
  
  cmd = sprintf("cut -f 2 %s | sort | uniq -c > %s", tmpFile, outputFile)
  system(cmd)
  
  system(sprintf("rm -rf %s", tmpFile))
}
