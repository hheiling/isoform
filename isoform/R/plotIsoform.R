plotIsoform <-
function(gene, iso, cexLabel=0.9)
{
  
  # --------------------------------------------------------- 
  # extract information
  # --------------------------------------------------------- 
  
  exon  = gene$exon
  nExn  = nrow(exon)
  
  loc   = cbind(exon$exon_id, exon[,3:5])
  names(loc) = c("exon_id1", "exon_id2", "count", "eLen")
  loc   = rbind(loc, gene$junc)
  
  par(mfrow=c(1,2), mar=c(1,1,1,1), xaxt="n", yaxt="n", bty="n")

  # --------------------------------------------------------- 
  # plot exons
  # --------------------------------------------------------- 

  starts = exon$start
  ends   = exon$end

  pos1   = starts
  pos2   = ends

  intronShrink = 1
  nExon = nrow(exon)

  if(nExon > 1){
    for(i in 2:nExon){
      pos1[i] = pos2[i-1] + intronShrink*(starts[i] - ends[i-1])
      pos2[i] = pos1[i] + (ends[i] - starts[i] + 1)
    }
  }

  plot(c(min(pos1), max(pos2)), c(-ncol(iso), 3.5), type="n", xlab="", ylab="")
  rect(pos1, 2.5, pos2, 3, col="brown", border=NA)
  text(min(pos1), 3.5, "[gene model]", pos=4)

  rate = exon$count/exon$eLen
  rect(pos1, 0.5, pos2, rate/max(rate)+0.5, col="seagreen", border=NA)
  text(min(pos1), 2, "[relative expression]", pos=4)

  # --------------------------------------------------------- 
  # plot detected isoforms
  # --------------------------------------------------------- 

  text(min(pos1), 0, "[possible isoforms]", pos=4)

  for(k in 1:ncol(iso)){
    wiso = which(iso[,k] == 1)
    rect(pos1[wiso], -k, pos2[wiso], -k+0.5, col="skyblue", border=NA)
  }

  # --------------------------------------------------------- 
  # plot exon junctions
  # --------------------------------------------------------- 

  mdd = matrix(0, nrow=nExn, ncol=nExn)
  for(i in 1:nrow(loc)){
    mdd[loc[i,1], loc[i,2]] = loc$count[i]/loc$eLen[i]
  }
  mdd = mdd/max(mdd)

  col = c("white", "seagreen3", "darkgreen", "orange", "red")
  bks = c(0, seq(0,1,length.out=5))

  image(mdd, col=col, breaks=bks) 

  abline(h=seq(-0.5, nExn-0.5, by=1)/(nExn-1), lty=1, col="white")
  abline(v=seq(-0.5, nExn-0.5, by=1)/(nExn-1), lty=1, col="white")

  text(x=(0:(nExn-1))/(nExn-1), y=(0:(nExn-1))/(nExn-1), 
  labels=1:nExn, col="lightyellow", font=2, cex=cexLabel)

  for(i in 1:(nrow(mdd)-1)){
    for(j in (i+1):nrow(mdd)){
      if(mdd[i,j] > 0){
        text(x=(i-1)/(nExn-1), y=(j-1)/(nExn-1), font=2, cex=cexLabel, 
             labels=sprintf("%d-%d", i, j), col="lightyellow")
      }
    }
  }

  lg = c("0.00-0.25", "0.25-0.50", "0.50-0.75", "0.75-1.00")
  legend("bottomright", legend=lg, fill=col[-1], bty="n")
}

