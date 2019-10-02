distSet <-
function(set1, set2){
  
  pos1 = min(set1$start, set2$start) - 1
  pos2 = max(set1$end,   set2$end)
  
  set1$start = set1$start - pos1 
  set2$start = set2$start - pos1
  
  set1$end   = set1$end - pos1 
  set2$end   = set2$end - pos1 
  
  ind1 = ind2 = rep(0, pos2 - pos1)
  
  for(i in 1:nrow(set1)){
    ind1[set1$start[i]:set1$end[i]] = 1
  }
  
  for(i in 1:nrow(set2)){
    ind2[set2$start[i]:set2$end[i]] = 1
  }
    
  length(which(ind1 != ind2))/length(which(ind1 + ind2 > 0))
  
}

