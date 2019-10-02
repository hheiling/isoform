effLen1 <-
function(rj, d, pdDist, lmax){
    
    eLen = 0
    
    if(rj >= d){
      lRange  = d:min(rj, lmax)
      eLen    = sum(pdDist[lRange]*(rj + 1 - lRange)) 
    }
    
    eLen
  }

