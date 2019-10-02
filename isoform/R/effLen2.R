effLen2 <-
function(rj, rh, rk, d, pdDist, lmax){
  eLen = 0
  
  if(rj < d || rk < d || rh+2*d > lmax){
    eLen   = 0
  }else{
    lRange  = (2*d+rh):min(rj+rh+rk, lmax)
    ljRange = pmin(rj, lRange - rh - d) - pmax(d, lRange - rh - rk) + 1
    eLen    = sum(pdDist[lRange]*(ljRange )) 
  }
}

