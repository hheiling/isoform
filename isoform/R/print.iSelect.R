`print.iSelect` <-
function(x, ...)
{

  ww = which(! names(x) %in% c("score", "mu", "pZI"))
  
  for(i in ww){
    cat(names(x)[i], ":\n", sep="")
    print(x[[i]])
    cat("\n")
  }
  
}

