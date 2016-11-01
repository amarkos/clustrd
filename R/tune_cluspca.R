tune_cluspca <- function(data, nclusrange = 2:7, ndimrange = 2:4, criterion = "asw", dst = "full", alpha = NULL, method = "RKM", center = TRUE, scale = TRUE, rotation = "none", nstart = 100, smartStart = NULL, seed = 1234){
  
  critval = matrix(0,max(length(nclusrange)),max(length(ndimrange)))
  
  m = 1
  n = 1
  for (k in nclusrange) {
    for (d in ndimrange) {
      #    if (k >= d) {
      ##    outclusCA[[k]] <- clusCA(data=data, nclus = k, ndim = d,nstart = nstart,smartStart = smartStart, seed = seed)
      print(paste('Running for',k,'clusters and',d,'dimensions...'))
      outcluspca <- cluspca(data=data, nclus = k, ndim = d,alpha = alpha,method = method,  center = TRUE, scale = TRUE, rotation=rotation, nstart = nstart, smartStart = smartStart, seed = seed)
      
      if (criterion == "asw")
      {
        critval[m,n] <- clusval(outcluspca, dst = dst)$asw
      }
      if (criterion == "ch")
      {
        critval[m,n] <- clusval(outcluspca, dst = dst)$ch
      }
      
      if (criterion == "crit")
      {
        critval[m,n] <- outcluspca$criterion
      }
      
      #    }
      n = n +1
    }
    n = 1
    m = m +1
  }
  
  #replace 0s with NAs
  critval[critval == 0] <- NA
  if (criterion == "crit")
  {
    indk.best <- which(critval == min(critval,na.rm =TRUE), arr.ind = TRUE)[1]
    indd.best <- which(critval == min(critval,na.rm =TRUE), arr.ind = TRUE)[2]
  } else {
    indk.best <- which(critval == max(critval,na.rm =TRUE), arr.ind = TRUE)[1]
    indd.best <- which(critval == max(critval,na.rm =TRUE), arr.ind = TRUE)[2]
  }
  
  k.best <- nclusrange[indk.best]
  d.best <- ndimrange[indd.best]
  outcluspcabest = cluspca(data = data, nclus = k.best, ndim = d.best,alpha = alpha, method = method,  center = TRUE, scale = TRUE, rotation= rotation, nstart = nstart, smartStart = smartStart, seed = seed)
  rownames(critval) = c(nclusrange)
  colnames(critval) = c(ndimrange)
  
  out <- list(cluspcaobj = outcluspcabest, nclusbest = k.best, ndimbest = d.best, critbest = round(critval[indk.best, indd.best],3), critgrid  = round(critval,3))
  out
}
