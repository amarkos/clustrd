tune_clusmca <- function(data, nclusrange = 2:7, ndimrange = 2:4, method = "clusCA", criterion = "asw", dst = "full", alpha = .5, nstart = 10, smartStart = NULL, seed = 1234){
  
 # outclusmca = list()
  critval = matrix(0,max(length(nclusrange)),max(length(ndimrange)))
  
  m = 1
  n = 1
 # if ((method == "clusCA") | (method == "iFCB")) {    #for Cluster CA or iFCB return missing when k <= d
    for (k in nclusrange) {
      for (d in ndimrange) {
        if (k > d) {
          ##    outclusCA[[k]] <- clusCA(data=data, nclus = k, ndim = d,nstart = nstart,smartStart = smartStart, seed = seed)
          print(paste('Running for',k,'clusters and',d,'dimensions...'))
          outclusmca <- clusmca(data = data, nclus = k, ndim = d,method = method, alpha = alpha, nstart = nstart,smartStart = smartStart, seed = seed)
          
          if (criterion == "asw")
          {
            critval[m,n] <- clusval(outclusmca, dst = dst)$asw
          }
          if (criterion == "ch")
          {
            critval[m,n] <- clusval(outclusmca, dst = dst)$ch
          }
          
          if (criterion == "crit")
          {
            critval[m,n] <- outclusmca$criterion
          }
          
        }
        n = n +1
      }
      n = 1
      m = m +1
    }
  # } else { #for MCAk return everything
  #   for (k in nclusrange) {
  #     for (d in ndimrange) {
  #       ##    outclusCA[[k]] <- clusCA(data=data, nclus = k, ndim = d,nstart = nstart,smartStart = smartStart, seed = seed)
  #       print(paste('Running for',k,'clusters and',d,'dimensions...'))
  #       outclusmca <- clusmca(data = data, nclus = k, ndim = d,method = method, alpha = alpha, nstart = nstart,smartStart = smartStart, seed = seed)
  #       
  #       if (criterion == "asw")
  #       {
  #         critval[m,n] <- clusval(outclusmca, dst = dst)$asw
  #       }
  #       if (criterion == "ch")
  #       {
  #         critval[m,n] <- clusval(outclusmca, dst = dst)$ch
  #       }
  #       
  #       if (criterion == "crit")
  #       {
  #         critval[m,n] <- outclusmca$criterion
  #       }
  #       n = n +1
  #     }
  #     n = 1
  #     m = m +1
  #   }
  # }
  
  #replace 0s with NAs
  critval[critval == 0] <- NA
  
  if (criterion == "crit")
  {
    indk.best <- which(critval == min(critval,na.rm =TRUE), arr.ind = TRUE)[1]
    indd.best <- which(critval == min(critval,na.rm =TRUE), arr.ind = TRUE)[2]
    #FIX: in case of tie returns the lowest (more parsimonious)
    #  print(indk.best)
    #print(indd.best)
  } else {
    indk.best <- which(critval == max(critval,na.rm =TRUE), arr.ind = TRUE)[1]
    indd.best <- which(critval == max(critval,na.rm =TRUE), arr.ind = TRUE)[2]
  }
  k.best <- nclusrange[indk.best]
  d.best <- ndimrange[indd.best]
  
  outclusmcabest = clusmca(data = data, nclus = k.best, ndim = d.best,method = method, alpha = alpha, nstart = nstart,smartStart = smartStart, seed = seed)
  rownames(critval) = c(nclusrange)
  colnames(critval) = c(ndimrange)
  
  crit.best = round(critval[indk.best, indd.best],3) 
  crit.grid  = round(critval,3)
  
  crit.grid[upper.tri(as.matrix(crit.grid),diag=TRUE)] <- "   "
  crit.grid = as.data.frame(crit.grid)
  out <- list(clusmcaobj = outclusmcabest, nclusbest = k.best, ndimbest = d.best, critbest = crit.best, critgrid  = crit.grid)
  out
}
