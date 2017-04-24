## Define a print method that will be automatically dispatched when print()
## is called on an object of class "tune_clusmca"
print.tune_clusmca <- function(x, ...) {
  
  size = x$clusobj$size
  csize = round((table(x$clusobj$cluID)/sum(table(x$clusobj$cluID)))*100,digits=1)
  k = x$nclusbest
  d = x$ndimbest

  tt = paste('(',csize,'%)',sep="")
  cs = paste(size, tt, sep = " ", collapse = ", ")
  cat(paste("\nThe best solution was obtained for ",k ," clusters of sizes ", paste(cs, collapse = ", ")," in ",d ," dimensions, for a cluster quality criterion value of ",round(x$critbest,3), ".", "\n", sep = ""))
  
  cat("\nCluster quality criterion values across the specified range of clusters (rows) and dimensions (columns):\n")
  print(x$critgrid)
  
  cat("\nCluster centroids:\n")
  xcent = data.frame(x$clusobj$centroid)
  for (i in 1:k) {
    rownames(xcent)[i] = paste("Cluster",i)
  }
  for (i in 1:ncol(xcent)) {
    colnames(xcent)[i] = paste("Component",i)
  }
  print(xcent)
  attc = data.frame(x$clusobj$attcoord)
  cat("\nVariable scores:\n")
  for (i in 1:ncol(attc)) {
    colnames(attc)[i] = paste("Component",i)
  }
  print(attc)

  cat(paste("\nObjective criterion value:",round(x$clusobj$criterion,3),"\n"))
  
  cat("\nAvailable output:\n", 
      sep = "\n")
  print(names(x))
  invisible(x)
  
}