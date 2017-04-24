summary.clusmca <- function(object, ...) {
  
  x = object 
  k = length(x$size)
  d = dim(x$centroid)[2]
  size = x$size
  csize = round((table(x$cluID)/sum(table(x$cluID)))*100,digits=1)
  tt = paste('(',csize,'%)',sep="")
  cs = paste(size, tt, sep = " ", collapse = ", ")
  
  cat(paste("Solution with ",k ," clusters of sizes ", paste(cs, collapse = ", ")," in ",d ," dimensions. ","\n", sep = ""))
  
  cat("\nCluster centroids:\n")
  xcent = data.frame(x$centroid)
  for (i in 1:k) {
    rownames(xcent)[i] = paste("Cluster",i)
  }
  for (i in 1:ncol(xcent)) {
    colnames(xcent)[i] = paste("Component",i)
  }
  print(xcent)
  attc = data.frame(x$attcoord)
  cat("\nVariable scores:\n")
  for (i in 1:ncol(attc)) {
    colnames(attc)[i] = paste("Component",i)
  }
  print(attc)
  
  cat("\nClustering vector:\n")
  print(x$cluID)
  
  cat(paste("\nObjective criterion value:",round(x$criterion,3),"\n"))
  
  cat("\nAvailable output:\n", 
      sep = "\n")
  print(names(x))
  invisible(x)
}