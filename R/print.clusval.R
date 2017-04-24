print.clusval <- function(x, ...) {
  cat(paste("\nAverage silhouette width:", round(x$asw,digits=3)))
  cat(paste("\nCalinski-Harabasz:", round(x$ch,digits=3)))  
}