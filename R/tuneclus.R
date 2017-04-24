tuneclus <- function(data, nclusrange = 2:7, ndimrange = 2:4, method = "RKM", criterion = "asw", dst = "full", alpha = NULL, alphak = NULL, center = TRUE, scale = TRUE, rotation = "none", nstart = 100, smartStart = NULL, seed = 1234){
  #wrapper for functions tune_cluspca(), tune_clusmca()
  
  method <- match.arg(method, c("RKM", "rkm","rKM","FKM", "fkm","fKM","clusCA", "clusca","CLUSCA","CLUSca", "ifcb","iFCB","IFCB","mcak", "MCAk", "MCAK","mcaK"), several.ok = T)[1]
  method <- tolower(method)
  
  if (method %in% c("rkm","fkm")) {
    out = tune_cluspca(data, nclusrange, ndimrange, method = method, criterion = criterion, dst = dst, alpha = alpha, center = center, scale = scale, rotation = rotation, nstart = nstart, smartStart = smartStart, seed = seed)
  } else if (method %in% c("clusca","ifcb","mcak")) {
    out = tune_clusmca(data, nclusrange, ndimrange, method = method, criterion = criterion, dst = dst,  alphak = alphak, nstart = nstart, smartStart = smartStart, seed = seed)
  }
  
  out
}