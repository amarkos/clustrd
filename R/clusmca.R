clusmca <- function(data,nclus,ndim,method="clusCA",alphak=.5,nstart=100,smartStart=NULL,gamma = TRUE,seed=1234){

  method <- match.arg(method, c("clusCA", "clusca","CLUSCA","CLUSca", "ifcb","iFCB","IFCB","mcak", "MCAk", "MCAK","mcaK"), several.ok = T)[1]
  method <- tolower(method)
  
  if(method=="clusca"){
    out=clusCA(data=data,nclus=nclus,ndim=ndim,nstart=nstart,smartStart=smartStart, gamma = gamma,seed=seed)
  }
  if(method=="ifcb"){
    out=iFCB(data=data,nclus=nclus,ndim=ndim,nstart=nstart,smartStart=smartStart, gamma = gamma,seed=seed)
  }
  if(method=="mcak"){
    out=MCAk(data=data,nclus=nclus,ndim=ndim,nstart=nstart,alphak = alphak,smartStart=smartStart, gamma = gamma,seed=seed)
  }
  return(out)
  }



