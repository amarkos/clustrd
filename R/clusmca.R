clusmca <- function(data,nclus,ndim,method="clusCA",alpha=.5,nstart=10,smartStart=NULL,gamma = TRUE,seed=1234){

#source("clusCA.R")
#source("iFCB.r")
#source("MCAk.r")
  
  if(method=="clusCA"){
    out=clusCA(data=data,nclus=nclus,ndim=ndim,nstart=nstart,smartStart=smartStart, gamma = gamma,seed=seed)
  }
  if(method=="iFCB"){
    out=iFCB(data=data,nclus=nclus,ndim=ndim,nstart=nstart,smartStart=smartStart, gamma = gamma,seed=seed)
  }
  if(method=="MCAk"){
    out=MCAk(data=data,nclus=nclus,ndim=ndim,nstart=nstart,alpha = alpha,smartStart=smartStart, seed=seed)
  }
  return(out)
  }



