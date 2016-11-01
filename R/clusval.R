clusval<-function(x,dst="full"){
  if(dst=="full"){
    if(class(x)=="cluspca"){
      data = scale(x$odata, center = x$center, scale = x$scale)
      oDist = daisy(data,metric="euclidean")
    }else{
      oDist=daisy(x$odata,metric="gower")
    }
  }else{
    oDist=daisy(x$obscoord,metric="euclidean")
  }
  
  clu_res=cluster.stats(d=oDist,x$cluID,wgap=F,sepindex=F,sepwithnoise=F)
  
  out=list()
  
  out$ch=clu_res$ch
  out$asw=clu_res$avg.silwidth
  #out$crit=x$criterion
  
  return(out)
}
