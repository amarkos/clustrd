iFCB<- function(data,nclus=3,ndim=2,nstart=100,smartStart=NULL,seed=1234){
  #  require("dummies")
  #  source("EmptyKmeans.r")  
  
  minobs = min(sapply(apply(data,2,unique),length))
  maxobs = max(sapply(apply(data,2,unique),length))
  
  q=ncol(data)
  n=nrow(data)
  
  dZ=as.matrix(dummy.data.frame(data,dummy.classes = "ALL"))
  print(dZ[1:5,1:5])
  dZ=as.matrix(dummy.data.frame((data),drop=F))
  print(dZ[1:5,1:5])
  ndZ = ncol(dZ)
  
  best_f=1000000
  fvec=c()
  for(b in 1:nstart){
    
    # Starting method
    if(is.null(smartStart)){
      myseed=seed+b
      set.seed(myseed)
      randVec= matrix(ceiling(runif(n)*nclus),n,1)
    }else{
      randVec=smartStart
    }
    
    
    C=dummy(randVec)
    
    w= -Inf
    ceps=0.00000001
    itmax=100
    it=0 ### inizializzazione contatore iterazioni
    imp=100000 
    f0 =10000000 ## inizializzazione criterio di arresto
    u=matrix(1,n,1);## vettore di 'uno' 
    #story.obscoord=list()
    
    while((it<=itmax)&&(imp>ceps)){
      
      it=it+1 
      #print("it")
      #print(it)
      
      Fmat = t(C) %*% dZ #crossprod(C,dZ) #
      P=Fmat/sum(Fmat)
      r= rowSums(P)
      c= colSums(P) #apply(P,2,sum)
      
      #r=t(t(r))
      #c=t(t(c))
      onec=matrix(1,nrow=ncol(P))
      nsSpc=sqrt(q)*   t(t(t(t(P)*as.vector(1/c)) - r %*% t(onec)) * as.vector(sqrt(c)))
      #
      
      nssvdres=svd(nsSpc)
      
      # print("decomposition done")
      
      nU  = nssvdres$u[,1:ndim]
      nV  = nssvdres$v[,1:ndim]
      nsv = nssvdres$d[1:ndim]
      
      G =  t(t(nU)*nsv)
      
      B =   (1/(sqrt(as.vector(c)*sqrt(q)))) * t(t(nV) * nsv)  
      
      Csize = colSums(C)
      Cw = as.vector(C %*% Csize)
      Y = Cw * (dZ / (n*sqrt(q))) %*% t(t(B) * as.vector(1/nsv))
      
      outK=try(kmeans(Y,centers=G,nstart=100),silent=T)
      if(is.list(outK)==F){
        #print("singleton")
        outK=EmptyKmeans(Y,centers=G)  
      #  break
      }
      
      G=outK$centers
      ngvec = outK$cluster
      C = dummy(ngvec)
      centerC = matrix(apply(C,2,sum),nrow=n,ncol=nclus,byrow=T)/n
      Zstar = dZ - matrix(c*(n*q)/n,nrow=n,ncol = ndZ,byrow=T)/n
      Cstar = C - centerC
      
      fA=Cstar- t(t(Zstar) * as.vector(sqrt(c))) %*% B %*% t(nU) 
      flossA=sum(diag(t(fA) %*% fA))  
      flossB=sum(diag(t(Y - Cstar %*% G) %*% (Y - Cstar %*% G)))    
      f=flossA+flossB 
      imp=f0-f
      fvec=c(fvec,f)
      f0=f
      #lambda=nsv
    }
    
    if(f<=best_f){
      
      #####gamma scaling
      distB = sum(diag(t(B)%*%  B))
      distG = sum(diag(t(G)%*% G))
      gamma = ((nclus/q)* distB/distG)^.25
      
      B = (1/gamma)*B
      G = gamma*G
      Y = gamma*Y
      #########################
      
      
   #   best_lam=lambda
      best_f=f
      best_ngvec=ngvec	
      best_Y=Y
      best_B=B
      best_G=G
      best_it=it
    }
    
    
    
  }## end FOR
  
  #####################
  ####################
  #lambda=best_lam
  f=best_f
  ngvec=best_ngvec	
  Y=best_Y	
  B=best_B
  G=best_G
  
  
  it=best_it
  # inert=sum(lambda[1:2]^2)
  #  t_inert=sum(lambda^2)
  
  ####################
  #####################
  if ((minobs==2) & (maxobs==2)) {
    B=B[seq(from=1,to=(2*q),by=2),]
  }
  
  
  #################################  
  #  C=dummy(ngvec)
  #  Dzi=diag(colSums(dZ)^-1)
  #  MZ=scale(dZ,scale=FALSE)
  #  MzDzi= t(C)%*%MZ%*%Dzi
  #  GB=G%*%t(B)
  #  MzDzi_GB=MzDzi-GB
  #  final_loss=sum(diag(t(MzDzi_GB)%*%MzDzi_GB))
  #################################  
  
  cluID = ngvec
  #library(plyr)
  ##reorder cluster membership according to cluster size
  csize = round((table(cluID)/sum( table(cluID)))*100,digits=2)
  aa = sort(csize,decreasing = TRUE)
  cluID = mapvalues(cluID, from = as.integer(names(aa)), to = as.integer(names(table(cluID))))
  #reorder centroids
  G = G[as.integer(names(aa)),]
  
  out=list() 
  out$obscoord=Y
  out$attcoord=B
  out$centroid=G
  out$cluID=cluID
  #  out$final_loss=final_loss
  out$criterion=f
  #  out$iters=it
  #  out$expl_inertia= (inert/t_inert)
  #  out$lambda=lambda
  # out$critvec=fvec
  out$csize=round((table(cluID)/sum( table(cluID)))*100,digits=1)
  out$odata=data.frame(lapply(data.frame(data),factor))
  class(out)="clusmca"
  return(out)
}
