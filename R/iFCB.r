iFCB<- function(data,nclus,ndim,nstart=100,smartStart=F){
  data = data.matrix(data)
  q=ncol(data)
  n=nrow(data)
  dZ=disjMake(data)
  Q=ncol(dZ)  
  oldf=1000000
  
  for(b in 1:nstart){
    
    gvec =matrix(ceiling(runif(n)*nclus),n,1)
    if(smartStart==T){
      gvec=kmeans(dZ,nclus,nstart=100)$cluster
    }
    C=matrix(0,n,nclus)
    for (i in 1:n){
      C[i,gvec[i]]=1
    }
    
    w= -Inf
    ceps=0.00000001
    itmax=100
    it=0 
    imp=100000 
    f0 =10000000 
    u=matrix(1,n,1)
    
    while((it<=itmax)&&(imp>ceps)){
      
      it=it+1 
      
      Fmat=t(C) %*% dZ
      P=Fmat/sum(Fmat)
      r= apply(P,1,sum)
      c= apply(P,2,sum)
      
      r=t(t(r))
      c=t(t(c))
      rvec=as.vector(r)
      cvec=as.vector(c)
      Dr=diag(rvec)
      Dc=diag(cvec)
      
      Dz=diag(diag(t(dZ)%*% dZ))
      DQ=Dz/(n*q)
      sqDQ=sqrt(DQ)
      invDQ=solve(DQ)
      invsqDQ=sqrt(invDQ)
      invDr=diag((1/rvec))
      invDc=diag(1/cvec)
      sqDr=diag(sqrt(rvec))
      sqDc=diag(sqrt(cvec))
      invsqDr=diag(sqrt(1/rvec)) 
      invsqDc=diag(sqrt(1/cvec)) 
      oner=matrix(1,nrow=nrow(P))
      onec=matrix(1,nrow=ncol(P))
      onen=matrix(1,n,1)
      eyen=diag(as.vector(onen))
      nsSpc= sqrt(q)*(P %*% invDc - r %*% t(onec)) %*% sqDc
      nssvdres=svd(nsSpc)
      nU=nssvdres$u[,1:ndim]
      nV=nssvdres$v[,1:ndim]
      nDsv=diag(nssvdres$d[1:ndim])
      invnDsv=diag(1/nssvdres$d[1:ndim])
      invnSqDsv=diag(1/sqrt(nssvdres$d[1:ndim]))
      G =  nU %*% nDsv
      B = invsqDc %*% nV %*% nDsv
      nsSCc = invsqDc %*% nV
      
      Csize=diag(t(C)%*% C)
      Cw = as.vector(C %*% Csize)
      Dw = diag(as.vector(Cw))
      Y = Dw %*% (dZ / (n*sqrt(q)))  %*% B %*% invnDsv
      
      outK = kmeans(Y,centers=G)
      ngvec = outK$cluster
      G=outK$centers
      
      C=matrix(0,n,nclus)
      for(i in 1:n){
        C[i,ngvec[i]]=1
      }
      
      centerOp=eyen-onen%*%t(onen)/n
      Cstar=centerOp %*% C
      Zstar=centerOp %*% dZ
      Dzstar=diag(diag(t(Zstar)%*% Zstar))
      
      flossA=sum(diag(t(Cstar-Zstar %*%sqDc %*% B %*% t(nU) ) %*% (Cstar-Zstar %*% sqDc %*% B %*% t(nU))))
      flossB=sum(diag(t(Y - Cstar %*% G) %*% (Y - Cstar %*% G)))
      
      f=flossA+flossB 
      imp=f0-f
      
      
      f0=f
      
    }
    
    if(f<=oldf){
      
      oldf=f
      oldngvec=ngvec	
      oldnewdZ=Y	
      oldattcoord=B  
      oldcentroids=G
      oldG=G
      oldC=C
      oldit=it
    }
    old_Y=Y
    old_B=B
    old_G=G
    
    old_it=it
  }## end FOR
  
  ####################
  f=oldf
  ngvec=oldngvec	
  Y=old_Y	
  B=old_B
  G=old_G
  it=old_it
  #####################
  otpt=list() 
  otpt$obscoord=Y
  otpt$attcoord=B
  otpt$centroid=G
  otpt$cluID=ngvec
  otpt$criterion=f
  otpt
}