clusCA <- function(data,nclus,ndim,nstart=100,smartStart=NULL,gamma = FALSE, seed=1234){
  K = nclus
  k = ndim
  nrs = nstart
  q = ncol(data)
  maxiter = 100
  maxinert=-1
  data=data.frame(data)
  Z=dummy.data.frame(data, dummy.classes = "ALL") # The original super indicator
  n=nrow(Z)
  Q=ncol(Z)
  
  Dzh=diag(as.vector(colSums(Z)^(.5))) 
  #Dzhi=pseudoinverse(Dzh)
  Dzhi= chol2inv(chol(Dzh))
  MZ=scale(Z,scale=FALSE)
  MZD = MZ %*% Dzhi
  
  # Do nrs random starts
  fvec=c()
  for (rs in 1:nrs){
    if(is.null(smartStart)){
      myseed=seed+rs
      set.seed(myseed)
      randVec= matrix(ceiling(runif(n)*nclus),n,1)
    }else{
      randVec=smartStart
    }
    
    Zki=dummy(randVec)
    Dk=t(Zki) %*% Zki
    Dks=Dk^(.5)
    #Dksi=pseudoinverse(Dks)
    Dksi= chol2inv(chol(Dks))
    DZkZD= sqrt(n/q)*Dksi%*% t(Zki) %*% MZD     # equation 6 on the paper
    svdDZkZD=svd(DZkZD)
    #this is for ndim = 1 to work
    if (k != 1) {
      Lk=diag(svdDZkZD$d[1:k])
    } else {
      Lk = data.matrix(svdDZkZD$d[1])
    }
    G = svdDZkZD$u
    Gi = G[,1:k]
    
    Gi=Dksi%*%Gi%*%Lk # CA row coordinates (section 2 in the paper right below formula (1))
    
    Bstar=svdDZkZD$v
    B=sqrt(n*q)*Dzhi %*% Bstar # as in eq. (8)
    
    Bns=B[,1:k]          
    Bi=Bstar[,1:k]           # attribute quantifications. Orthonormal
    
    #              
    
    #  inertia = sum(t(Lk)%*% Lk) # explained inertia in k dimensions
    
    Yi=sqrt((n/q)) * MZD%*%Bi  # The coordinates for the subjects. as in eq. (10).
    
    GDGbef=sum(diag(t(Gi) %*% Dk %*% Gi))   # Objective value before K means This is equivalent to formula (6), but then in the paper
    
    objbef=GDGbef
    ######## END first fixed C step: Given random C, B and G are optimal
    #  obj_improv=c()
    #  kmeanobj=1000        # initialize  obj value
    improv=10
    iter=0   
    
    objective=sum(diag((t(Yi)%*%Zki%*%chol2inv(chol(t(Zki)%*%Zki))%*%t(Zki)%*%Yi)))
    while ((improv > 0.0001) && (iter<maxiter)){
      iter = iter + 1
      outK = try(kmeans(Yi,centers=Gi,nstart=100),silent=T)
      #empty clusters
      if(is.list(outK) == F){
        outK=EmptyKmeans(Yi,centers=Gi) 
        #  break 
      }
      Zki=outK$cluster
      Gi=outK$centers
      
      Zki = dummy(Zki)
      Dk = t(Zki) %*% Zki        
      Dks = Dk^(.5)                  # New Dch weights
      #Dksi = pseudoinverse(Dks)
      Dksi = chol2inv(chol(Dks))
      #END Fixed B step: Given B, C is optimal.
      #Now: Fix C and recalculate B
      DkZkZD=sqrt(n/q)*Dksi%*% t(Zki)%*%MZD # Make the new matrix for the SVD.
      
      outDkZkZD=svd(DkZkZD)       # New SVD, CA analysis
      
      G=outDkZkZD$u
      Gi=G[,1:k]
      # this is for ndim = 1 to work
      if (k != 1) {
        Lk=diag(outDkZkZD$d[1:k])
      } else {
        Lk = data.matrix(outDkZkZD$d[1])
      }
      
      Gi = Dksi %*% Gi %*% Lk
      Bstar=outDkZkZD$v
      Bi=Bstar[,1:k]
      B=sqrt(n*q)*Dzhi%*%Bi
      Bns=B[,1:k]
      #Attribute quantifications
      Yi= sqrt((n/q)) * MZD %*% Bi # Subject coordinates
      lambda = outDkZkZD$d
      objective=c(objective,sum(diag(t(Lk)%*% Lk)))
      
      # B rescaled in such a way that B'DzB=nqI.
      improv=sum(diag(t(Gi)%*% t(Zki)%*%Zki%*%Gi))-objbef
      objbef=sum(diag(t(Gi) %*% t(Zki) %*% Zki %*% Gi))
      
      #varsi=sum(diag(t(Gi)%*% Dk %*% Gi)) #no need to calc inside the loop
      
    }
    #FIX: 16/09/2016, calculated outside the loop
    varsi=sum(diag(t(Gi)%*% Dk %*% Gi))
    
    fvec = c(fvec, objbef)
    if (varsi > maxinert){ #gamma
      if (gamma == TRUE) { 
        distB = sum(diag(t(Bns)%*%  Bns))
        distG = sum(diag(t(Gi)%*% Gi))
        g = ((K/Q)* distB/distG)^.25
        
        Bsol = (1/g)*Bns
        Gsol = g*Gi
        Ysol = g*Yi
      } else {
        Bsol = Bns
        Gsol = Gi
        Ysol = Yi
      }
      Csol = Zki
      maxinert = varsi
    }
  }
  
  wone = which(Csol==1,arr.ind=T)
  cluster = matrix(0,n,1)
  cluster[wone[,1]] = wone[,2]
  
  ##reorder cluster membership according to cluster size
  #csize = round((table(cluster)/sum( table(cluster)))*100,digits=2)
  #aa = sort(csize,decreasing = TRUE)
  size = table(cluster)
  aa = sort(size,decreasing = TRUE)
  cluster = mapvalues(cluster, from = as.integer(names(aa)), to = as.integer(names(table(cluster))))
  #reorder centroids
  Gsol = Gsol[as.integer(names(aa)),]
  
  out=list()
  out$obscoord=Ysol # observations coordinates
  out$attcoord=Bsol # attributes coordinates
  out$centroid=Gsol # centroids
  cluster = as.integer(cluster)
  names(cluster) = rownames(data) 
  out$cluster=cluster   # cluster membership
  out$criterion=maxinert # criterion
  #  out$iters=iters # number of iterations
  out$size=as.integer(aa) #round((table(cluster)/sum( table(cluster)))*100,digits=1)
  out$odata=data.frame(lapply(data.frame(data),factor))
  out$nstart = nstart
  class(out)="clusmca"
  return(out)
}
