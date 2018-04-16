MCAk <- function(data, nclus = 3, ndim = 2, alphak = .5, nstart = 100, smartStart=NULL,gamma = TRUE, seed=1234){
  out=list()
  data=data.frame(data)
  if (alphak == 1) { #Tandem approach MCA + k-means
   
    n = nrow(data)
    #asymmetric map, biplot
    A = mjca(data)$colcoord[,1:ndim]
    Fm = mjca(data)$rowpcoord[,1:ndim]
    
    if(is.null(smartStart)){
      set.seed(seed)
      randVec= matrix(ceiling(runif(n)*nclus),n,1)
    }else{
      randVec=smartStart
    }
    
    U = dummy(randVec)
    center = pseudoinverse(t(U) %*% U) %*% t(U) %*% Fm 
    outK = try(kmeans(Fm, centers = center, nstart = 100), silent = TRUE)
    
    center=outK$centers
    
    if(is.list(outK) == FALSE){
      outK = EmptyKmeans(Fm,centers=center)  
     # break
    }
    
    index = outK$cluster
    cluster = as.numeric(index)
    cluster = as.integer(cluster)
    names(cluster) = rownames(data) 
    size = table(cluster)
    aa = sort(size,decreasing = TRUE)
    
    if (gamma == TRUE) {
      distB = sum(diag(t(A)%*%  A))
      distG = sum(diag(t(center)%*% center))
      g = ((nclus/ncol(data))* distB/distG)^.25
      
      A = (1/g)*A
      center = g*center #is this needed
      Fm = g*Fm
    }
    
    cluster = mapvalues(cluster, from = as.integer(names(aa)), to = as.integer(names(table(cluster))))
    #reorder centroids
    center = center[as.integer(names(aa)),]
    out$obscoord=Fm # observations coordinates 
    out$attcoord=A # attributes coordinates 
    out$centroid = center
    out$cluster = cluster
    out$criterion = 1 # criterion
    out$size=as.integer(aa)  #round((table(cluster)/sum( table(cluster)))*100,digits=1)
    out$odata=data.frame(lapply(data.frame(data),factor))
    out$nstart = nstart
    class(out)="clusmca"
    return(out)
  } else {
    
    zz = as.matrix(dummy.data.frame(data,dummy.classes = "ALL"))
    # data = data.matrix(data)
    data=data.frame(data)
    n = nrow(data)
    zitem = ncol(data)            
    zncati = sapply(data.frame(data), function(x) length(unique(x))) #apply(data,2,max)
    oner = matrix(1,n,1)
    muz  = colMeans(zz)
    z = zz - oner %*% muz  
    rr = colSums(t(zz) %*% zz)
    invsqDr = diag(1/sqrt(rr))
    M = z
    Pz = z %*% invsqDr
    PPz = t(Pz) %*% Pz
    
    svdPz = svd(PPz)
    invsqD = diag((1/sqrt(svdPz$d))) 
    Fm = Pz %*%  svdPz$v %*% invsqD
    F0 = Fm[, 1:ndim]
    
    oldf = 1e+06
    
    for(b in 1:nstart){
      Fv={}
      Fm = F0
      # Starting method
      if(is.null(smartStart)){
        myseed = seed+b
        set.seed(myseed)
        randVec = matrix(ceiling(runif(n)*nclus),n,1)
      }else{
        randVec = smartStart
      }
      
      U = dummy(randVec)
      
      center = pseudoinverse(t(U) %*% U) %*% t(U) %*% Fm 
      itmax = 100
      it = 0
      ceps = 1e-04
      imp = 1e+05
      f0 = 1e+05
      
      while((it <= itmax ) && ( imp > ceps ) ){
        it=it+1
        
        ####################################################################
        ## STEP 1: update of U #############################################
        ####################################################################
        #use Lloyd's k-means algorithm to get the results of Hwang and Takane (2006)
        #outK=try(kmeans(Fm,centers=center,algorithm="Lloyd",nstart=100),silent=T)
        outK = try(kmeans(Fm,centers=center,nstart=100),silent=T)
        
        if(is.list(outK)==F){
          outK = EmptyKmeans(Fm,centers=center)  
          #  break
        }
        center=outK$centers
        index = outK$cluster
        U = dummy(index)
        U0=U-oner %*% colMeans(U)
        uu = colSums(t(U)%*% U)
        invsqDru=diag(1/sqrt(c(rr,uu)))
        
        ####################################################################
        ## STEP 2: update of Fm and Wj ######################################
        ####################################################################
        MU = cbind(alphak*M,(1-alphak)*U0)
        Pzu = MU %*% invsqDru
        
        Pzu[is.nan(Pzu)] <- 0
        PPzu = t(Pzu) %*% Pzu
        svdPzu = svd(PPzu)
        invsqD = diag(1/sqrt(svdPzu$d))
        Fm = Pzu %*% svdPzu$u %*% invsqD
        Fm = Fm[, 1:ndim]
        ft1 = 0
        k = 1
        kk = 0
        kk = kk+zncati[1]
        Tm = z[,k:kk]
        W = pseudoinverse(t(Tm)%*% Tm)%*% t(Tm) %*% Fm 
        A = W
        ft1 = ft1+sum(diag((t(Fm) %*% Fm)-(t(Fm) %*% Tm %*% W)))
        k = kk+1
        for(j in 2:zitem){
          kk = kk+zncati[j]
          Tm = z[,k:kk]
          W = pseudoinverse(t(Tm)%*% Tm)%*% t(Tm) %*% Fm 
          A = rbind(A,W)
          ft1 = ft1 + sum(diag((t(Fm) %*% Fm)-(t(Fm) %*% Tm %*% W))) ## MCA
          k=kk+1
        }
        ft2 = sum(diag((t(Fm) %*% Fm) - (t(Fm) %*% U %*% center)))
        #check again
        f =  alphak*ft1 + (1-alphak)*ft2
        
        imp=f0-f
        f0=f
        Fv = cbind(Fv,f)
      } # end WHILE
      
      if (f <= oldf){
        
        #####gamma scaling
        if (gamma == TRUE) {
          distB = sum(diag(t(A)%*%  A))
          distG = sum(diag(t(center)%*% center))
          g = ((nclus/zitem)* distB/distG)^.25
          
          A = (1/g)*A
          center = g*center #is this needed
          Fm = g*Fm
        }
        #########################
        
        oldF = Fm
        oldindex = index
        oldf = f				
        Uold = U				
        Aold = A
        centerold = center
      }
    } ##end of FOR
    
   
    Fm = oldF
    index = oldindex
    cluster = as.numeric(index)
    f = oldf
    U = Uold
    A = Aold
    # it=itold
    center = centerold
    
    size = table(cluster)
    aa = sort(size,decreasing = TRUE)
    
    cluster = mapvalues(cluster, from = as.integer(names(aa)), to = as.integer(names(table(cluster))))
    #reorder centroids
    center = center[as.integer(names(aa)),]
    
    out$obscoord = Fm # observations coordinates 
    out$attcoord = A # attributes coordinates 
    out$centroid = center # centroids
    cluster = as.integer(cluster)
    names(cluster) = rownames(data) 
    out$cluster = cluster #as.numeric(index) # cluster membership
    out$criterion = f # criterion
    out$size = as.integer(aa)  #round((table(cluster)/sum( table(cluster)))*100,digits=1)
    out$odata = data.frame(lapply(data.frame(data),factor))
    out$nstart = nstart
    class(out) = "clusmca"
    return(out)
  }  
}
