FactorialKM <- function(data,nclus,ndim,nstart=100,smartStart=F){
  data = data.matrix(data) 
  # data standardization
  n = dim(data)[1]
  data = data.frame(scale(data, center = TRUE, scale = TRUE))
  
  # step 0 - initialization
  # starting values for A
  A = princomp(data,cor=TRUE)$loadings[,c(1:ndim)]
  
  oldf = 1000000                   
  
  for (b in 1:nstart){
    ceps=0.00000001
    itmax=100
    it=0 
    
    imp=100000 
    f0 =10000000 
    while ((it <= itmax) && (imp > ceps)){
      it = it+1
      #step 1 - clustering
      P = data.matrix(data)%*%A
      
      if(smartStart==T){
        outkstart=kmeans(P,nclus,nstart=100)
      } else {
        outkstart=kmeans(P,nclus)        
      }
      
      cluID=as.factor(outkstart$cluster)
      U = diag(nlevels(cluID))[cluID,] #dummy cluster membership
      UU = solve(t(U)%*%U) #projector
      Y = UU%*%t(U)%*%P #factor centroid
      
      # step 2 - update loadings
      m = (t(data.matrix(data))%*%(U%*%UU%*%t(U)-diag(n))%*%data.matrix(data))
      A = svd(m)$v[,c(1:ndim)]
      
      #checking convergence
      Fi = norm(P-U%*%Y)
      imp = f0 - Fi
      f0 = Fi
    }
    
    if (Fi <= oldf){
      oldf = Fi
      Uold = U
      Aold = A
      Yold = Y
      Fold = P
      indexold = data.frame(cluID)
    }
    
  }
  f = oldf
  U = Uold
  A = Aold
  Y = Yold
  P = Fold
  index = indexold
  
  #assign output
  out=list()
  out$obscoord=P
  out$attcoord=A
  out$centroid=Y
  out$cluID=as.numeric(data.matrix(cluID))
  out$criterion=f
  out  
}