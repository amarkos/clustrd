cluspca <- function(data,nclus,ndim,alpha=NULL,method="RKM",center = TRUE, scale = TRUE, rotation="none",nstart=100,smartStart=NULL,seed=1234)
{
  #NOTE: FactorialKM needs smartstart k-means or else to perform well
  #FIX: K=2, d=2 does not work for RKM
  method <- match.arg(method, c("RKM", "rkm","rKM","FKM", "fkm","fKM"), several.ok = T)[1]
  method <- toupper(method)
  
  #  If alpha = .5 gives RKM, alpha=1 PCA and alpha =0  FKM. 
  if (is.null(alpha) == TRUE)
  {  
    if (method == "RKM") {
      alpha = .5
    } else if (method == "FKM") {
      alpha = 0
    }
  }
  
  odata = data
  data =  scale(data, center = center, scale = scale)
  
  data = data.matrix(data)
  n = dim(data)[1]
  m = dim(data)[2]
  conv=1e-6  # convergence criterion
  func={}; AA = {}; FF = {}; YY = {}; UU={}
  for (run in c(1:nstart)) {
    
    # Starting method
    if(is.null(smartStart)){
      myseed=seed+run
      set.seed(myseed)
      randVec= matrix(ceiling(runif(n)*nclus),n,1)
    }else{
      randVec=smartStart
    }
    
    U = dummy(randVec)
    #update A
    P = U%*%pseudoinverse(t(U)%*%U)%*%t(U)
    A = eigen(t(data)%*%((1-alpha)*P-(1-2*alpha)*diag(n))%*%data)$vectors
    A = A[,1:ndim]
    
    #update Y
    G = data%*%A
    Y = pseudoinverse(t(U)%*%U)%*%t(U)%*%G
    f = alpha*ssq(data - G%*%t(A))+(1-alpha)*ssq(data%*%A-U%*%Y)
    f = as.numeric(f) #fixes convergence issue 01 Nov 2016
    fold = f + 2 * conv*f
    iter = 0
    #iterative part
    while (f<fold-conv*f) {
      fold=f
      iter=iter+1
      outK = try(kmeans(G,centers=Y,nstart=100),silent=T)
      
      if(is.list(outK)==FALSE){
        outK = EmptyKmeans(G,centers=Y)  
        #  break
      }
      
      v = as.factor(outK$cluster)
      U = diag(nlevels(v))[v,] #dummy cluster membership
      
      # update A
      P = U%*%pseudoinverse(t(U)%*%U)%*%t(U)
      A = eigen(t(data)%*%((1-alpha)*P-(1-2*alpha)*diag(n))%*%data)$vectors
      A = A[,c(1:ndim)]
      G = data %*% A
      #update Y
      Y = pseudoinverse(t(U)%*%U)%*%t(U)%*%G
      # criterion
      f = alpha*ssq(data - G%*%t(A))+(1-alpha)*ssq(data%*%A-U%*%Y)
    }
    func[run] = f
    #fpXunc[run]=fpX
    FF[[run]] = G
    AA[[run]] = A
    YY[[run]] = Y
    UU[[run]] = U
  }
  
  
  ##reorder according to cluster size
  mi = which.min(func)
  U=UU[[mi]]
  cluID = apply(U,1,which.max)
  
  #csize = round((table(cluID)/sum( table(cluID)))*100,digits=2)
  size = table(cluID)
  aa = sort(size,decreasing = TRUE)
  cluID = mapvalues(cluID, from = as.integer(names(aa)), to = as.integer(names(table(cluID))))
  #reorder centroids
  centroid = YY[[mi]]
  centroid = centroid[as.integer(names(aa)),]
  #######################
  
  ### rotation options ###
  if (rotation == "varimax") { #with Kaiser Normalization
    AA[[mi]] = varimax(AA[[mi]])$loadings[1:m,1:ndim]
    FF[[mi]] = data%*%AA[[mi]]
    #update center
    centroid = pseudoinverse(t(U)%*%U)%*%t(U)%*%FF[[mi]] 
    centroid = centroid[as.integer(names(aa)),]
  } else if (rotation == "promax") {
    AA[[mi]] = promax(AA[[mi]])$loadings[1:m,1:ndim]
    FF[[mi]] = data%*%AA[[mi]]
    #update center
    centroid = pseudoinverse(t(U)%*%U)%*%t(U)%*%FF[[mi]] 
    centroid = centroid[as.integer(names(aa)),]
  }
  
  #  distB = sum(diag(t(AA[[mi]])%*% AA[[mi]]))
  #  distG = sum(diag(t(centroid)%*% centroid))
  #  gamma = ((nclus/m)* distB/distG)^.25
  
  #  AA[[mi]] = (1/gamma)*AA[[mi]]
  #  centroid = gamma*centroid
  #  FF[[mi]] = gamma*FF[[mi]]
  
  ##########################
  
  #assign output
  out=list()
  mi = which.min(func)
  out$obscoord = FF[[mi]]
  rownames(out$obscoord) = rownames(data)
  out$attcoord = data.matrix(AA[[mi]])#[1:m,1:ndim]
  rownames(out$attcoord) = colnames(data)
  out$centroid = centroid #YY[[mi]]
  names(cluID) = rownames(data) 
  out$cluID = cluID #apply(U,1,which.max)
  out$criterion = func[mi]
  out$size = as.integer(aa) #round((table(cluID)/sum(table(cluID)))*100,digits=1)
  out$odata = data.frame(odata)
  out$scale = scale
  out$center = center
  out$nstart = nstart
  class(out) = "cluspca"
  return(out)   
}

ssq = function(a) {
  t(as.vector(c(as.matrix(a))))%*%as.vector(c(as.matrix(a)))
}