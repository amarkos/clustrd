outOfIndependence=function(data,Gvec,labs,nolabs=FALSE,fixmarg=TRUE,firstfew=0,minx=-2.5,maxx=2.5,segSize=4,textSize=6){
  
  value = NULL
  newplace = NULL
  lbls = NULL
  
  data=data.frame(data)
  data=dummy.data.frame(data, dummy.classes = "ALL")
  
  K=max(Gvec)
  C=matrix(0,nrow(data),max(Gvec))
  
  for(j in 1:max(Gvec)){
    C[which(Gvec==j),j]=1
  }
  
  P=t(data) %*% C
  n=nrow(data)
  
  P=P/sum(P)
  
  c=apply(P,2,sum)
  c=t(t(c))
  r=apply(P,1,sum)
  r=t(t(r))
  
  invsqDc=diag(as.vector(1/sqrt(c)))
  invsqDr=diag(as.vector(1/sqrt(r)))
  eP= r%*% t(c)
  devP=invsqDr %*% (P-eP) %*% invsqDc
  
  
  ###### HERE STARTS THE FOR LOOP
  dfP=list()
  sortOp=list()
  bp=list()
  
  colorPal=rainbow(K)
  
  for(jj in 1:K){
    #topfew=which(abs(devP[,jj]*sqrt(n))>1)
    #print(labs[topfew])
    dfP[[jj]]=data.frame(value=devP[,jj]*sqrt(n),place=1:nrow(devP),lbls=labs)
    sortOp[[jj]]=sort(abs(dfP[[jj]]$value),decreasing=T,index.return=T)
    #   sortOp2=sort(abs(dfP2$value),decreasing=T,index.return=T)
    #   sortOp3=sort(abs(dfP3$value),decreasing=T,index.return=T)  
    
    dfP[[jj]]=dfP[[jj]][sortOp[[jj]]$ix,]
    dfP[[jj]]$newplace=nrow(devP):1
    xran=c(min(dfP[[jj]]$value)-.5,max(dfP[[jj]]$value)+.5)
    if(firstfew>0){
      dfP[[jj]]=dfP[[jj]][1:firstfew,]#names(dfP[[jj]])
      dfP[[jj]]$newplace=firstfew:1
    }
    
    bbp=ggplot(data=dfP[[jj]], aes(x=value,y=newplace),labels=lbls)
    
    if(fixmarg==TRUE){
      bbp=bbp+geom_segment(data=dfP[[jj]],aes(x=0,xend=value,y=newplace,yend=newplace),colour=colorPal[jj],size=segSize,alpha=.25)
      bbp=bbp+theme(legend.position="none")+xlab("")+ylab("")+xlim(c(minx,maxx))
      bbp=bbp+theme(axis.text.x  = element_text(size=textSize),axis.text.y  = element_text(size=textSize))
     # bbp=bbp+xlab(paste("Standardized residuals")) + ylab(paste("Variable categories"))  
      if(firstfew==0){bbp=bbp+theme(axis.line=element_blank(),axis.ticks = element_blank())}
      
    }
    else{
      
      bbp=bbp+geom_segment(data=dfP[[jj]],aes(x=0,xend=value,y=newplace,yend=newplace),colour=colorPal[jj],size=segSize,alpha=.25)
      bbp=bbp+theme(legend.position="none")+xlab("")+ylab("")+xlim(xran)
      bbp=bbp+theme(axis.text.x  = element_text(size=textSize),axis.text.y  = element_text(size=textSize))
   #   bbp=bbp+xlab(paste("Standardized residuals")) + ylab(paste("Variable categories"))          
      if(firstfew==0){bbp=bbp+theme(axis.line=element_blank(),axis.ticks = element_blank())}
    }
    if(nolabs==F){
      bbp=bbp+geom_text(data=dfP[[jj]],aes(label=lbls),size=textSize)
    }
    bp[[jj]]=bbp
  }
  #   
  #  
  out=list()
  out$G=bp
  
  out
}
