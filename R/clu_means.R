clu_means<-function(x, id, disp=TRUE, center=TRUE, scale=TRUE){
 # require("dplyr")
#  require("ggplot2")
#  require("cowplot")
#  require("GGally")
#  require("RColorBrewer")
  clu = NULL
  funs = NULL
  
  x =  data.frame(scale(as.matrix(x), center = center, scale = scale))
  
  p=ncol(x)
  gm=apply(x,2,mean)
  
  id=factor(id)
  csize=as.vector(table(id)/sum(table(id)))
  
  x$clu=id
  
  clum=(x %>% group_by(clu) %>% summarise_each(funs(mean)))
  
  am=rbind(clum[,-1],gm)
  ymax=max(apply(am,2,max))
  ymin=min(min(apply(am,2,min)),0)
  bm=data.frame(t(am))
  names(bm)=c(paste("C",1:nrow(clum),sep=""),"all")
  bm$names=row.names(bm)
  
  az=all(bm$all<0.00001)
  
  myplots=list()
  par_bm=data.frame(t(bm[-ncol(bm)]))
  gnam=paste(names(bm)[-ncol(bm)]," (",round(csize*100,digits=1),"%",")",sep="")
#  cnm=paste(cnames,": ",round(csize*100,2),"%",sep="")
  
  gnam[length(gnam)] = "all"
  par_bm$clusters=gnam
  par_bm$csize=c(csize,1/length(csize))
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  mypal=gg_color_hue(length(csize))
  mypal=c("black",mypal)
  
 # if (scale == T) {
#    pco=ggparcoord(par_bm[1:(dim(par_bm)[1]-1),],columns=1:p,groupColumn=p+1,scale="globalminmax",mapping = ggplot2::aes(size = 5*csize)) 
#  } else {
    pco=ggparcoord(par_bm,,columns=1:p,groupColumn=p+1,scale="globalminmax",mapping = ggplot2::aes(size = 3*csize))
#  }
  pco=pco+scale_size_identity()
  pco=pco+scale_colour_manual(values=mypal)
  
#  if (scale == T) {
#    pco=pco+geom_vline(xintercept=1:p,alpha=.5) + xlab("variables") + ylab("z-score") 
#  } else {
    pco=pco+geom_vline(xintercept=1:p,alpha=.1) + xlab("") + ylab("mean") 
#  }
  
  # for(jj in levels(id)){
  #   
  #   jj=as.numeric(jj)
  #   ys=names(bm)[jj]
  #   a=ggplot(data=bm, aes_string(x="names",y=ys,fill="names"))+geom_bar(stat="identity",width=.25)
  #   a=a+xlab("")+ylab("")+guides(fill=FALSE)+ylim(ymin,ymax)
  #   if(az==F){
  #     a=a+geom_linerange(aes(x=names,ymin=0,ymax=all),data=bm)
  #   }else{
  #     a=a+geom_hline(yintercept=0)
  #   }
  #   a=a+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=8))
  #   a=a+theme(axis.text.y = element_text(size=8))+coord_flip()
  #   myplots[[jj]]<-a
  #   
  #   # print(myplots[[jj]])
  # }
  # 
  # cnames=paste("C",1:nrow(clum),sep="")
  # cnm=paste(cnames,": ",round(csize*100,2),"%",sep="")
  # 
  # 
  # b=suppressWarnings(plot_grid(plotlist=myplots,ncol=2,labels=cnm,
  #                              label_size=8))
  if(disp==T){
    #    print(b)
    print(pco)
  }else{
    ggsave(filename = "pcoord_cludesc.pdf",pco,height=5, width=7)
  }
  out=list()
  out$pcoordplot=pco
  return(out)
  
}
