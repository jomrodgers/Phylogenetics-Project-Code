#modified source code for additional plotting capability

ln.mean<-function(x){
  if(x[1]==x[2]) return(x[1])
  else {
    a<-x[2]
    b<-log(x[1])-log(x[2])
    return(a/b*exp(b)-a/b)
  }
}





plot.multirateBM_heat<-function(x,digits=1,...){
  cols<-setNames(heat.colors(1000),
                 1:1000)
  est.sig2<-apply(x$tree$edge,1,function(e,x) 
    ln.mean(x[e]),x=x$sig2)
  ln.sig2<-log(est.sig2)
  min.sig2<-min(ln.sig2)
  max.sig2<-max(ln.sig2)
  edge.states<-vector()
  for(i in 1:length(est.sig2)){
    edge.states[i]<-round((ln.sig2[i]-min.sig2)/
                            (max.sig2-min.sig2)*999)+1
  }
  tree<-paintBranches(x$tree,edge=x$tree$edge[1,2],
                      state=edge.states[1])
  for(i in 2:length(edge.states))
    tree<-paintBranches(tree,edge=tree$edge[i,2],
                        state=edge.states[i])
  plotTree(x$tree,plot=FALSE,...)
  lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  obj<-plot(tree,colors=cols,lwd=3,split.vertical=TRUE,
            xlim=c(-0.3,1.05)*diff(lastPP$x.lim),
            add=TRUE,...)
  lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  h<-max(nodeHeights(x$tree))
  LWD<-diff(par()$usr[1:2])/dev.size("px")[1]
  lines(x=rep(-0.25*h+LWD*15/2,2),y=c(2,Ntip(x$tree)-1))
  nticks<-10
  Y<-cbind(seq(2,Ntip(x$tree)-1,length.out=nticks),
           seq(2,Ntip(x$tree)-1,length.out=nticks))
  X<-cbind(rep(-0.25*h+LWD*15/2,nticks),
           rep(-0.25*h+LWD*15/2+0.02*h,nticks))
  for(i in 1:nrow(Y)) lines(X[i,],Y[i,])
  ticks<-exp(seq(min(log(x$sig2)),max(log(x$sig2)),
                 length.out=nticks))
  add.color.bar(Ntip(x$tree)-3,
                heat.colors(1000),
                title=expression(paste("evolutionary rate ( ",sigma^2,")")),
                lims=NULL,digits=3,
                direction="upwards",
                subtitle="",lwd=15,
                x=-0.25*h,
                y=2,prompt=FALSE)
  text(x=X[,2],y=Y[,2],signif(ticks,digits),pos=4,cex=0.7)
}



plot.multirateBM_hcl<-function(x,digits=1,...){
  cols<-setNames(rev(hcl.colors(1000)),
                 1:1000)
  est.sig2<-apply(x$tree$edge,1,function(e,x) 
    ln.mean(x[e]),x=x$sig2)
  ln.sig2<-log(est.sig2)
  min.sig2<-min(ln.sig2)
  max.sig2<-max(ln.sig2)
  edge.states<-vector()
  for(i in 1:length(est.sig2)){
    edge.states[i]<-round((ln.sig2[i]-min.sig2)/
                            (max.sig2-min.sig2)*999)+1
  }
  tree<-paintBranches(x$tree,edge=x$tree$edge[1,2],
                      state=edge.states[1])
  for(i in 2:length(edge.states))
    tree<-paintBranches(tree,edge=tree$edge[i,2],
                        state=edge.states[i])
  plotTree(x$tree,plot=FALSE,...)
  lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  obj<-plot(tree,colors=cols,lwd=3,split.vertical=TRUE,
            xlim=c(-0.3,1.05)*diff(lastPP$x.lim),
            add=TRUE,...)
  lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  h<-max(nodeHeights(x$tree))
  LWD<-diff(par()$usr[1:2])/dev.size("px")[1]
  lines(x=rep(-0.25*h+LWD*15/2,2),y=c(2,Ntip(x$tree)-1))
  nticks<-10
  Y<-cbind(seq(2,Ntip(x$tree)-1,length.out=nticks),
           seq(2,Ntip(x$tree)-1,length.out=nticks))
  X<-cbind(rep(-0.25*h+LWD*15/2,nticks),
           rep(-0.25*h+LWD*15/2+0.02*h,nticks))
  for(i in 1:nrow(Y)) lines(X[i,],Y[i,])
  ticks<-exp(seq(min(log(x$sig2)),max(log(x$sig2)),
                 length.out=nticks))
  add.color.bar(Ntip(x$tree)-3,
                rev(hcl.colors(1000)),
                title=expression(paste("evolutionary rate ( ",sigma^2,")")),
                lims=NULL,digits=3,
                direction="upwards",
                subtitle="",lwd=15,
                x=-0.25*h,
                y=2,prompt=FALSE)
  text(x=X[,2],y=Y[,2],signif(ticks,digits),pos=4,cex=0.7)
}

plot.multirateBM_topo<-function(x,digits=1,...){
  cols<-setNames(topo.colors(1000),
                 1:1000)
  est.sig2<-apply(x$tree$edge,1,function(e,x) 
    ln.mean(x[e]),x=x$sig2)
  ln.sig2<-log(est.sig2)
  min.sig2<-min(ln.sig2)
  max.sig2<-max(ln.sig2)
  edge.states<-vector()
  for(i in 1:length(est.sig2)){
    edge.states[i]<-round((ln.sig2[i]-min.sig2)/
                            (max.sig2-min.sig2)*999)+1
  }
  tree<-paintBranches(x$tree,edge=x$tree$edge[1,2],
                      state=edge.states[1])
  for(i in 2:length(edge.states))
    tree<-paintBranches(tree,edge=tree$edge[i,2],
                        state=edge.states[i])
  plotTree(x$tree,plot=FALSE,...)
  lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  obj<-plot(tree,colors=cols,lwd=3,split.vertical=TRUE,
            xlim=c(-0.3,1.05)*diff(lastPP$x.lim),
            add=TRUE,...)
  lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  h<-max(nodeHeights(x$tree))
  LWD<-diff(par()$usr[1:2])/dev.size("px")[1]
  lines(x=rep(-0.25*h+LWD*15/2,2),y=c(2,Ntip(x$tree)-1))
  nticks<-10
  Y<-cbind(seq(2,Ntip(x$tree)-1,length.out=nticks),
           seq(2,Ntip(x$tree)-1,length.out=nticks))
  X<-cbind(rep(-0.25*h+LWD*15/2,nticks),
           rep(-0.25*h+LWD*15/2+0.02*h,nticks))
  for(i in 1:nrow(Y)) lines(X[i,],Y[i,])
  ticks<-exp(seq(min(log(x$sig2)),max(log(x$sig2)),
                 length.out=nticks))
  add.color.bar(Ntip(x$tree)-3,
                topo.colors(1000),
                title=expression(paste("evolutionary rate ( ",sigma^2,")")),
                lims=NULL,digits=3,
                direction="upwards",
                subtitle="",lwd=15,
                x=-0.25*h,
                y=2,prompt=FALSE)
  text(x=X[,2],y=Y[,2],signif(ticks,digits),pos=4,cex=0.7)
}



