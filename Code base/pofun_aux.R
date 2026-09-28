
rupo<-function(n,N,S,T,DRAW=TRUE) {
  
  #start state
  m<-matrix(0,n,n)
  m[upper.tri(m)]<-1
  rownames(m)<-colnames(m)<-1:n
  mc <- my.transitive.closure(m)
  mr <- transitive.reduction(mc)
  wold<-sum(mc-mr)

  #showDAG(mr,edge.color='black',vertex.label.cex=0.6,
  #        vertex.color=NA,vertex.size=0,edge.arrow.size=0)
  #title('MCMC start state')
  #readline("Here is the start state. Type return to continue.")
  
  POu<-vector('list',N)
  
  pb<-txtProgressBar(min=0, max=T/S, style=3)
  
  for (t in 1:T) {
    v<-sample(1:n,2)
    mp<-m
    mp[v[1],v[2]]<-1-mp[v[1],v[2]]
    
    if (mc[v[1],v[2]]!=mr[v[1],v[2]]) {
      m<-mp #if we toggled an edge that doesnt change the closure, accept
    } else {
      mpc <- my.transitive.closure(mp)
      if (all(mpc+t(mpc)<2)) { #check no loops
        mpr <- transitive.reduction(mpc)
        wnew<-sum(mpc-mpr)
        if (runif(1)<2^(wold-wnew)) { #accept/reject step
          wold<-wnew
          mc<-mpc
          mr<-mpr 
          m<-mp #accept the proposal and save closure and reduction
        } 
      }
    }
    
    if (t%%S==0){
      #save the marginal process on t-reductions
      POu[[t/S]]<-mc 
      if (DRAW) {
        showDAG(mr,edge.color='black',vertex.color=NA,
                vertex.size=20,edge.arrow.size=0.25)
      }
      setTxtProgressBar(pb, t/S)
      flush.console()
    }
  }
  return(POu)
}

intersect.SO<-function(y){
  N=length(y)
  items=unique(unlist(y))
  n=length(items)
  h=matrix(0,n,n)
  row.names(h)<-colnames(h)<-items
  for (i in 1:N) {
    o=match(y[[i]],as.numeric(rownames(h))) 
    ni=length(o)
    for (j in 1:(ni-1)) {
      for (k in (j+1):ni) {
        h[o[j],o[k]]<-1
      }
    }
  }
  return(h>0 & t(h)==0)
}
