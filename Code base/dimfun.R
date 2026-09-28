crown.PO<-function(n){
  m=floor(n/2)
  s1=combn(1:m, m-1, simplify = FALSE) 
  s2=combn(1:m, 1, simplify = FALSE) 
  hc=matrix(0,n,n); rownames(hc)<-colnames(hc)<-1:n
  for (i in 1:m) {
    for (j in 1:m) {
      hc[i,m+j]=0+all(s2[[j]] %in% s1[[i]])
    }
  }
  return(hc)
}

check.dec<-function(hc,hc.decompose){
  hc2<-intersect.TO(hc.decompose)
  return(all(hc==hc2))
}

intersect.TO<-function(tos) {
  n<-dim(tos)[1]
  v<-c(as.data.frame(tos)) 
  hc<-order2partial(v,n)
  return(hc)
}

decompose<-function(hc) {
  n=dim(hc)[1]
  if (sum(hc)==n*(n-1)/2) return(le.expand(hc)) #it is a total order
  #if (is.vsp(hc)) return(matrix(NA,n,2))
  noi=apply(hc,2,sum)==0
  noo=apply(hc,1,sum)==0
  float=which(noi&noo); has.float<-length(float)>0
  if (has.float) {
    oldnames=as.numeric(colnames(hc)[-float])
    hc=hc[-float,-float]; 
    n=dim(hc)[1]; 
    rownames(hc)<-colnames(hc)<-1:n
  }
  le<-le.expand(hc); m<-dim(le)[2]
  for (d in 2:(floor(n/2))) {             #increase dim d from 2 to n/2
    s=combn(1:m, d, simplify = FALSE)   #take all subsets of size d
    ns=length(s)
    for (i in 1:ns) {
      j=s[[i]]
      v<-c(as.data.frame(le[,j]))         #convert selected lists to list() format
      hc2<-order2partial(v,n)             #convert lists to PO
      if (all(hc==hc2)) {                 #if PO matches hc then we have a decomposition into lists
        if (has.float) {
          old.le<-apply(le[,j],2,function(x) oldnames[x])
          c1=rbind(as.matrix(rev(float)),old.le[,1,drop=FALSE])
          c2=rbind(old.le[,2:d,drop=FALSE],float%*%matrix(rep(1,d-1),1,d-1))
          deco<-unname(cbind(c1,c2))
          return(deco)
        } else {
          return(le[,j])
        }
      }
    }
  }
  stop('no decomposition found')
  return(NA)
}

dimension<-function(hc) {
  n=dim(hc)[1]
  if (sum(hc)==n*(n-1)/2) return(1) #it is a total order
  if (is.vsp(hc)) return(2)
  noi=apply(hc,2,sum)==0
  noo=apply(hc,1,sum)==0
  float=which(noi&noo); has.float<-length(float)>0
  if (has.float) {hc=hc[-float,-float]; n=dim(hc)[1]; rownames(hc)<-colnames(hc)<-1:n}
  le<-le.expand(hc); m<-dim(le)[2]
  for (d in 2:(floor(n/2))) {             #increase dim d from 2 to n/2
    s=combn(1:m, d, simplify = FALSE)   #take all subsets of size d
    ns=length(s)
    for (i in 1:ns) {
      j=s[[i]]
      v<-c(as.data.frame(le[,j]))         #convert selected lists to list() format
      hc2<-order2partial(v,n)             #convert lists to PO
      if (all(hc==hc2)) {                 #if PO matches hc then we have a decomposition into lists
        return(d)
      }
    }
  }
  stop('no decomposition found')
  return(NA)
}


