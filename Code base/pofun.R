showDAG<-function(m=NULL,...) {
  mr<-transitive.reduction(m)
  g<-graph_from_adjacency_matrix(m,mode ="directed")#as(m,"graphNEL")
  gr<-graph_from_adjacency_matrix(mr,mode ="directed")#as(m,"graphNEL")
  #h<-graph_from_adjacency_matrix(my.transitive.closure(m),mode ="directed") #ed 5-4-22 for mnem
  h<-gr #seems like a bug in layout function puts vertices on top of one another
  plot(g,layout=layout_with_sugiyama(h)$layout,...);   	 
}

showDAGTRTC<-function(m=NULL,mr=NULL,mc=NULL,mt=NULL,t=NULL) {
  
  #draw DAG's (adjacency matrix format)
  
  #example
  #p<-0.5
  #n<-15
  #for (k in 1:10) {
  #	m<-randDGDAG(n,p,DAG=TRUE)
  #	mc <- my.transitive.closure(m) #ed 5-4-22 for mnem
  #	mr <- transitive.reduction(mc) 
  #	showDAGTRTC(m=m,mr=mr,mc=mc,t=k)
  #}
  
  pics<-0	
  if (!is.null(m)) {
    pics<-pics+1
    g<-graph_from_adjacency_matrix(m,mode ="directed")#as(m,"graphNEL")
  }
  if (!is.null(mc)) {
    pics<-pics+1
    gc<-graph_from_adjacency_matrix(mc,mode ="directed")#as(mc,"graphNEL")
  }
  if (!is.null(mr)) {
    pics<-pics+1
    gr<-graph_from_adjacency_matrix(mr,mode ="directed")#as(mr,"graphNEL")
  }
  if (!is.null(mt)) {
    pics<-pics+1
    gt<-graph_from_adjacency_matrix(mt,mode ="directed")#as(mt,"graphNEL")
  }
  if (pics==0) {
    warning('3 Null graphs passed to showDAGTRTC')
    return()
  } else {
    par(mfrow=c(1,pics)); 
    if (!is.null(m)) {
      plot(g,main=c('random DAG #',as.character(t))); 
    }
    if (!is.null(mr)) {
      plot(gr,main='reduction')
    }
    if (!is.null(mc)) {
      plot(gc,main='closure')
    }
    if (!is.null(mt)) {
      plot(gt,main='truth')
    }
    
  }
}

showDAGs<-function(B,T,PO,poi=1:T) {
  #display the DAGs in PO arranged in a lattice and labeled by year
  n.poi=length(poi)
  n.row=min(4,floor(sqrt(n.poi)))
  par(mfrow=c(n.row,ceiling(length(poi)/n.row)),mai=c(0,0,0,0))
  for (t in poi) {
    mr<-transitive.reduction(PO[[t]])         #for drawing we want reduction
    showDAG(mr,edge.arrow.size=0.5/length(poi),vertex.color=NA,vertex.frame.color=NA,
            vertex.label.cex=1,vertex.size=20,edge.color='black',vertex.label.family="sans")
    text(0.8,-1,paste('[',as.character(B+t-1),']',sep=''),col=2,cex=1)
  }
}

randDGDAG<-function(n=10,p=0.5,DAG=TRUE) {
  
  #generate a random DAG (Brightwell dbn)
  #output is adjacency matrix
  
  #p<-0.5
  #n<-15
  #m<-randDGDAG(n,p,DAG=TRUE)
  #mc <- my.transitive.closure(m) #ed 5-4-22 for mnem
  #mr <- transitive.reduction(mc) 
  #showDAGTRTC(m=m,mr=mr,mc=mc)
  
  ft<-matrix(0+(runif(n*n)<p),n,n)
  if (DAG) {ft[lower.tri(ft,diag=TRUE)]<-0} else {diag(ft)<-0}
  ft
}

rZPO<-function(n,K=floor(n/2),rho=NA,b=rep(0,n)) {
  #simulate one random PO via a random Z-matrix
  if (is.na(rho)) {rho=rRprior()}
  Sig=matrix(rho,K,K); diag(Sig)=1; 
  U=mvrnorm(n,rep(0,K),Sig)+b; 
  O=latent2order(U); 
  h=order2partial(O,n); 
  return(h)
} 

suborder<-function(h,s){
  return(h[s,s])
}

nle<-lecount

noconflict<-function(mc,la) {
  #test if there are conflicts 
  stop('old noconflict had a bug - see noconflict2')
  n.order<-length(la)
  if (n.order==0) {return(TRUE)} 
  for (k in 1:n.order) {
    o<-la[[k]]$o
    ol<-la[[k]]$ll
    ind<-match(1:n,as.numeric(rownames(mc)));
    mcs<-mc[ind,ind]
    for (i in 2:ol) {
      if (mcs[o[i],o[i-1]]==1) {return(FALSE)} 
    }
  }
  TRUE
}


seq2dag<-function(o,n,p=1) {
  m<-matrix(0,n,n)
  for (k in 2:length(o)) {
    m[o[k-1],o[k]]<-p
  }
  m
}

dagwidth<-function(tc) {
  #compute dag width from closure; uses library(igraph)
  #tc<-my.transitive.closure(tc) #ed 5-4-22 for mnem
  op=1-(tc+t(tc))
  opNEL<-as(op,'graphNEL') #XXX TODO this may not be right now
  opIG<-igraph.from.graphNEL(opNEL)
  clique.number(opIG)
}


dagdepth<-function(tr) {
  #compute dag depth from transitive reduction
  if (length(tr)==1) {return(1)}
  n<-dim(tr)[1]
  cs<-apply(tr,2,sum)
  if (sum(cs)==0) {return(1)}
  csi<-(cs==0)
  bs<-apply(tr,1,sum)
  bsi<-(bs==0)
  free<-which(bsi&csi)
  k<-length(free)
  if (k==n) {return(1)}
  if (k>0) { 
    tr<-tr[-free,-free]
    cs<-apply(tr,2,sum)
    csi<-(cs==0)
    bs<-apply(tr,1,sum)
    bsi<-(bs==0)
  }
  tops<-which(csi)
  bots<-which(bsi)
  if (length(bots)>length(tops)) {tops<-bots}
  return(1+dagdepth(tr[-tops,-tops]))
}

latent2order<-function(z) {
  #input: z is n by k
  #output: k perms of 1:n equal rank in z columns
  n<-dim(z)[1]
  k<-dim(z)[2]
  resm<-n-apply(z,2,rank)+1
  resl<-list()
  for (i in 1:k) {
    resm[resm[,i],i]<-1:n
    resl[[i]]<-resm[,i]
  }
  resl
}

order2partial<-function(v,n=NULL) {
  #this
  #output is the transitive closure of 
  #the intersection of the list of complete orders v
  if (is.null(n)) {n<-max(v)}
  w<-lapply(v,seq2dag,n)
  x<-lapply(w,my.transitive.closure) #ed 5-4-22 for mnem
  z<-matrix(0,n,n); 
  colnames(z)<-rownames(z)<-1:n
  for (y in x) {z<-z+y}
  0+(z==length(v))
}

unifLE<-function(tc,le=NULL) {
  #given PO=tc sample one LE uniformly at random
  if ((dim(tc)[1]==1) && (dim(tc)[2]==1)) {
    le<-c(le,as.numeric(rownames(tc)[1]))    
    return(le)
  } else {
    tops<-which(apply(tc,2,sum)==0)
    n.tops=length(tops)
    if (n.tops==1) {
      v=tops
    } else {
      weights=rep(NA,n.tops)
      for (k in 1:n.tops) {
        weights[k]=nle(tc[-tops[k],-tops[k]]) #do I need the transivie closure?
      }
      v<-sample(x=tops,size=1,prob=weights)
    }
    le<-c(le,as.numeric(rownames(tc)[v]))
    tc<-tc[-v,-v,drop=FALSE]
    return(unifLE(tc=tc,le=le))
  }
}

PunifLE<-function(tc,le=NULL,p=0) {
  #given PO=tc sample one LE from the 'lkddown' observation model in which individuals can be promoted up the rank
  if ((dim(tc)[1]==1) && (dim(tc)[2]==1)) {
    le<-c(le,as.numeric(rownames(tc)[1]))    
    return(le)
  } else {
    if (runif(1)<p) {
	    v<-sample(x=1:(dim(tc)[1]),size=1)
    } else {
      tops<-which(apply(tc,2,sum)==0)
      n.tops=length(tops)
      if (n.tops==1) {
        v=tops
      } else {
        weights=rep(NA,n.tops)
        for (k in 1:n.tops) {
          weights[k]=nle(tc[-tops[k],-tops[k],drop=FALSE]) #do I need the transitive closure?
        }
        v<-sample(x=tops,size=1,prob=weights)
      }
    }
    le<-c(le,as.numeric(rownames(tc)[v]))
    tc<-tc[-v,-v,drop=FALSE]
    return(PunifLE(tc,le,p))
  }
}


p1<-function(tc,first.or.last='first') {
  #given PO=tc calculate the probability each node is 'first' or 'last'
  if ((dim(tc)[1]==1) && (dim(tc)[2]==1)) {
    weights=1
  } else {
    weights=0*tc[1,] #to pick up column names
    if (first.or.last=='first') {
      loc<-which(apply(tc,2,sum)==0) #top nodes
    } else {
      loc<-which(apply(tc,1,sum)==0) #bottom nodes
    }
    n.loc=length(loc)
    if (n.loc==1) {
      weights[loc]=1
    } else {
      for (k in 1:n.loc) {
        weights[loc[k]]=nle(tc[-loc[k],-loc[k],drop=FALSE])
      }
    }
  }
  return(weights/sum(weights))
}

my.transitive.closure<-function(h) {
  hc=transitive.closure(h)
  diag(hc)<-0
  return(hc)
}

is.bucket<-function(tr) {
  #test for bucket order - input can be tr or tc
  if (length(tr)==1) {return(TRUE)}
  n<-dim(tr)[1]
  cs<-apply(tr,2,sum)
  if (sum(cs)==0) {return(TRUE)} #if empty it is a bucket
  csi<-(cs==0) #true for top indices
  bs<-apply(tr,1,sum)
  bsi<-(bs==0) #true for bottom indices
  free<-which(bsi&csi)
  k<-length(free)
  if (k>0) {return(FALSE)} #if there are any free nodes (and not all free) then it is at best a parallel bucket
  tops<-which(csi)
  ft=tr[tops[1],] #take the relns of one of the tops
  if (all(apply(tr[tops,,drop=FALSE],1,function(x){all(x==ft)}))) { #check the tops all have the same relns to non-tops
    return(is.bucket(tr[-tops,-tops,drop=FALSE]))
  } else {return(FALSE)}
}

is.vsp<-function(tc) {
  tcg<-graph_from_adjacency_matrix(tc,mode ="directed")
  b<-matrix(c(0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0),4,4,byrow=TRUE)
  bg<-graph_from_adjacency_matrix(b,mode ="directed")
  return(!subgraph_isomorphic(bg,tcg,induced=TRUE))
}

##########################################################################
#enumerate LE's of a PO may be useful for Mallows noise etc

le.expand<-function(tr) {
  
  #enumerate linear extensions of the partial order
  #with transitive reduction tr (binary matrix)
  #for i=1,...,nle(h), the entries in les[,i] are one LE
  #they use the row names of h so if we set 
  #rownames(h)<-colnames(h)<-1:dim(h)[1]
  #then we get vanilla LE's but for PO analysis LE's use bishop
  #labels which need not be 1:dim(h)[1]
  
  if (length(tr)==1) {les=matrix(as.numeric(rownames(tr)[1])); return(les)}
  n<-dim(tr)[1]
  cs<-apply(tr,2,sum)
  csi<-(cs==0)
  tops<-which(csi)
  les=c()
  for (i in tops) {
    trr<-tr[-i,-i,drop=FALSE]
    les.i=le.expand(trr)
    n.i=dim(les.i)[2]
    v=rep(as.numeric(rownames(tr)[i]),n.i)
    les=cbind(les,rbind(v,les.i))
  }
  rownames(les)<-NULL
  return(les)
}

enum.seq<-function(x) {
  
  #this just enumerates all permutations of the entries in x and puts them in a n! x n matrix 
  
  n=length(x)
  if (n==1) return(x)
  seq=c()
  for (i in 1:n) {
    seq=rbind(seq,cbind(rep(x[i],factorial(n-1)),enum.seq(x[-i])))
  }
  return(seq)
}

#hamming distance calculation for the list o respect to PO (tr)
hamming <-function(tr,o){
  if (length(tr)==1) {les=c(as.numeric(o!=rownames(tr)[1])); return(les)}
  cs<-apply(tr,2,sum)
  csi<-(cs==0)
  tops<-which(csi)
  les=c() 
  
  for (i in tops) {
    trr<-tr[-i,-i,drop=FALSE]
    oo <- o[-1,drop=FALSE]
    les.i=hamming(trr,oo)
    n.i=nle(trr)
    v=rep(as.numeric(o[1]!=rownames(tr)[i]),n.i)
    les=c(les,(v+les.i))
  }
  return(les)
}
