rRprior<-function(fac=1/6,tol=1e-4) {
  while (1-(rho<-rbeta(1,1,fac))<tol) {} #for numerical stability
  return(rho)
}

rTprior<-function(a=1,b=1) {
  theta=rbeta(1,a,b)
  return(theta)
}

rUprior<-function(active,NF,T,theta,rho) {
  
  NB=dim(active)[1]
  Sig=matrix(rho,NF,NF); diag(Sig)=1;
  
  #We dont need U for inactive bishops but we simulate it anyway and "NA" it out at the end
  U=array(dim=c(NB,NF,T))
  U[,,1]=mvrnorm(NB,rep(0,NF),Sig)
  
  for (t in 2:T){
    U[,,t]=theta*U[,,(t-1)]+mvrnorm(NB,rep(0,NF),Sig*(1-theta^2)) #"*" for eqm U[b,,t]~mvrnorm(1,rep(0,NF),Sig)
  }

  for (t in 1:T) {for (b in 1:NB) {if (!active[b,t]) U[b,,t]<-NA}}

  return(U)
}

ZtoPO<-function(Z,active,B,T,display.Z=FALSE,PO=NA,years=1:T) {
  #Z=state$Z; display.Z=FALSE;PO=NA;years=1:T

  if (is.na(PO[1])) PO=vector("list",T)
  
  for (t in years) {
    active.t<-which(active[,t])
    if (length(active.t)==1) stop('sorry havnt coded to handle case of just one bishop active in a year 1x1 PO incidence matrix')
    NB.t=length(active.t); 
    vp<-latent2order(Z[active.t,,t])
    mp<-order2partial(vp,NB.t); 
    colnames(mp)<-rownames(mp)<-active.t #mp is closure
    PO[[t]]<-mp
  }
  
  if (display.Z) showDAGs(B,T,PO,years=1:T)

  return(PO)
}

rBprior<-function(beta=NA,i=NA,DB=NA,SigB=1,WARN=TRUE) {
  #this is not correct for ordering if note$constrainbeta=TRUE
  if (WARN) warning('rBprior doesnt sort even if note$constrainbeta=TRUE - this would be a problem if used in MCMC')
  if (is.na(beta[1])) {
    beta=rnorm(DB,0,SigB)
  } else {
    beta[i]=rnorm(length(i),0,SigB)
  }
  return(beta)
}

ComputeZ<-function(years=1:T,bishops=1:NB,U,rank,beta) {
  #just compute Z for the years and bishops you want
  #We need the whole of U,rank and beta as the function stands
  #TODO? for speed just pass the part of U, rank, beta needed and rewrite below? Indexing...

  Z=NA*U
  NB=length(bishops) #we may want this for less than all bishops
  for (t in years) {
    for (b in bishops) {
      Z[b,,t]=U[b,,t]+beta[rank[b,t]]
    }
  }
  return(Z)
}

dBprior.unconstrained<-function(beta,Sigb=1) {
  #log beta prior
  db=sum(dnorm(beta,0,Sigb,log=TRUE))
  return(db)  
}

dBprior.constrained<-function(beta,Sigb=1) {
  #log beta prior
  if (any(diff(beta)>0)) return(-Inf)
  db=sum(dnorm(beta,0,Sigb,log=TRUE))
  return(db)  
}

dUprior<-function(active,U,theta,rho,bishops=1:NB,years=1:T) {

  NF=dim(U)[2]
  T=dim(U)[3]
  NB=dim(U)[1]

  #if years is not 1:T (ie not all) we are conditioning on U in the other years so we need contributions from years which follow YOI
  if (length(years)<T) {
    years=unique(c(years,years+1)) #add in the years after
    bad.year.index=which(years>T) 
    if (length(bad.year.index)>0) years=years[-bad.year.index]
  } 
  #this is delicate but doesnt compute unnecessary year-tranitions eg if years=c(2,5) we need 1->2->3 and 4->5->6 but not 3->4

  Sig=matrix(rho,NF,NF); diag(Sig)=1;
  
  #tot=rep(0,NB) #do all bishops but only calculate for BOI in "bishops" index vector
  tot=matrix(0,NB,T)
  #alot of loops! This is due to different bishops starting at different times
  for (t in years) {     #go through years of interest
    for (b in bishops) { #go through bishops of interest
      if (active[b,t]) { #was this bishop active in this year?
        if (t==1 || !active[b,t-1]) { #this is the first year this bishop appears in U
          #tot[b]=tot[b]+dmvnorm(U[b,,t],rep(0,NF),Sig,log=TRUE) #use the VAR(1) initialisation density
          tot[b,t]=dmvnorm(U[b,,t],rep(0,NF),Sig,log=TRUE)
          #-NF*log(2*pi)/2-log(det(Sig))/2-t(u)%*%solve(Sig)%*%u/2 if u=U[b,,t]
        } else {
          #tot[b]=tot[b]+dmvnorm(U[b,,t]-theta*U[b,,t-1],rep(0,NF),Sig*(1-theta^2),log=TRUE) #use the update density
          tot[b,t]=dmvnorm(U[b,,t]-theta*U[b,,t-1],rep(0,NF),Sig*(1-theta^2),log=TRUE)
        }
      }
    }
  }

  return(list(lpU.tot=sum(tot),lpU.all=tot))
}

dRprior<-function(rho,fac=1/6,tol=1e-4) {
  if (rho>1-tol) {return(-Inf)}
  dR=dbeta(rho,1,fac,log=TRUE)-pbeta(1-tol,1,fac,log=TRUE)
  return(dR)
}

dTprior<-function(theta,a=1,b=1) {
  dth=dbeta(theta,a,b,log=TRUE)
  return(dth)
}

dPrior<-function(active,U,beta,theta,rho,tau,p,q,cla,B,E,PhPar,dBprior) {
  tot=dTprior(theta)+dRprior(rho)+dBprior(beta)+
    dUprior(active,U,theta,rho)$lpU.tot+dYprior(cla,tau,B,E)+dPprior(p,PhPar)+dQprior(q,PhPar)
  return(tot)
}

rPprior<-function(PhPar) { #p-prior mean 1/10 - could go a bit harder - 7-4-22 was 1,
  if (PhPar$model=='lkmallow') {
    a=PhPar$p.m$a; b=PhPar$p.m$b
    rpp=runif(1,min=a,max=b)
  } else {
    a=PhPar$p$a; b=PhPar$p$b 
    rpp=rbeta(1,shape1=a,shape2=b)
  }
  return(rpp)
}

dPprior<-function(p,PhPar) {
  if (PhPar$model=='lkmallow') {
    a=PhPar$p.m$a; b=PhPar$p.m$b
    dpp=dunif(p,min=a,max=b,log=TRUE)
  } else {
    a=PhPar$p$a; b=PhPar$p$b 
    dpp=dbeta(p,shape1=a,shape2=b,log=TRUE)
  }
  return(dpp)
}

rQprior<-function(PhPar) { 
  a=PhPar$q$a; b=PhPar$q$b 
  return(rbeta(1,shape1=a,shape2=b))
}

dQprior<-function(q,PhPar) {
  a=PhPar$q$a; b=PhPar$q$b 
  return(dbeta(q,shape1=a,shape2=b,log=TRUE))
}

rYprior<-function(cla,B,E) {
  NL=length(cla)
  tau=rep(NA,NL)
  for (i in 1:NL) {
    a=max(B,cla[[i]]$tl)
    b=min(E,cla[[i]]$tu)
    tau[i]=floor(runif(1,min=a,max=b+1))-B+1
  }
  return(tau)
}

dYprior<-function(cla,tau,B,E,i=NA) {
  #basicly just check tau is admissable
  if (is.na(i)) {i=1:length(cla)}
  tau=tau+B-1
  tot=0
  for (j in i) {
    a=max(B,cla[[j]]$tl)
    b=min(E,cla[[j]]$tu)
    if (tau[j]<a|tau[j]>b) {tot=-Inf; break}
  }
  return(tot)
}

year2list<-function(tau,T) {
  #y2l[[t]] gives lists assigned to year t by tau
  y2l=vector("list",T)
  for (t in 1:T) y2l[[t]]=which(tau==t)
  return(y2l)
}

loglkd2<-function(r=1,mc,la,model,q=1) {
    
  lkddown.fac<-function(r,mc,o) {
    
    #evaluate the log-probability for the
    #placement of o[1] - the first person in the order o
    #given the PO with reduction mr
    
    #mr,o, an order of length 1, just one way to place it
    if (length(o)==1) {return(0)}
    
    #the sub-DAG for the order
    mla<-mc[o,o]
    #the sub-DAG with the first element removed
    mlb<-mc[o[-1],o[-1],drop=FALSE]
    
    #first person may have been placed at random
    fac <- r/length(o)
    
    #if the first person is in a place that does
    #not violate the proposed PO (given by mc) then
    #they may have been placed using the distribution
    #over linear extensions
    if (sum(mla[,1])==0) {
      fac<-fac+(1-r)*nle(mlb)/nle(mla)
    }
    #return the log-likelihood for this placement
    #plus the log-likelihood for the subsequent placements
    return(log(fac)+lkddown.fac(r,mc,o[-1]))
    
  }
  
  lkdup.fac<-function(r,mc,o) {
    
    #evaluate the log-probability for the
    #placement of o[n] - the last person in the order o
    #given the PO with reduction mr
    
    #mr,o, an order of length 1, just one way to place it
    n<-length(o)
    if (n==1) {return(0)}
    
    #the sub-DAG for the order
    mla<-mc[o,o]
    #the sub-DAG with the last element removed
    mlb<-mc[o[-n],o[-n],drop=FALSE]
    
    #last person may have been placed at random
    fac <- r/length(o)
    
    #if the last person is in a place that does
    #not violate the proposed PO (given by mc) then
    #they may have been placed using the distribution
    #over linear extensions
    if (sum(mla[n,])==0) {
      fac<-fac+(1-r)*nle(mlb)/nle(mla)
    }
    #return the log-likelihood for this placement
    #plus the log-likelihood for the subsequent placements
    return(log(fac)+lkdup.fac(r,mc,o[-n]))
    
  }
  
  lkdnat.fac<-function(r,mc,o) {
    
    #evaluate the log-probability for the
    #placement of o[1] - the first person in the order o
    #given the PO with reduction mr
    
    #mr,o, an order of length 1, just one way to place it
    n<-length(o)
    if (n==1) {return(0)}
    
    #the sub-DAG for the order
    mla<-mc[o,o]
    
    #first person may have been placed at random
    fac <- r/length(o)
    
    #if the last person is in a place that does
    #not violate the proposed PO (given by mc) then
    #they may have been placed at random from the legal placements
    if (sum(mla[,1])==0) {
      fac<-fac+(1-r)/sum(apply(mla,2,sum)==0)
    }
    #return the log-likelihood for this placement
    #plus the log-likelihood for the subsequent placements
    return(log(fac)+lkdnat.fac(r,mc,o[-1]))
    
  }
  
  lkd.mallow<-function(r,mc,o) {
    #evaluate the log-probability for mallows model
    n<-length(o)
    if (n==1) {return(0)}
    
    hamming.dist <- hamming(mc[o,o],o) 

    #disagreed elements in any n-length LE; 0(zero disagree), 1,2,...,n (all disagree)
    #derangements
    n.derange=c(1,0)
    for(i in 3:(n+1)){n.derange=c(n.derange,(i-2)*sum(n.derange[i-c(1,2)]))}
    
    #freq of numbers of disagreed elements
    #number of choosing n.derange out of n elements
    f.derange <- c(1,choose(n,c(1:n)))
    
    #normalizing constant
    n.const <- sum(n.derange*f.derange*exp(-r*c(0:n)))
    
    return(log(mean(exp(-r*hamming.dist)/n.const)))
  }
  
  n.order<-length(la)
  llkda<-matrix(0,1,n.order)
  for (k in 1:n.order) {
    ind=match(la[[k]]$o,as.numeric(rownames(mc))) #10/12/19 - handles names that are not packed in time series
    if (model=='lkddown') {llkda[1,k]<-lkddown.fac(r,mc,ind)}
    if (model=='lkdup') {llkda[1,k]<-lkdup.fac(r,mc,ind)}
    if (model=='lkdnat') {llkda[1,k]<-lkdnat.fac(r,mc,ind)}
    if (model=='bidir') {o=la[[k]]$o; llkda[1,k]<-log(QP.LEProb(mc[ind,ind],o,r,q))}
    if (model=='lkmallow') {llkda[1,k]<-lkd.mallow(r,mc,ind)}
    if (model=='prior') {llkda[1,k]<-0}
  }
  return(llkda)
  
}

legaldag2<-function(ed,NB) {
  #loopfree
  m<-matrix(0,NB,NB)
  for (wl in ed) {m<-m+seq2dag(wl$o,NB)}
  m<-0+(m>0)
  mc <- as(my.transitive.closure(m),"matrix") #ed 5-4-22 for mnem
  bad<-mc+t(mc)
  mc[bad>1]<-0
  return(mc)
}

Uadjust<-function(U,active,T,y2l,cla) {
  #adjust an existing U so it doesnt clash any lists
  #just very crude shifting 
  NB=dim(active)[1]
     
  for (t in 1:T){
    mc=legaldag2(cla[y2l[[t]]],NB)
    g<-graph_from_adjacency_matrix(transitive.reduction(mc),mode ="directed")
    height=layout_with_sugiyama(g)$layout[,2]
    height=-1.5+3*(height-min(height))/(1+max(height)-min(height))
    for (b in 1:NB) {
      if (active[b,t]) U[b,,t]=U[b,,t]-mean(U[b,,t])+height[b]
    }
  }
  return(U)
}

log.lkd.ser<-function(h,cla,y2l,years=1:T,p=0,model='lkddown',q=1,cl=NA,f0=NA) {
  tot=vector('list',length(years))
  for (t in years) {
    if (length(y2l[[t]])>0) {
      tot[[t]]=loglkd2(r=p,h[[t]],cla[y2l[[t]]],model,q=q)
    }
  }
  return(tot)
}

#parallel log.lkd  
log.lkd.par<-function(h,cla,y2l,years=1:T,p=0,model='lkddown',q=1,cl,f0) {
  n.years=length(years)
  tot=vector('list',n.years)  
  X=vector('list',n.years)
  for (i in 1:n.years) {
    X[[i]]$h=h[[years[i]]]
    X[[i]]$la=cla[y2l[[years[i]]]]
    X[[i]]$model=model
    X[[i]]$p=p
    X[[i]]$q=q
  } 
  tot[years]=clusterApplyLB(cl, X, f0)
  return(tot)
}

noconflict2<-function(h,y2l,cla,years=1:T) { #h must be trans.closed
  #test if there are conflicts 
  for (t in years) {
    lat=cla[y2l[[t]]]
    mc=h[[t]]
    n.order<-length(lat)
    if (n.order>0) { 
      for (k in 1:n.order) {
        o<-lat[[k]]$o
        ol<-lat[[k]]$ll
        ind=match(o,as.numeric(rownames(mc)))
        for (i in 1:(ol-1)) {
          for (j in (i+1):ol) {
            if (mc[ind[j],ind[i]]==1) return(list(outcome=FALSE,t=t,bl=y2l[[t]][k],b=o[c(j,i)]))
          } 
        }
      }
    }
  }
  return(list(outcome=TRUE))
}

###########################################################################
#QP functions for queue jumping noise both ways

QPunifLE<-function(tc,p=0,q=0) {
  
  #given PO=tc sample one LE from the mixed error model p is prob error and
  #q is prob for (top down) direction of each entry in thelist build
  
  if ((dim(tc)[1]==1) && (dim(tc)[2]==1)) {
    le<-as.numeric(rownames(tc)[1])   
    return(le)
  } else {
    if (runif(1)<q) { #next entry is top down
      v=queue.next(tc,'top',p)
      act=as.numeric(rownames(tc)[v])
      tc<-tc[-v,-v,drop=FALSE]
      completion<-QPunifLE(tc,p,q)
      le<-c(act,completion)
    } else {
      v=queue.next(tc,'bottom',p)
      act=as.numeric(rownames(tc)[v])
      tc<-tc[-v,-v,drop=FALSE]
      completion<-QPunifLE(tc,p,q)
      le<-c(completion,act)
    }
  }
  return(le)
}

queue.next<-function(tc,direct,p) {
  #choose next actor to add to LE from direct(tion) - as in QPunifLE()
  
  if (runif(1)<p) {      #noise doesnt care about direction
    v<-sample(x=1:(dim(tc)[1]),size=1)
  } else {
    if (direct=='top') { #select next entry using top down (lkddown) queue jumping
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
    } else {             #select next entry using top down (lkddown) queue jumping
      bots<-which(apply(tc,1,sum)==0)
      n.bots=length(bots)
      if (n.bots==1) {
        v=bots
      } else {
        weights=rep(NA,n.bots)
        for (k in 1:n.bots) {
          weights[k]=nle(tc[-bots[k],-bots[k],drop=FALSE]) #do I need the transitive closure?
        }
        v<-sample(x=bots,size=1,prob=weights)
      }
    }
  }
  return(v)    
}

QP.LEProb<-function(tc,le,p,q) {
  
  #if tc is a transitively closed PO and le is one list then calculate the likelihood for the list
  #given the PO tc, and p (err prob) and q (prob choose top down at a given le entry insertion)
  #call this on the suborder so dim(tc)[1]==length(le)
  
  #returns a probability not a log probability as it is a sum
  
  if (!(dim(tc)[1]==length(le))) 
    stop('err in QP.LEProb dim(tc)[1]!=length(le)')
  
  if ( length(setdiff(le,as.numeric(rownames(tc))))>0 | length(setdiff(as.numeric(rownames(tc)),le))>0 ) 
    stop('err in QP.LEProb le and tc rownames dont match')
  
  n=length(le)
  if (n==1) return(1)
  
  ntc=NA
  if (q>0) {
    top.i=which(as.numeric(rownames(tc))==le[1])
    if (length(top.i)!=1) stop('err in QP.LEProb length(top.i)!=1')
    le.not=le[-1]
    tc.not=tc[-top.i,-top.i,drop=FALSE]
    prob.top=p/n
    if (all(tc[,top.i]==0)) {ntc=nle(tc); prob.top=prob.top+(1-p)*nle(tc.not)/ntc}
    top.fac=q*prob.top*QP.LEProb(tc.not,le.not,p,q)
  } else {top.fac=0}
  
  if (q<1) {
    bot.i=which(as.numeric(rownames(tc))==le[n])
    if (length(bot.i)!=1) stop('err in QP.LEProb length(bot.i)!=1')
    le.nob=le[-n]
    tc.nob=tc[-bot.i,-bot.i,drop=FALSE]
    prob.bot=p/n
    if (all(tc[bot.i,]==0)) {if (is.na(ntc)) {ntc=nle(tc)}; prob.bot=prob.bot+(1-p)*nle(tc.nob)/ntc}
    bot.fac=(1-q)*prob.bot*QP.LEProb(tc.nob,le.nob,p,q)
  } else {bot.fac=0}
  
  return(top.fac+bot.fac)
}
