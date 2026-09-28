
#exploring poset dimension
#GKN 12/12/24
#see https://doi.org/10.1016/0012-365X(73)90025-3 for some useful stuff
#there is actually julia code which does this and alot more - probably alot faster!
#https://juliapackages.com/p/posets

rm(list=ls())
#gc()

library(MASS)
library(mvtnorm)
library(mnem)    
library(igraph) 
library(Rgraphviz)
library(graph) 
library(coda)
library(lecount)
library(partitions)

#
setwd("C:/Users/nicholls.NICHOLLS2389/OneDrive - Nexus365/Documents/GitHub/Partial-order-HMM")

source(file="pofun.R")
source(file="modelfun.R")
source(file='vsp/vspfun.R')
source(file="dimfun.R")

#example
set.seed(103)
n=8; K=4
hc=rZPO(n,K,rRprior(),b=rep(0,n))
hr=transitive.reduction(hc); 
is.vsp(hc) #it isnt a VSP (ranseed 103)

showDAG(hc,main='transitive closure')   
showDAG(hr,main='transitive reduction')
#get all the LEs of hc and check their intersection gives us back the original PO hc
le<-le.expand(hc); dim(le)[2]; nle(hc); check.dec(hc,le)

#now get the dimension
dimension2(hc) #just the dimension - this example has dim 2 but isnt a VSP
(hc.decompose<-decompose2(hc)) #the decomposition - actually doesnt handle the case when hc is a VSP of dim=2 (but gets dim correct)
check.dec(hc,hc.decompose)
dimension3(hc)

#another example showing fast handling of floats
hc=matrix(c(0,0,1,1,0,0, 0,0,0,1,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0),6,6,byrow = TRUE) #Forbidden subgraph plus 2 floats
rownames(hc)<-colnames(hc)<-1:6
showDAG(hc)  

dimension2(hc) #should be 2
(hc.decompose<-decompose2(hc)) #add the floats to the decomposition of the FSG
check.dec(hc,hc.decompose)
dimension3(hc)

#look at relation between dimension, depth and width n=8 is already very slow - about 8 mins n=8, N=100
N=100; n=8; K=floor(n/2); T=V=D=W=L=rep(NA,N); H=vector('list',N)
set.seed(101)
for (t in 1:N) {
  H[[t]]<-rZPO(n,K,runif(1),b=rep(0,n))
  #print(nle(H[[t]]))
  #hc.decompose<-decompose2(H[[t]])
  V[t]=is.vsp(H[[t]])
  print(c(t,T[t]<-unname(system.time({D[t]<-dimension3(H[[t]])})['user.self']),D[t])) #dim(hc.decompose)[2]
  L[t]=dagdepth(H[[t]])
  #W[t]=dagwidth(H[[t]])
}
#D3=D; T3=T; T3S=sum(T)
sum(T)
100*table(D)/N #percent different dimensions
ftable(D,L) #notice dimension tends to increase with width note bound D \le W from theory


#make a crown PO - this PO achieves the max dimension floor(n/2) - makes sense when n even 
n=8; m=floor(n/2)
s1=combn(1:m, m-1, simplify = FALSE) 
s2=combn(1:m, 1, simplify = FALSE) 
hc=matrix(0,n,n); rownames(hc)<-colnames(hc)<-1:n
for (i in 1:m) {
  for (j in 1:m) {
    hc[i,m+j]=0+all(s2[[j]] %in% s1[[i]])
  }
}
showDAG(hc)
(n.le<-nle(hc))
choose(n.le,floor(n/2)) #the true dim of crown is floor(n/2) so dim-alg must search over more than n choose n/2 sets of LEs for decomposition
if (n<8) decompose2(hc) #should get dim=3 for n=6 - n=8 crown > 11 billion sets of LEs of size 4 ( > 60 million of size 3) 

if (n==8) {
  #actually it is easy to work out the decomposition of the n=8 crown by hand!
  hc.decompose=matrix(c(1,2,3,5,4,6,7,8, 4,1,2,6,3,7,8,5, 3,4,1,7,2,8,5,6, 2,3,4,8,1,5,6,7),8,4)
  hc.decompose
  showDAG(intersect.TO(hc.decompose))
  check.dec(hc,hc.decompose)
}

#can we identify a big crown as a subgraph?
h=hc
rownames(h)<-colnames(h)<-sample(1:n) 
tcg<-graph_from_adjacency_matrix(hc,mode ="directed")
bg<-graph_from_adjacency_matrix(h,mode ="directed")
subgraph_isomorphic(bg,tcg,induced=TRUE)
#yes