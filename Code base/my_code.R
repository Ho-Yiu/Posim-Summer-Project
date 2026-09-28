#basic code building up the project

#check whether or not mc is a dag
check.dag <- function(mc){!any(mc+t(mc)==2)}

#convert a sequence into a transitive closure
seqtodag <- function(LE,n){
  my.transitive.closure(seq2dag(LE,n))
}

#given a list of mcs, find the average and make an binary matrix that has entry one iff the entry in the average is above p
consensus.po <- function(po,p=.5){
  b <- length(po); out <- po[[1]]; n <- dim(out)[1]
  if (b==1) {return(po[[1]])}
  for (i in 2:b){
    out <- out + po[[i]]
  }
  matrix(as.numeric(out/b>p),n,n)
}

#Simulate an LE of mc under the top observation model
rFsm <- function(mc,noise.prob=0){
  if (noise.prob < 0 || noise.prob > 1){
    cat(c("noise.prob=",noise.prob,"\n"))
    stop("noise.prob must be in between 0 and 1")
  }
  if (is.null(dim(mc))) {return(mc)}
  n <- dim(mc)[1]; indegree <- colSums(mc); active <- rep(TRUE,n); out <- numeric(n);
  label <- as.numeric(colnames(mc)); id <- 1:n
  top <- indegree == 0
  for (i in 1:n){
    if (runif(1)<noise.prob){candidates <- active}else{candidates <- active & top}
    chosen.id <- my.sample(id[candidates],1)
    out[i] <- label[chosen.id]
    active[chosen.id] <- FALSE
    successors <- mc[chosen.id, ]!=0
    indegree[successors] <- indegree[successors] - 1
    top <- indegree == 0
  }
  out
}

#Find the log p value of LE and mc under maxLE, LE is a sequence
lpFsm <- function(LE, mc, p=0) {
  if (nrow(mc) != ncol(mc)) {stop("mc must be square")}
  rn <- rownames(mc)
  cn <- colnames(mc)
  if (is.null(rn) || is.null(cn) || !identical(rn, cn)) {stop("mc must have identical row and column names")}
  
  LE_chr <- as.character(LE)
  if (!setequal(LE_chr, rn)) {stop("LE must be a permutation of the row/column names of mc")}
  
  # Map labels in LE to matrix positions
  idx <- match(LE_chr, rn)
  n <- length(idx)
  A <- mc != 0
  diag(A) <- FALSE
  indegree <- colSums(A)
  active <- rep(TRUE, n)
  lp <- 0
  
  for (k in seq_len(n)) {
    x <- idx[k]                 # matrix position for current label
    tops <- active & (indegree == 0)
    n.tops <- sum(tops)
    
    if (!tops[x]) {
      if (p==0){
        #if the model is noise free, return -Inf when next node doesn't belongs to tops
        return(-Inf)
      }else{
        #if not noise free and next node is not in tops, add this to lp
        lp <- lp + log(p/(n-k+1))
      }
    }else{
      #if next node is in tops, add this to lp regardless of noise
      lp <- lp + log(p/(n-k+1) + (1-p)/n.tops)
    }
    
    active[x] <- FALSE
    successors <- A[x, ] & active
    indegree[successors] <- indegree[successors] - 1 #save time by only doing the calculation on the sucessors since the rest of indegree are unchanged 
  }
  
  lp
}

#given a adjacent matrix of a sequence, this function returns a vector out of it
dag2seq <- function(le){
  n <- nrow(le)
  v <- colSums(le)
  r <- numeric(n)
  
  for (j in 0:(n - 1)) {
    idx <- which(v == j)
    if (length(idx) == 1) {
      r[j + 1] <- idx
    } else {
      stop("dag2seq expects exactly one vertex per level; got multiple.")
    }
  }
  r
}

#give you the index of element of y wrt the element in x
which2 <- function(x,y){
  n <- dim(x)
  v <- numeric(n)
  for (i in 1:n){v[i] <- which(names(x) == names(y)[i])}
  x[v]
}

#check if 2 mcs are consistent
consistent <- function(mc1,mc2){!any(mc1+t(mc2) == 2)}

#This finds the maximal elements of an mc
findtop <- function(mc){
  if(is.null(dim(mc))){
    1
  }else{
    which(colSums(mc) == 0)
  }
}

#This is the same as sample except it will only return the element given when length(x)==1
my.sample <- function (x, size, replace = FALSE, prob = NULL) {
  if (length(x)==1) {return(x)} else{return(sample(x, size, replace, prob))}
}

#do I still use this code???
bubble_step <- function(le, tc) {
  n <- length(le)
  if (n < 2L || runif(1) < 0.5) return(le)
  i <- sample.int(n - 1L, 1L); a <- le[i]; b <- le[i + 1L]
  
  # Swap only when a and b are incomparable.
  if (tc[a, b] == 0 && tc[b, a] == 0) le[c(i, i + 1L)] <- le[c(i + 1L, i)]
  le
}

#generate an LE of mc uniformly at random 
runifLE2 <- function(tc, N, burn = 1000L, thin = 10L, le = rmaxLE(tc)) {
  for (i in seq_len(burn))
    le <- bubble_step(le, tc)
  
  out <- matrix(NA_integer_, N, length(le))
  
  for (i in seq_len(N)) {
    for (j in seq_len(thin)) le <- bubble_step(le, tc)
    out[i, ] <- le
  }
  return(out)
}

#observe a bucket order with the geometric observation model
rgeom.bucket <- function(LE, p=.5) {
  if (!is.numeric(p) || length(p) != 1L || is.na(p) || p <= 0 || p >= 1)
    stop("p must be a single numeric value strictly between 0 and 1")
  n <- length(LE)
  if (n==0) return(list())
  group <- cumsum(c(TRUE, runif(n - 1L) < p))
  return(unname(split(LE, group)))
}

#given a bucket order b, return an LE of b uniformly at random
runif.bucket.LE <- function(b,size=1) {
  if (!is.list(b)) stop("Argument b of runif.bucket.LE is not a list")
  lens <- lengths(b); out <- vector("list",size)
  for (i in seq_len(size)){
    v <- integer(sum(lens)); pos <- 1L
    for (x in b) {
      n <- length(x)
      if (n == 0L) next
      v[pos:(pos + n - 1L)] <- x[sample.int(n)]
      pos <- pos + n
    }
    out[[i]] <- v
  }
  return(out)
}

#This function generate n.samples bucket order out of a given mc according to the geometric model
rgeom.bucket.from.po <- function(mc,p=.5,n.samples=1,draw=TRUE,suborder=FALSE,noise.prob=0){
  if (is.null(mc)) return(rep(list(c(1)),n.samples))
  n <- dim(mc)[1]; out <- vector("list",n.samples)
  if (!suborder) {
    out <- lapply(out,function(x){rgeom.bucket(rFsm(mc=mc,noise.prob=noise.prob),p)})
  }else{
    out <- lapply(out,function(x){rgeom.bucket(rFsm(mc=mc,noise.prob=noise.prob)[sort(sample(1:n,sample(2:n,1)))],p)})
  }
  out
}

#This code checks the depth of a dag
dag.depth <- function(mc,n){
  if (!all(dim(mc)==n)) {stop("dimension of mc in dag.depth in not the same as n")}
  if (!check.dag(mc)) {stop("mc in dag.depth has to be a dag")}
  indegree <- colSums(mc); active <- rep(T,n); out <- 0
  while(sum(active)>0){
    tops <- which(indegree==0 & active)
    indegree <- indegree - colSums(mc[tops, , drop=FALSE])
    active[tops] <- FALSE
    out <- out + 1
  }
  out
}

#This function finds the hamming distance between 2 mcs
hamming.dist <- function(mc1,mc2){sum(abs(mc1-mc2))}

#mcmc code
rposterior.bucket <- function(n=0L,b,n.samples=1,burn=1000,thin=burn,draw=TRUE,m=matrix(0,n,n,dimnames = list(1:n,1:n)),return.lsm=F,r=.5, weights = rep(0,n), noise.prob=0){
  #weights: ratio of the logarithm of no of posets there are of different depth, p: the p in the noise model
  
  if (class(b)!="list") stop("Argument b in rposterior.bucket2 must be a list")
  if (class(n)!="integer" | !(n > 1)) stop("Argument n of rposterior.bucket2 must be an integer greater than 1")
  if (length(weights) != n) stop("Argument n of rposterior.bucket2 must be a vector of length n")
  
  #n <- 5L; n.samples=3; burn <- 1000; thin <- burn; draw <- TRUE
  #a <- 0
  #initialization
  B <- length(b); Y <- pos <- vector("list",B); M <- numeric(B); out <- vector("list",n.samples)
  mc <- my.transitive.closure(m); mr <- transitive.reduction(m)
  mp <- m; mpc <- mc; mpr <- mr; m.depth <- mp.depth <- 1
  mcmr <- sum(mc - mr)*log(2)
  lsm <- lsmp <- numeric(B)
  c <- 1 #counter for counting how many mc has been recorded for the output
  a <- 0; aa <- 0 #counter for the number of accepted proposals, and number of accepted matrices proposal respectively
  
  for (i in seq_len(B)){
    Y[[i]] <- runif.bucket.LE(b[[i]],1)[[1]]
    pos[[i]] <- sort(Y[[i]])
    M[i] <- length(Y[[i]])
    lsm[i] <- lpFsm(Y[[i]],mc[pos[[i]],pos[[i]]],noise.prob)
  }
  
  #Don't touch - this is the progress bar!
  pb<-txtProgressBar(min=0, max=1, style=3)
  
  t <- burn + thin*n.samples
  
  if (return.lsm) {like.trace <- numeric(n.samples)}
  
  #mcmc
  for (i in seq_len(burn + thin*n.samples)){
    mp <- m; mpc <- mc; mpr <- mr; lsmp <- lsm; mcmrp <- mcmr; mp.depth <- m.depth
    
    if (runif(1)<r){
      u <- sample(1:n,2,replace=FALSE)
      
      if (!(mc[u[1],u[2]] == mr[u[1],u[2]])){
        aa <- aa+1
        }else{
        
        mp[u[1],u[2]] <- 1-m[u[1],u[2]]
        mpc <- my.transitive.closure(mp)
        
        if(check.dag(mpc)){
          cons <- TRUE
          j <- 1
          
          if (noise.prob==0){
            while(j<=B && cons){
              cons <- consistent(unname(mpc[pos[[j]],pos[[j]]]),seqtodag(rank(Y[[j]]),M[j]))
              j <- j+1
            }
          }
          
          if (cons){
            for (k in 1:B){
              ry <- sort(Y[[k]])
              lsmp[k] <- lpFsm(Y[[k]],mpc[ry,ry],noise.prob)
            }
            mpr <- transitive.reduction(mpc)
            mcmrp <- (sum(mpc - mpr))*log(2)
            mp.depth <- dag.depth(mpc,n)
            
            if (log(runif(1)) < sum(lsmp) - sum(lsm) + mcmr - mcmrp + weights[m.depth] - weights[mp.depth]){
              m <- mp; mc <- mpc; mr <- mpr; lsm <- lsmp; mcmr <- mcmrp; m.depth <- mp.depth
              a <- a+1; aa <- aa+1
            }
          }
        }
      }
    }else{
      
      y.id <- my.sample(1:B,1)
      yp <- runif.bucket.LE(b[[y.id]],1)[[1]]
      lsmp[y.id] <- lpFsm(yp,mc[pos[[y.id]],pos[[y.id]]],noise.prob)
      
      if (noise.prob==0){
        if (consistent(mc[pos[[y.id]],pos[[y.id]]],seqtodag(rank(yp),M[y.id]))){
          if (log(runif(1)) < sum(lsmp) - sum(lsm)){
            Y[[y.id]] <- yp
            pos[[y.id]] <- sort(yp)
            lsm[y.id] <- lsmp[y.id]
            a <- a+1
          }
        }
      }else{
        if (log(runif(1)) < sum(lsmp) - sum(lsm)){
          Y[[y.id]] <- yp
          pos[[y.id]] <- sort(yp)
          lsm[y.id] <- lsmp[y.id]
          a <- a+1
        }
      }
      
    }
    
    #recording the mcmc
    if (i > burn && (i-burn)%%thin==0){
      if (return.lsm) {like.trace[c] <- sum(lsm)}
      out[[c]] <- matrix(as.integer(mc),n,n,dimnames=dimnames(m)); c <- c+1
      if (draw){
        showDAG(transitive.reduction(mc),vertex.color = "white", vertex.size = 20, edge.arrow.size = 0.5)
        title(paste(c(i,"th sample from posterior"),collapse=""))
      }
    }
    
    
    if (i%%max(1,floor(t/100))==0) setTxtProgressBar(pb, i/(t))
    
  }
  
  close(pb)
  
  if (return.lsm){
    list(samples = out, acceptance = c(a/t,aa/t,aa/a), trace = like.trace, weights = weights)
  }else{
    list(samples = out, acceptance = c(a/t,aa/t,aa/a), weights = weights)
  }
  
}

#mcmc code
rposterior.bucket.mean <- function(n=0L,b,n.samples=1,burn=1000,thin=burn,draw=TRUE,m=matrix(0,n,n,dimnames = list(1:n,1:n)),return.lsm=F,r=.5, weights = rep(0,n), noise.prob=0){
  #weights: ratio of the logarithm of no of posets there are of different depth, p: the p in the noise model
  #this only returns the mean partial order for the q-consensus
  
  if (class(b)!="list") stop("Argument b in rposterior.bucket2 must be a list")
  if (class(n)!="integer" | !(n > 1)) stop("Argument n of rposterior.bucket2 must be an integer greater than 1")
  if (length(weights) != n) stop("Argument n of rposterior.bucket2 must be a vector of length n")
  
  #n <- 5L; n.samples=3; burn <- 1000; thin <- burn; draw <- TRUE
  #a <- 0
  #initialization
  B <- length(b); Y <- pos <- vector("list",B); M <- numeric(B); out <- matrix(0,n,n,dimnames=list(1:n,1:n))
  mc <- my.transitive.closure(m); mr <- transitive.reduction(m)
  mp <- m; mpc <- mc; mpr <- mr; m.depth <- mp.depth <- 1
  mcmr <- sum(mc - mr)*log(2)
  lsm <- lsmp <- numeric(B)
  c <- 1 #counter for counting how many mc has been recorded for the output
  a <- 0; aa <- 0 #counter for the number of accepted proposals, and number of accepted matrices proposal respectively
  
  for (i in seq_len(B)){
    Y[[i]] <- runif.bucket.LE(b[[i]],1)[[1]]
    pos[[i]] <- sort(Y[[i]])
    M[i] <- length(Y[[i]])
    lsm[i] <- lpFsm(Y[[i]],mc[pos[[i]],pos[[i]]],noise.prob)
  }
  
  #Don't touch - this is the progress bar!
  pb<-txtProgressBar(min=0, max=1, style=3)
  
  t <- burn + thin*n.samples
  
  #if (return.lsm) {like.trace <- numeric(n.samples)}
  
  #mcmc
  for (i in seq_len(burn + thin*n.samples)){
    mp <- m; mpc <- mc; mpr <- mr; lsmp <- lsm; mcmrp <- mcmr; mp.depth <- m.depth
    
    if (runif(1)<r){
      u <- sample(1:n,2,replace=FALSE)
      
      if (!(mc[u[1],u[2]] == mr[u[1],u[2]])){
        aa <- aa+1
      }else{
        
        mp[u[1],u[2]] <- 1-m[u[1],u[2]]
        mpc <- my.transitive.closure(mp)
        
        if(check.dag(mpc)){
          cons <- TRUE
          j <- 1
          
          if (noise.prob==0){
            while(j<=B && cons){
              cons <- consistent(unname(mpc[pos[[j]],pos[[j]]]),seqtodag(rank(Y[[j]]),M[j]))
              j <- j+1
            }
          }
          
          if (cons){
            for (k in 1:B){
              ry <- sort(Y[[k]])
              lsmp[k] <- lpFsm(Y[[k]],mpc[ry,ry],noise.prob)
            }
            mpr <- transitive.reduction(mpc)
            mcmrp <- (sum(mpc - mpr))*log(2)
            mp.depth <- dag.depth(mpc,n)
            
            if (log(runif(1)) < sum(lsmp) - sum(lsm) + mcmr - mcmrp + weights[m.depth] - weights[mp.depth]){
              m <- mp; mc <- mpc; mr <- mpr; lsm <- lsmp; mcmr <- mcmrp; m.depth <- mp.depth
              a <- a+1; aa <- aa+1
            }
          }
        }
      }
    }else{
      
      y.id <- my.sample(1:B,1)
      yp <- runif.bucket.LE(b[[y.id]],1)[[1]]
      lsmp[y.id] <- lpFsm(yp,mc[pos[[y.id]],pos[[y.id]]],noise.prob)
      
      if (noise.prob==0){
        if (consistent(mc[pos[[y.id]],pos[[y.id]]],seqtodag(rank(yp),M[y.id]))){
          if (log(runif(1)) < sum(lsmp) - sum(lsm)){
            Y[[y.id]] <- yp
            pos[[y.id]] <- sort(yp)
            lsm[y.id] <- lsmp[y.id]
            a <- a+1
          }
        }
      }else{
        if (log(runif(1)) < sum(lsmp) - sum(lsm)){
          Y[[y.id]] <- yp
          pos[[y.id]] <- sort(yp)
          lsm[y.id] <- lsmp[y.id]
          a <- a+1
        }
      }
      
    }
    
    #recording the mcmc
    if (i > burn && (i-burn)%%thin==0){
      #if (return.lsm) {like.trace[c] <- sum(lsm)}
      out <- out+mc; c <- c+1
      #if (draw){
      # showDAG(transitive.reduction(mc),vertex.color = "white", vertex.size = 20, edge.arrow.size = 0.5)
      #title(paste(c(i,"th sample from posterior"),collapse=""))
      #}
    }
    
    
    if (i%%max(1,floor(t/100))==0) setTxtProgressBar(pb, i/(t))
    
  }
  
  close(pb)
  
  list(mean = out/n.samples, acceptance = c(a/t,aa/t,aa/a), weights = weights)
  
}

my.rupo <- function(n,h,n.samples=1L,burn=1000,thin=burn){
  
  if (!all(dim(h)==n)) {stop("dimension of h in my.rupo is not the same as n")}
  if (!check.dag(h)) {stop("h in my.rupo must be a dag")}
  
  #initial state
  m <- matrix(0,n,n,dimnames = list(1:n,1:n))
  mc <- my.transitive.closure(m); mr <- transitive.reduction(m)
  mcmr <- 0
  out <- vector("list",n.samples); outed <- 0
  t <- burn + n.samples*thin
  accepted <- 0
  
  pb<-txtProgressBar(min=0, max=1, style=3)
  
  #mcmc
  for (i in 1:t){
    v <- sample(1:n,2)
    if (mc[v[1],v[2]] != mr[v[1],v[2]]){
      #if we toggled an edge that does not change the colsure, accept
    }else{
      mp <- m; mp[v[1],v[2]] <- 1-mp[v[1],v[2]]
      mpc <- my.transitive.closure(mp); mpr <- transitive.reduction(mp)
      if (check.dag(mpc)){
        if (consistent(mpc,h)){
          mcmrp <- sum(mpc - mpr)
          if (log(runif(1)) < (mcmr - mcmrp)*log(2)){
            m <- mp; mc <- mpc; mr <- mpr; mcmr <- mcmrp
            accepted <- accepted + 1
          }
        }
      }
    }
    
    if (i>burn && (i-burn)%%thin==0){
      outed <- outed + 1
      out[[outed]] <- matrix(as.integer(mc),n,n,dimnames=list(1:n,1:n))
    }
    
    setTxtProgressBar(pb, i/(t))
    
  }
  
  return(list(
    samples = out,
    acceptance = accepted/t
  ))
  
}

bucket2mc <- function(b,n){
  if (!is.list(b)) stop("b must be a bucket order")
  
  m <- matrix(0,n,n,dimnames=list(1:n,1:n))
  B <- length(b)
  
  if (B >= 2){
    for (i in 1:(B-1)){
      m[b[[i]],unlist(b[(i+1):B])] <- 1
    }
  }
  
  m
}

intersection.po <- function(po){
  if (!is.list(po)) stop("po must be a list of matrices")
  if (length(po)==0) stop("po must contain at least one matrix")
  
  out <- po[[1]]
  
  if (length(po)>=2){
    for (i in 2:length(po)){
      out <- out * po[[i]]
    }
  }
  
  out
}

bucket2mc.full <- function(b,n){
  if (!is.list(b)) stop("b must be a list of bucket orders")
  
  out <- vector("list",length(b))
  
  for (i in 1:length(b)){
    m <- matrix(0,n,n,dimnames=list(1:n,1:n))
    B <- length(b[[i]])
    
    if (B >= 2){
      for (j in 1:(B-1)){
        m[b[[i]][[j]],unlist(b[[i]][(j+1):B])] <- 1
      }
    }
    
    out[[i]] <- m
  }
  
  out
}

mean.po <- function(list.po){
  if (!is.list(list.po)) stop("list.po has to be a list in list.po")
  len <- length(list.po)
  s <- list.po[[1]]
  for (i in 2:len){
    s <- s + list.po[[i]]
  }
  s/len
}

consensus.po2 <- function(mc,p=1/3){
  n <- dim(mc)[1]
  for (i in 1:n){
    for (j in 1:n){
      mc[i,j] <- (mc[i,j] > p)
    }
  }
  my.transitive.closure(mc)
}

hamming.dist2 <- function(mc1,mc2){sum(abs(mc1-mc2)^2)}