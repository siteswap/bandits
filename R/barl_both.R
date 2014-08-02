
library(hash)

L <- 2 # levers
# matrix of params
mp <- data.frame(aq=rep(5,L),bq=rep(95,L),ap=rep(1,L),bp=rep(1,L))

# Init hash
# h <- hash()

# Quality of lever l, with t trials ahead.
barl_both.q <- function(l,mp,t){
  
  # 1 step quality
  q1 <- mp[l,"aq"]/(mp[l,"bq"]+mp[l,"aq"]) * mp[l,"ap"]/(mp[l,"bp"]+mp[l,"ap"]) 
  
  if(t==0){
    
    return( q1 ) # Immediate value of lever l
    
  }else{
    # TODO - this code is dreadful
    
    # Quality no clicks
    mp$bp[l]  <- mp$bp[l] + 1
    k <- toString(c(mp,t-1))
    if(has.key(k,h)){
      sq <- unname(values(h,keys=k))
    }else{
      sq <- max(sapply(1:L,FUN=barl_both.q,mp=mp,t=t-1))
      h[k] <- sq
    }
    
    # Quality click, no acq
    mp$bp[l]  <- mp$bp[l] - 1
    mp$ap[l]  <- mp$ap[l] + 1
    mp$bq[l]  <- mp$bq[l] + 1
    k <- toString(c(mp,t-1))
    if(has.key(k,h)){
      fq <- unname(values(h,keys=k))
    }else{
      fq <- max(sapply(1:L,FUN=barl_both.q,mp=mp,t=t-1))
      h[k] <- fq
    }
    # prob of a click
    pc <- mp[l,"ap"]/(mp[l,"bp"]+mp[l,"ap"])
    
    # Quality acq
    mp$bq[l]  <- mp$bq[l] - 1
    mp$aq[l]  <- mp$aq[l] + 1
    k <- toString(c(mp,t-1))
    if(has.key(k,h)){
      aq <- unname(values(h,keys=k))
    }else{
      aq <- max(sapply(1:L,FUN=barl_both.q,mp=mp,t=t-1))
      h[k] <- aq
    }
    
    return( q1 + q1*aq + pc*fq + (1 - q1 - pc)*sq )
  }
}

barl_both.q.all <- function(mp,t){
  sapply(1:L,FUN=barl_both.q,mp=mp,t=t)/(t+1)
}


# clear(h)


