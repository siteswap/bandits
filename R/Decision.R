library(plyr)
library(hash)
library(R.matlab)

adtk.mab <- function(arms=5,campaignLen=1000,DFUN=adtk.ts_acqs,MFUN=adtk.m4,trueVals,prior) {
  
  qp <- trueVals$q*trueVals$p
  oracle <- max(qp)
  
  res <- data.frame(a=rep(0,arms),c=rep(0,arms),n=rep(0,arms))  # results
  armChoices <- rep(0,campaignLen)
  metrics <- data.frame(aqr=rep(0,campaignLen),regret=rep(0,campaignLen))
  
  mab <- list(arms=arms,armChoices=armChoices,res=res,trueVals=trueVals,round=0,len=campaignLen)
  adtk.mabplot( mab, method="arms" )
  
  for(round in 1:campaignLen){
    
    mab <- list(arms=arms,armChoices=armChoices,res=res,trueVals=trueVals,round=round,len=campaignLen,prior=prior)
    arm <- DFUN(mab)
    armChoices[round] <- arm
    
    res[arm,] <- res[arm,] + MFUN(trueVals[arm,]) # c(a,c,1)
    
    # Calc loss functions / KL div
    metrics[round,"aqr"] <- sum(res$a)/sum(res$n)
    metrics[round,"regret"] <- -1*mean(qp[armChoices[1:round]] - oracle)
    # Regret hould be less noisy than realized aqr
    
    mab <- list(arms=arms,armChoices=armChoices,res=res,trueVals=trueVals,round=round)
    if(0==(round %% 100)){
      adtk.mabplot( mab, method="arms" )
    }
  }
  
  mab <- list(arms=arms,armChoices=armChoices,
              res=res,trueVals=trueVals,round=campaignLen,metrics=metrics)
  
  return(mab)
}

# Documented model 4.
adtk.m4 <- function(vals){
  c <- rbinom(1,1,vals$p)
  a <- rbinom(1,c,vals$q)
  c(a,c,1)
}


# Thompson - choose according to how frequently arm is maximum.


# Achieve this through sampling rather than integral
adtk.ts <- function(arms,k,n,s1,s2){
  
  # Sample from posterior
  samples <- 100
  s <- matrix(rbeta(arms*samples,shape1=(s1 + k),shape2=(s2 + n - k)),nrow=arms)
  best <- apply(X=s,MARGIN=2,FUN=function(x){which(x==max(x))})
  df <- rbind(data.frame(g=best,c=1),data.frame(g=1:arms,c=0))
  probs <- ddply(df,~g,summarise,freq=sum(c))$freq/samples
  
  sample(1:arms,size=1,prob=probs)
}


adtk.ts_acqs <- function(mab){ 
  
  # Use method of moments to get estimate parameters
  # for equivalent beta distribution.
  pr <- mab$prior
  means <- c(pr$ap/(pr$ap + pr$bp), pr$aq/(pr$aq + pr$bq))
  moment2 <- c(beta(pr$ap + 2, pr$bp)/beta(pr$ap, pr$bp),
                beta(pr$aq + 2, pr$bq)/beta(pr$aq, pr$bq) )
  sig <- prod(moment2) - prod(means^2)
  mu <- prod(means)
  
  if(sig > mu*(1-mu) ){ warning("Method of moments estimation not valid") }

  alpha <- mu*(mu*(1-mu)/sig - 1)
  beta <- (1-mu)*(mu*(1-mu)/sig - 1)
  
  adtk.ts(mab$arms,mab$res$a,mab$res$n,alpha,beta) 
  } 

# TODO - needs further unit testing / validation.
adtk.gi_clicks <- function(mab){ 
  
  # Use method of moments to get estimate parameters
  # for equivalent beta distribution.
  pr <- mab$prior
  t <- mab$len - mab$round + 1
  a <- mab$res$a
  c <- mab$res$c
  n <- mab$res$n
  indices <- mapply(gi.v,a=pr$ap+c,b=pr$bq+n-c,t=t)
  equalBestArms <- which(indices==max(indices))
  if(length(equalBestArms)==1){
    return(equalBestArms)
  } else {
    return(sample(equalBestArms,size=1))
  }
} 

# Requires MATLAB server connection
# Create singleton
adtk.gi_both.isOpen <- FALSE # TODO - protect with exsits()
adtk.gi_both.matlab <- 0
adtk.gi_both <- function(mab){ 
  
  # Use method of moments to get estimate parameters
  # for equivalent beta distribution.
  pr <- mab$prior
  t <- mab$len - mab$round + 1
  a <- mab$res$a
  c <- mab$res$c
  n <- mab$res$n
  
  if(!adtk.gi_both.isOpen){
    Matlab$startServer()
    adtk.gi_both.matlab <- Matlab()
    adtk.gi_both.isOpen <- open(adtk.gi_both.matlab)
    evaluate(adtk.gi_both.matlab, "addpath '/home/alex/git/bandits/matlab' ")
    evaluate(adtk.gi_both.matlab, "global h;", "h = containers.Map;")
  }
  # close(matlab)
  
  M <- cbind(rep(t,length(n)),pr$ap+c,pr$bp+(n-c),pr$aq+a,pr$bq+(c-a))  
  setVariable(adtk.gi_both.matlab,M=M)
  evaluate(adtk.gi_both.matlab, "gi=GIBoth(M);")
  gi <- getVariable(adtk.gi_both.matlab,"gi")
  indices <- gi$gi
  
  equalBestArms <- which(indices==max(indices))
  if(length(equalBestArms)==1){
    return(equalBestArms)
  } else {
    return(sample(equalBestArms,size=1))
  }
} 

adtk.ts_clicks <- function(mab){ adtk.ts(mab$arms,mab$res$c,mab$res$n,mab$prior$ap,mab$prior$bp) }

adtk.ts_both <- function(mab){ # TODO - much duplicate code 
  
  # The 
  a <- mab$res$a
  c <- mab$res$c
  n <- mab$res$n
  
  # For each arm:
  samples <- 100
  arms <- mab$arms
  p <- rbeta(arms*samples,shape1=(c+mab$prior$ap),shape2=(n-c+mab$prior$bp))
  q <- rbeta(arms*samples,shape1=a+mab$prior$aq,shape2=(c-a+mab$prior$bq))
  pq <- p*q
  s <- matrix(pq,nrow=arms)
  # return vectors of probs where pq is higher of all arms.
  best <- apply(X=s,MARGIN=2,FUN=function(x){which(x==max(x))})
  df <- rbind(data.frame(g=best,c=1),data.frame(g=1:arms,c=0))
  probs <- ddply(df,~g,summarise,freq=sum(c))$freq/samples
  
  sample(1:arms,size=1,prob=probs)
}

# UCB1 - P. Auer, N. Cesa-Bianchi, and P. Fischer
# Experimentally, this explore too much and does not converge.
# Need to review derivation of rule.
adtk.ucb <- function(mab){
  
  # Action that maximizes xbar_j + sqrt(2*ln(n)/n_j)
  n <- mab$round
  kj <- mab$res$c # optimizing CTR
  nj <- mab$res$n
  
  # Pull each lever atleast once
  if(any(nj==0)){ return( which(nj==0)[1] )} 
  
  x <- kj/nj + sqrt(2*log(n)/nj)
  
  # Randomize choice when levers are equal
  equalBestArms <- which(x==max(x))
  if(length(equalBestArms)==1){
    r <- equalBestArms
  } else {
    r <- sample(equalBestArms,size=1)
  }
  # print(c(r,x))
  r
}

# Bayes adaptive RL
adtk.barl <- function(mab){
 
  t <- mab$len - mab$round
  betaVal <- mab$prior$bq + mab$res$n - mab$res$a
  alphaVal <- mab$prior$aq + mab$res$a
  mp <- data.frame(a=alphaVal,b=betaVal)
  valfun <- barl.q.all(mp,t)
  # Where arms have equal value, choose randomly
  # This means errors average out when taken repeatedly
  equalBestArms <- which(valfun==max(valfun))
  if(length(equalBestArms)==1){
    return(equalBestArms)
  } else {
    return(sample(equalBestArms,size=1))
  }
}

# Bayes adaptive RL
adtk.barl_both <- function(mab){
  
  t <- mab$len - mab$round
  betap <- mab$prior$bp + mab$res$n - mab$res$c
  alphap <- mab$prior$ap + mab$res$c
  betaq <- mab$prior$bq + mab$res$c - mab$res$a
  alphaq <- mab$prior$aq + mab$res$a
  mp <- data.frame(aq=alphaq,bq=betaq,ap=alphap,bp=betap)
  valfun <- barl_both.q.all(mp,t)
  # Where arms have equal value, choose randomly
  # This means errors average out when taken repeatedly
  equalBestArms <- which(valfun==max(valfun))
  if(length(equalBestArms)==1){
    return(equalBestArms)
  } else {
    return(sample(equalBestArms,size=1))
  }
}

adtk.PLOT_OFF <- TRUE

adtk.mabplot <- function(ts,method="aqr"){
  
  if(adtk.PLOT_OFF){return()}
  
  if(method=="aqr"){
    
    truev <- ts$trueVals$q * ts$trueVals$p
    plot(ts$metrics$aqr,type='l',ylim=c(0,max(truev)*1.2),ylab="mean achieved")
    
    abline(h=max(truev),col="red")  
    
  }else if (method=="arms"){
    
    nplot <- min(ts$arms,5)
    xvals <- seq(0,1,0.001)
    par(mfrow=c(nplot,1))
    for(a in 1:nplot){
      plot(x=xvals,y=dbeta(x=xvals,
                           shape1=1 + ts$res[a,"a"],
                           shape2=1 + ts$res[a,"n"] - ts$res[a,"a"]),
           type='l',
           ylab="density")
      truev <- ts$trueVals$q * ts$trueVals$p
      abline(v=truev[a],col="red")  
      if(a==1){ title(paste("First 5 arms, round: ",ts$round)) }
    }
    
  }
}


#########################
### Bayes Adaptive RL ###
#########################


# Init hash
barl.h <- hash()
# clear(barl.h)
# mp <- data.frame(a=rep(1,L),b=rep(1,L))

# Quality of lever l, with t trials ahead.
barl.q <- function(l,mp,t){
  
  h <- barl.h # create shorter reference
  L <- dim(mp)[1]
  # 1 step quality
  q1 <- mp[l,"a"]/(mp[l,"a"]+mp[l,"b"]) 
  
  if(t==0){
    
    return( q1 ) # Immediate value of lever l
    
  }else{
    
    # Quality under success
    mp$a[l]  <- mp$a[l] + 1
    # sq <- max(q(1,mp,t-1),q(2,mp,t-1))
    k <- toString(c(mp,t-1))
    if(has.key(k,h)){
      sq <- unname(values(h,keys=k))
    }else{
      sq <- max(sapply(1:L,FUN=barl.q ,mp=mp,t=t-1))
      h[k] <- sq
    }
    
    # Quality under failure
    mp$a[l]  <- mp$a[l] - 1
    mp$b[l]  <- mp$b[l] + 1
    k <- toString(c(mp,t-1))
    if(has.key(k,h)){
      fq <- unname(values(h,keys=k))
    }else{
      fq <- max(sapply(1:L,FUN=barl.q ,mp=mp,t=t-1))
      h[k] <- fq
    }
    
    return( q1 + q1*sq + (1 - q1)*fq )
  }
}

barl.q.all <- function(mp,t){
  L <- dim(mp)[1]
  sapply(1:L,FUN=barl.q ,mp=mp,t=t)/(t+1)
}


# Init hash
barl_both.h <- hash()
# mp <- data.frame(aq=rep(5,L),bq=rep(95,L),ap=rep(1,L),bp=rep(1,L))

# Quality of lever l, with t trials ahead.
barl_both.q <- function(l,mp,t){
  
  h <- barl_both.h # create shorter reference
  L <- dim(mp)[1]
  # 1 step quality
  q1 <- mp[l,"aq"]/(mp[l,"bq"]+mp[l,"aq"]) * mp[l,"ap"]/(mp[l,"bp"]+mp[l,"ap"]) 
  
  if(t==0){
    
    return( q1 ) # Immediate value of lever l
    
  }else{
    # TODO - tidy up and standardize with barl.q
    
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
    # prob of a click, no acq
    pc <- mp[l,"bq"]/(mp[l,"bq"]+mp[l,"aq"]) * mp[l,"ap"]/(mp[l,"bp"]+mp[l,"ap"]) 
    
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
  L <- dim(mp)[1]
  sapply(1:L,FUN=barl_both.q,mp=mp,t=t)/(t+1)
}

