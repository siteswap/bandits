library(plyr)


adtk.mab <- function(arms=5,campaignLen=1000,DFUN=adtk.ts_acqs,MFUN=adtk.m4,trueVals=c()) {
  
  if( length(trueVals)==0 ){ 
    trueVals <- data.frame(q=rbeta(arms,shape1=5,shape2=10),
                           p=rbeta(arms,shape1=5,shape2=5))
  }  
  
  qp <- trueVals$q*trueVals$p
  oracle <- max(qp)
  
  res <- data.frame(a=rep(0,arms),c=rep(0,arms),n=rep(0,arms))  # results
  armChoices <- rep(0,campaignLen)
  metrics <- data.frame(aqr=rep(0,campaignLen),regret=rep(0,campaignLen))
  
  mab <- list(arms=arms,armChoices=armChoices,res=res,trueVals=trueVals,round=0,len=campaignLen)
  adtk.mabplot( mab, method="arms" )
  
  for(round in 1:campaignLen){
    
    mab <- list(arms=arms,armChoices=armChoices,res=res,trueVals=trueVals,round=round,len=campaignLen)
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
adtk.ts <- function(arms,k,n){
  
  # Sample from posterior
  samples <- 100
  s <- matrix(rbeta(arms*samples,shape1=(1 + k),shape2=(1 + n - k)),nrow=arms)
  best <- apply(X=s,MARGIN=2,FUN=function(x){which(x==max(x))})
  df <- rbind(data.frame(g=best,c=1),data.frame(g=1:arms,c=0))
  probs <- ddply(df,~g,summarise,freq=sum(c))$freq/samples
  
  sample(1:arms,size=1,prob=probs)
}

adtk.ts_acqs <- function(mab){ adtk.ts(mab$arms,mab$res$a,mab$res$n) }
adtk.ts_clicks <- function(mab){ adtk.ts(mab$arms,mab$res$c,mab$res$n) }


adtk.ts_both <- function(mab){ # TODO - much duplicate code 
  
  # The 
  a <- mab$res$a
  c <- mab$res$c
  n <- mab$res$n
  
  # For each arm:
  samples <- 100
  p <- rbeta(arms*samples,shape1=(c+1),shape2=(n-c+1)) # Flat prior on clicks
  s1 <- 5   # Strong prior on q TODO - parameterize this
  s2 <- 100 # Strong prior on q TODO - parameterize this
  q <- rbeta(arms*samples,shape1=a+s1,shape2=(c-a+s2))
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
  kj <- mab$res$a
  nj <- mab$res$n
  
  # Pull each lever atleast once
  if(any(nj==0)){ return( which(nj==0)[1] )} 
  
  x <- kj/nj + sqrt(2*log(n)/nj)
  
  # Randomize choie when levers are equal
  sample(which(x==max(x)),size=1) 
}

# Bayes adaptive RL
adtk.barl <- function(mab){
 
  
  t <- mab$len - mab$round
  betaVal <- 1+ mab$res$n - mab$res$a
  alphaVal <- 1 + mab$res$a
  mp <- data.frame(a=alphaVal,b=betaVal)
  valfun <- q.all(mp,t)
  # Where arms have equal value, choose randomly
  # This means errors average out when taken repeatedly
  sample(which(valfun==max(valfun)),size=1)
}

adtk.mabplot <- function(ts,method="aqr"){
  
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


