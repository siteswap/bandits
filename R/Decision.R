

adtk.mab <- function(arms=5,campaignLen=1000,a=6,b=10,DFUN=adtk.ts) {
  
  trueVals <- rbeta(arms,shape1=a,shape2=b)       # levers
  res <- data.frame(a=rep(0,arms),c=rep(0,arms),n=rep(0,arms))  # results
  armChoices <- rep(0,campaignLen)
  loss <- data.frame(loss1=rep(0,campaignLen))
  
  mab <- list(arms=arms,armChoices=armChoices,res=res,trueVals=trueVals,round=0)
  adtk.mabplot( mab )
  
  for(round in 1:campaignLen){
  
    mab <- list(arms=arms,armChoices=armChoices,res=res,trueVals=trueVals,round=round)
    arm <- DFUN(mab)
    armChoices[round] <- arm
    res[arm,"a"] <- res[arm,"a"] + rbinom(1,size=1,prob=trueVals[arm])
    res[arm,"n"] <- res[arm,"n"] + 1
    
    # Calc loss functions / KL div
    loss[round,"loss1"] <- sum(res$a)/sum(res$n)

    mab <- list(arms=arms,armChoices=armChoices,res=res,trueVals=trueVals,round=round)
    if(0==(round %% 100)){
      adtk.mabplot( mab )
    }
  }
  
  mab <- list(arms=arms,armChoices=armChoices,
              res=res,trueVals=trueVals,round=campaignLen,loss=loss)
  
  return(mab)
}

# Thompson - choose according to how frequently arm is maximum.
# Achieve this through sampling rather than integral
adtk.ts <- function(mab){
  
  # Sample from posterior
  res <- mab$res
  arms <- mab$arms
  samples <- 100
  s <- matrix(rbeta(arms*samples,shape1=(1 + res$a),shape2=(1 + res$n - res$a)),nrow=arms)
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

adtk.mabplot <- function(ts,method=""){
  
  if(method=="loss"){
    
    par(mfrow=c(1,1))
    plot(ts$loss$loss1,type='l',ylim=c(0,1),ylab="mean achieved")
    abline(h=max(ts$trueVals),col="red")  
    
  }else{
  
    nplot <- min(ts$arms,5)
    xvals <- seq(0,1,0.001)
    par(mfrow=c(nplot,1))
    for(a in 1:nplot){
      plot(x=xvals,y=dbeta(x=xvals,
                           shape1=1 + ts$res[a,"a"],
                           shape2=1 + ts$res[a,"n"] - ts$res[a,"a"]),
           type='l',
           ylab="density")
      abline(v=ts$trueVals[a],col="red")  
      if(a==1){ title(paste("First 5 arms, round: ",ts$round)) }
    }
  
  }
}


# 1. one only tries to maximise number of clicks, ignoring conversions

# 2. one only tries to maximise number of conversions, ignoring clicks

# 3. one tries to maximise number of conversions, considering clicks as well, assuming however, 
#     that the two rates are independent (formulate the model carefully and extract a confidence interval for UCB). 

