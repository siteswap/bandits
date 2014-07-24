
adtk.mab <- function (arms=5,campaignLen=1000,a=6,b=10,decision="ts") {
  
  trueVals <- rbeta(arms,shape1=a,shape2=b)       # levers
  pst <- data.frame(a=rep(1,arms),b=rep(1,arms))  # posteriors
  armChoices <- rep(0,campaignLen)
  loss <- data.frame(loss1=rep(0,campaignLen))
  
  adtk.mabplot( list(arms=arms,armChoices=armChoices,pst=pst,trueVals=trueVals,round=0) )
  
  for(round in 1:campaignLen){
  
    if(decision=="ts"){
      arm <- adtk.ts(pst,arms)
    }
    
    armChoices[round] <- arm
    result <- rbinom(1,size=1,prob=trueVals[arm])
    
    # Update posteriors (assuming m2)
    pst[arm,"a"]  <- pst[arm,"a"] + result
    pst[arm,"b"]  <- pst[arm,"b"] + (!result)*1
    
    # Calc loss functions / KL div
    loss[round,"loss1"] <- (sum(pst$a) - arms)/(sum(pst$a) + sum(pst$b) - 2*arms)

    if(0==(round %% 100)){
      adtk.mabplot( list(arms=arms,armChoices=armChoices,pst=pst,trueVals=trueVals,round=round) )
    }
  }
  
  list(arms=arms,armChoices=armChoices,
        pst=pst,trueVals=trueVals,round=campaignLen,loss=loss)
}

# Thompson - choose according to how frequently arm is maximum.
# Achieve this through sampling rather than integral
adtk.ts <- function(pst,arms){
  
  # Sample from posterior
  samples <- 100
  s <- matrix(rbeta(arms*samples,shape1=pst$a,shape2=pst$b),nrow=arms)
  best <- apply(X=s,MARGIN=2,FUN=function(x){which(x==max(x))})
  df <- rbind(data.frame(g=best,c=1),data.frame(g=1:arms,c=0))
  probs <- ddply(df,~g,summarise,freq=sum(c)/samples)$freq
  
  sample(1:arms,size=1,prob=probs)
}

adtk.mabplot <- function(ts,method=""){
  
  if(method=="loss"){
    
    par(mfrow=c(1,1))
    plot(ts$loss$loss1,type='l',ylim=c(0,1))
    abline(h=max(ts$trueVals),col="red")  
    
  }else{
  
    nplot <- min(ts$arms,5)
    xvals <- seq(0,1,0.001)
    par(mfrow=c(nplot,1))
    for(a in 1:nplot){
      plot(x=xvals,y=dbeta(x=xvals,
                           shape1=ts$pst[a,"a"],
                           shape2=ts$pst[a,"b"]),
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

