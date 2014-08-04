
library(reshape)
library(ggplot2)

compareStrats <- function(arms,initVals,campaignLen,trials,strats,strat_names,prior){
  
  results <- array(0,c(length(strats),campaignLen,trials))
  
  for(t in 1:trials){  
    for(s in 1:length(strats)){
      F <- strats[[s]]
      results[s,,t] <- adtk.mab(DFUN=F,trueVals=initVals,campaignLen=campaignLen,arms=arms,prior=prior)$metrics$regret
    }
  }
  
  # Careful of joint winners
  # winners <- apply(X=results[,100,],MARGIN=2,function(x){ which(x==max(x)) }) 
  
  means <- apply(X=results,MARGIN=c(1,2),FUN=mean)
  rownames(means) <- strat_names
  sdevs <- apply(X=results,MARGIN=c(1,2),FUN=sd)
  rownames(sdevs) <- strat_names
  
  
  df <- melt(means,id=c())
  df$se <- melt(sdevs,id=c())$value
  return(df)
}

#############
# Validation 1 - set p==1
#############

# With pin-point distribution on acquisition rate, all algos are doing the same thing - learning click rate.

arms <- 2
prior <- data.frame(aq=1000,bq=1,ap=1,bp=1) 
initVals <- data.frame(q=rbeta(arms,shape1=prior$aq,shape2=prior$bq),p=rbeta(arms,shape1=prior$ap,shape2=prior$bp))
campaignLen <- 200
trials <- 10
strats <- c(adtk.ts_clicks,adtk.ts_acqs,adtk.ts_both) # Careful, takes the values as the values at time of assignment
strat_names <- c("ts_clicks","ts_acqs","ts_both")


df <- compareStrats(arms,initVals,campaignLen,trials,strats,strat_names,prior)

ggplot(data=df, aes(x=X2,y=value,colour=X1)) + 
  geom_line() +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1, alpha=.5) +
  ggtitle("Validation \n (p==1)")


#############
# Validation 2 - set q==1
#############

# Clicks should never learn
# acqs and both should be equivalent

prior <- data.frame(aq=1,bq=1,ap=1000,bp=1) 
initVals <- data.frame(q=rbeta(arms,shape1=prior$aq,shape2=prior$bq),p=rbeta(arms,shape1=prior$ap,shape2=prior$bp))

df <- compareStrats(arms,initVals,campaignLen,trials,strats,strat_names,prior)

ggplot(data=df, aes(x=X2,y=value,colour=X1)) + 
  geom_line() +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1, alpha=.5) +
  ggtitle("Validation \n (q==1)")


###############
# Experiment 1 - compare Thompson Sampling of clicks vs acquisitions vs both
###############

# Remember - clicks will learn the best p fast.
# When this is not same as best r, it will never converge, 
# but when it does, it will look really good. 
# Should test over many values of init parameters.

arms <- 10
prior <- data.frame(aq=10,bq=15,ap=1,bp=1) 
initVals <- data.frame(q=rbeta(arms,shape1=prior$aq,shape2=prior$bq),p=rbeta(arms,shape1=prior$ap,shape2=prior$bp))
campaignLen <- 300
trials <- 10
strats <- c(adtk.ts_clicks,adtk.ts_acqs,adtk.ts_both) # Careful, takes the values as the values at time of assignment
strat_names <- c("ts_clicks","ts_acqs","ts_both")


df <- compareStrats(arms,initVals,campaignLen,trials,strats,strat_names,prior)

ggplot(data=df, aes(x=X2,y=value,colour=X1)) + 
  geom_line() +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1, alpha=.5) +
  ggtitle( paste("Thompson Sampling \n clicks vs acquisitions vs both \n r:",
                 toString(round(initVals$q*initVals$p,2)),
                 "\n p:",toString(round(initVals$p,2)),
                 "\n q:",toString(round(initVals$q,2))
                 ))




###############
# Experiment 2 - Thompson Sampling vs Bayes Adaptive (acquisitions only)
###############

arms <- 2 # Stick to 2 arms for these short campaigns
prior <- data.frame(aq=1,bq=1,ap=1000,bp=1) 
initVals <- data.frame(q=rbeta(arms,shape1=prior$aq,shape2=prior$bq),p=1) 
campaignLen <- 10
trials <- 10
strats <- c(adtk.ts_acqs,adtk.barl)
strat_names <- c("ts_acqs","ba_acqs")


df <- compareStrats(arms,initVals,campaignLen,trials,strats,strat_names,prior)

ggplot(data=df, aes(x=X2,y=value,colour=X1)) + 
  geom_line() +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1, alpha=.5) +
  ggtitle( paste("Thompson Sampling vs Bayes Adaptive \n (acquisitions only) \n r:",
               toString(round(initVals$q*initVals$p,2)),
               "\n p:",toString(round(initVals$p,2)),
               "\n q:",toString(round(initVals$q,2))
                ))



###############
# Experiment 3 - Thompson Sampling vs Bayes Adaptive (both)
###############

# REMEMBER - clicks will learn the best p fast.
# When this is not same as best r, it will never converge, 
# but when it does, it will look really good. 
# Should test over many values of init parameters.


arms <- 2 # Stick to 2 arms for these short campaigns
prior <- data.frame(aq=3,bq=6,ap=1,bp=1) 
initVals <- data.frame(q=rbeta(arms,shape1=prior$aq,shape2=prior$bq),p=rbeta(arms,shape1=prior$ap,shape2=prior$bp))
campaignLen <- 10
trials <- 100
strats <- c(adtk.ts_acqs,adtk.ts_both,adtk.barl,adtk.barl_both)
strat_names <- c("ts_acqs","ts_both","ba_acqs","ba_both")

df <- compareStrats(arms,initVals,campaignLen,trials,strats,strat_names,prior)

ggplot(data=df, aes(x=X2,y=value,colour=X1)) + 
  geom_line() +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1, alpha=.5) +
  ggtitle( paste("Thompson Sampling vs Bayes Adaptive \n (both) \n r:",
               toString(round(initVals$q*initVals$p,2)),
               "\n p:",toString(round(initVals$p,2)),
               "\n q:",toString(round(initVals$q,2))
               ))



##########################################
### Illustration of Bayes optimal rule ###
##########################################

barl.illustration <- function(){

  # Experimentally, the performance of an arm does not
  # impact the relative value of other arms.
  # If the best arm is very good, it has the effect of 
  # 'stretching out' t as it requires several rounds of failure
  # before other arms are played.
  # It does not seem to change the 'best arm' which is the only
  # choice that matters.
  # Is this same as Gittins index? Can we reduce the computation?
  
  par(mfrow=c(2,3))
  maxt <- 8
  results <- c(1,0,0,0,0,0,0,0,0,0,0,0)
  mp <- data.frame(a=c(1,4,7),b=c(1,3,5))
  
  getVals <- function (mp,maxt) {
    vals <- matrix(data=rep(0,maxt*3),nrow=maxt)
    for(r in 0:(maxt-1)){
      qvals <- barl.q.all(mp,t=r)
      vals[r+1,] <- qvals / sum(qvals) # sum(qvals)
    }
    vals
  }
  
  for(t in 1:6){
    
    vals <- getVals(mp,maxt=(maxt-t+1))
    plot(x=t:maxt,vals[,1],type='l',ylim=c(min(vals),max(vals)),xlim=c(1,maxt),
         main=paste("Index value of each arm \n",toString(mp)))
    lines(x=t:maxt,vals[,2],col="red")
    lines(x=t:maxt,vals[,3],col="blue")
    abline(v=t)
    
    end <- vals[maxt-t+1,]
    lever <- which(end==max(end))[1]
    s <- results[t]
    mp[lever,] <- mp[lever,] + c(s,1-s)
  }
  
}

