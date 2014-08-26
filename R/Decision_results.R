
library(reshape)
library(ggplot2)
library(gridExtra)

compareStrats <- function(arms,initVals,campaignLen,trials,strats,strat_names,prior){
  
  RESAMPLE_OUTER <- TRUE
  results <- array(0,c(length(strats),campaignLen,trials))
  
  for(t in 1:trials){  
    
    if(RESAMPLE_OUTER){ # Resample
      initVals <- data.frame(q=rbeta(arms,shape1=prior$aq,shape2=prior$bq),
                             p=rbeta(arms,shape1=prior$ap,shape2=prior$bp))
    }
    
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
  meanSE <- boostrap.se(results,campaignLen,trials,strats)
  rownames(meanSE) <- strat_names
  
  
  df <- melt(means,id=c())
  df$se <- melt(sdevs,id=c())$value
  
  if(RESAMPLE_OUTER){ 
    df$se <- melt(meanSE,id=c())$value
  }
  colnames(df) <- c("strategy","round","regret","se")
  return(df)
}

# Calc standard error of mean
boostrap.se <- function(r,campaignLen,trials,strats){
  
  means <- array(0,c(length(strats),campaignLen,200))
  # size <- floor(trials/3) # TODO - ideal size for resample?
  
  for(i in 1:200){ # TODO - create config object to pass around.
    resample <- r[,, sample(1:trials,size=trials,replace=T) ]
    means[,,i] <- apply(X=resample,MARGIN=c(1,2),FUN=mean)
  }
  
  apply(X=means,MARGIN=c(1,2),FUN=sd)
}


#############
# Validation 1 - set q==1
#############

# With pin-point distribution on acquisition rate, all algos are doing the same thing - learning click rate.

arms <- 2
prior <- data.frame(aq=1000,bq=1,ap=1,bp=1) 
initVals <- data.frame(q=rbeta(arms,shape1=prior$aq,shape2=prior$bq),p=rbeta(arms,shape1=prior$ap,shape2=prior$bp))
campaignLen <- 200
trials <- 100
strats <- c(adtk.ts_clicks,adtk.ts_acqs,adtk.ts_both) # Careful, takes the values as the values at time of assignment
strat_names <- c("ts_clicks","ts_acqs","ts_both")


df <- compareStrats(arms,initVals,campaignLen,trials,strats,strat_names,prior)

ggplot(data=df, aes(x=round,y=regret,colour=strategy)) + 
  geom_line() +
  geom_errorbar(aes(ymin=regret-se, ymax=regret+se), width=.1, alpha=.5) +
  ggtitle("Validation \n (p==1)")


#############
# Validation 2 - set p==1
#############

# Clicks should never learn
# acqs and both should be equivalent

prior <- data.frame(aq=1,bq=1,ap=1000,bp=1) 
initVals <- data.frame(q=rbeta(arms,shape1=prior$aq,shape2=prior$bq),p=rbeta(arms,shape1=prior$ap,shape2=prior$bp))

df <- compareStrats(arms,initVals,campaignLen,trials,strats,strat_names,prior)

ggplot(data=df, aes(x=round,y=regret,colour=strategy)) + 
  geom_line() +
  geom_errorbar(aes(ymin=regret-se, ymax=regret+se), width=.1, alpha=.5) +
  ggtitle("Validation \n (q==1)")


###############
# Experiment 1 - compare Thompson Sampling of clicks vs acquisitions vs both
###############

# Remember - clicks will learn the best p fast.
# When this is not same as best r, it will never converge, 
# but when it does, it will look really good. 
# Should test over many values of init parameters.

arms <- 10
prior <- data.frame(aq=5,bq=20,ap=5,bp=5) 
initVals <- data.frame(q=rbeta(arms,shape1=prior$aq,shape2=prior$bq),p=rbeta(arms,shape1=prior$ap,shape2=prior$bp))
campaignLen <- 500
trials <- 100
strats <- c(adtk.ts_clicks,adtk.ts_acqs,adtk.ts_both) # Careful, takes the values as the values at time of assignment
strat_names <- c("ts_clicks","ts_acqs","ts_both")


df1 <- compareStrats(arms,initVals,campaignLen,trials,strats,strat_names,prior)
# df1_r <- df1[df1$round %in% seq(1,campaignLen,by=25),]

p1 <- ggplot(data=df1, aes(x=round,y=regret,colour=strategy)) + 
  geom_line() +
  geom_errorbar(aes(ymin=regret-se, ymax=regret+se), width=.2, alpha=.5) +
  ggtitle( paste(
                 "\n Arms: ",toString(arms),
                 "\n prior(q):",toString(prior[,c(1,2)]),
                 "\n prior(p):",toString(prior[,c(3,4)])
                  ))


prior <- data.frame(aq=1,bq=1,ap=1,bp=1) 
initVals <- data.frame(q=rbeta(arms,shape1=prior$aq,shape2=prior$bq),p=rbeta(arms,shape1=prior$ap,shape2=prior$bp))
df2 <- compareStrats(arms,initVals,campaignLen,trials,strats,strat_names,prior)

p2 <- ggplot(data=df2, aes(x=round,y=regret,colour=strategy)) + 
  geom_line() +
  geom_errorbar(aes(ymin=regret-se, ymax=regret+se), width=.1, alpha=.5) +
  ggtitle( paste(
                 "\n Arms: ",toString(arms),
                 "\n prior(q):",toString(prior[,c(1,2)]),
                 "\n prior(p):",toString(prior[,c(3,4)])
  ))

prior <- data.frame(aq=1,bq=1,ap=2,bp=10) 
initVals <- data.frame(q=rbeta(arms,shape1=prior$aq,shape2=prior$bq),p=rbeta(arms,shape1=prior$ap,shape2=prior$bp))
df3 <- compareStrats(arms,initVals,campaignLen,trials,strats,strat_names,prior)

p3 <- ggplot(data=df3, aes(x=round,y=regret,colour=strategy)) + 
  geom_line() +
  geom_errorbar(aes(ymin=regret-se, ymax=regret+se), width=.1, alpha=.5) +
  ggtitle( paste(
                 "\n Arms: ",toString(arms),
                 "\n prior(q):",toString(prior[,c(1,2)]),
                 "\n prior(p):",toString(prior[,c(3,4)])
  ))


prior <- data.frame(aq=2,bq=10,ap=1,bp=1) 
initVals <- data.frame(q=rbeta(arms,shape1=prior$aq,shape2=prior$bq),p=rbeta(arms,shape1=prior$ap,shape2=prior$bp))
df4 <- compareStrats(arms,initVals,campaignLen,trials,strats,strat_names,prior)

p4 <- ggplot(data=df4, aes(x=round,y=regret,colour=strategy)) + 
  geom_line() +
  geom_errorbar(aes(ymin=regret-se, ymax=regret+se), width=.1, alpha=.5) +
  ggtitle( paste(
                 "\n Arms: ",toString(arms),
                 "\n prior(q):",toString(prior[,c(1,2)]),
                 "\n prior(p):",toString(prior[,c(3,4)])
  ))

# Something where acqs is really bad
prior <- data.frame(aq=1,bq=20,ap=1,bp=1) 
initVals <- data.frame(q=rbeta(arms,shape1=prior$aq,shape2=prior$bq),p=rbeta(arms,shape1=prior$ap,shape2=prior$bp))
df5 <- compareStrats(arms,initVals,campaignLen,trials,strats,strat_names,prior)

p5 <- ggplot(data=df5, aes(x=round,y=regret,colour=strategy)) + 
  geom_line() +
  geom_errorbar(aes(ymin=regret-se, ymax=regret+se), width=.1, alpha=.5) +
  ggtitle( paste(
    "\n Arms: ",toString(arms),
    "\n prior(q):",toString(prior[,c(1,2)]),
    "\n prior(p):",toString(prior[,c(3,4)])
  ))

# Same again
arms <- 30
prior <- data.frame(aq=1,bq=20,ap=1,bp=1) 
initVals <- data.frame(q=rbeta(arms,shape1=prior$aq,shape2=prior$bq),p=rbeta(arms,shape1=prior$ap,shape2=prior$bp))
df6 <- compareStrats(arms,initVals,campaignLen,trials,strats,strat_names,prior)

p6 <- ggplot(data=df6, aes(x=round,y=regret,colour=strategy)) + 
  geom_line() +
  geom_errorbar(aes(ymin=regret-se, ymax=regret+se), width=.1, alpha=.5) +
  ggtitle( paste(
    "\n Arms: ",toString(arms),
    "\n prior(q):",toString(prior[,c(1,2)]),
    "\n prior(p):",toString(prior[,c(3,4)])
  ))


grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2, 
             main = "Thompson Sampling \n clicks vs acquisitions vs both ")


# ggplot(data=df, aes(x=X2,y=value,colour=X1)) + 
#   geom_line() +
#   geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1, alpha=.5) +
#   ggtitle( paste("Thompson Sampling \n clicks vs acquisitions vs both ",
#                  "\n Arms: ",toString(arms),
#                  "\n r:",toString(round(initVals$q*initVals$p,2)),
#                  "\n p:",toString(round(initVals$p,2)),
#                  "\n q:",toString(round(initVals$q,2))
#                  ))




###############
# Experiment 2 - Thompson Sampling vs Bayes Adaptive (acquisitions only)
###############

arms <- 3 # Stick to 2 arms for these short campaigns
prior <- data.frame(aq=1,bq=1,ap=1000,bp=1) 
initVals <- data.frame(q=rbeta(arms,shape1=prior$aq,shape2=prior$bq),p=1) 
campaignLen <- 10
trials <- 10
strats <- c(adtk.ts_acqs,adtk.barl)
strat_names <- c("ts_acqs","ba_acqs")


df <- compareStrats(arms,initVals,campaignLen,trials,strats,strat_names,prior)

ggplot(data=df, aes(x=round,y=regret,colour=strategy)) + 
  geom_line() +
  geom_errorbar(aes(ymin=regret-se, ymax=regret+se), width=.1, alpha=.5) +
  ggtitle( paste("Thompson Sampling vs Bayes Adaptive \n (acquisitions only) ",
               "\n Arms: ",toString(arms),
               "\n r",toString(round(initVals$q*initVals$p,2)),
               "\n p:",toString(round(initVals$p,2)),
               "\n q:",toString(round(initVals$q,2))
                ))



###############
# Experiment 3 - TS vs BS by acqs and both
###############

# REMEMBER - clicks will learn the best p fast.
# When this is not same as best r, it will never converge, 
# but when it does, it will look really good. 
# Should test over many values of init parameters.


arms <- 2 # Stick to 2 arms for these short campaigns
prior <- data.frame(aq=2,bq=6,ap=1,bp=1) 
initVals <- data.frame(q=rbeta(arms,shape1=prior$aq,shape2=prior$bq),p=rbeta(arms,shape1=prior$ap,shape2=prior$bp))
campaignLen <- 15
trials <- 1000
strats <- c(adtk.ts_both,adtk.barl,adtk.barl_both,adtk.ts_clicks)
strat_names <- c("ts_both","ba_acqs","ba_both","ts_clicks")

df <- compareStrats(arms,initVals,campaignLen,trials,strats,strat_names,prior)

ggplot(data=df, aes(x=round,y=regret,colour=strategy)) + 
  geom_line() +
  geom_errorbar(aes(ymin=regret-se, ymax=regret+se), width=.1, alpha=.5) +
  ggtitle( paste("Thompson Sampling vs Bayes Adaptive RL",
    "\n Arms: ",toString(arms),
    "\n prior(q):",toString(prior[,c(1,2)]),
    "\n prior(p):",toString(prior[,c(3,4)])
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
  
  par(mfrow=c(3,2))
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
         ylab="Q value",xlab="Time Horizon",
         main=paste("Q value of arms",
                    "\n Alpha:",toString(mp$a),
                    "\n Beta:",toString(mp$b)
                    ))
    lines(x=t:maxt,vals[,2],col="red")
    lines(x=t:maxt,vals[,3],col="blue")
    abline(v=t)
    
    end <- vals[maxt-t+1,]
    lever <- which(end==max(end))[1]
    s <- results[t]
    mp[lever,] <- mp[lever,] + c(s,1-s)
  }
  
}

#########################################
### Benchmarking Gittins Index policy ###
#########################################


arms <- 10 # Stick to 2 arms for these short campaigns
prior <- data.frame(aq=1000,bq=1,ap=1,bp=1) 
initVals <- data.frame(q=rbeta(arms,shape1=prior$aq,shape2=prior$bq),p=rbeta(arms,shape1=prior$ap,shape2=prior$bp))
campaignLen <- 40
trials <- 100
strats <- c(adtk.gi_clicks,adtk.ts_clicks)
strat_names <- c("gi_clicks","ts_clicks")

df <- compareStrats(arms,initVals,campaignLen,trials,strats,strat_names,prior)

ggplot(data=df, aes(x=round,y=regret,colour=strategy)) + 
  geom_line() +
  geom_errorbar(aes(ymin=regret-se, ymax=regret+se), width=.1, alpha=.5) +
  ggtitle( paste("Gittins Index vs Thompson Sampling - Clicks only",
                 "\n Arms: ",toString(arms),
                 "\n prior(q):",toString(prior[,c(1,2)]),
                 "\n prior(p):",toString(prior[,c(3,4)])
  ))

