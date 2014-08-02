
library(reshape)
library(ggplot2)

arms <- 3
initVals <- data.frame(q=rbeta(arms,shape1=5,shape2=95),p=rbeta(arms,shape1=1,shape2=1))
campaignLen <- 100
trials <- 10

ts_clicks <- array(0, c(campaignLen, trials))
ts_acqs <- array(0, c(campaignLen, trials))
ts_both <- array(0, c(campaignLen, trials))

# TODO - parrallelize
for(v in 1:trials){  
  # 1. one only tries to maximise number of clicks, ignoring conversions
  ts_clicks[,v] <- adtk.mab(DFUN=adtk.ts_clicks,trueVals=initVals,campaignLen=campaignLen,arms=arms)$metrics$regret
  # 2. one only tries to maximise number of conversions, ignoring clicks
  ts_acqs[,v] <- adtk.mab(DFUN=adtk.ts_acqs,trueVals=initVals,campaignLen=campaignLen,arms=arms)$metrics$regret
  # 3. one tries to maximise number of conversions, considering clicks as well, assuming however, 
  #     that the two rates are independent (formulate the model carefully and extract a confidence interval for UCB). 
  ts_both[,v] <- adtk.mab(DFUN=adtk.ts_both,trueVals=initVals,campaignLen=campaignLen,arms=arms)$metrics$regret
}

r <- rbind(t(ts_clicks),t(ts_acqs),t(ts_both))
df <- as.data.frame(r)

means <- ddply(df, .(g=rep(as.factor(c("clicks","acqs","both")),each=10)), 
            function(x){apply(X=x,MARGIN=2,FUN=mean)})
m <- as.data.frame(t(means[,-1]))
colnames(m) <- c("clicks","acqs","both")
m$t <- 1:campaignLen



sd <- ddply(df, .(g=rep(as.factor(c("clicks","acqs","both")),each=10)), 
               function(x){apply(X=x,MARGIN=2,FUN=sd)})
s <- as.data.frame(t(sd[,-1]))
colnames(s) <- c("clicks","acqs","both")
s$t <- 1:campaignLen

all <- melt(m,id=c("t"))
all$se <- melt(s,id=c("t"))$value



ggplot(data=all, aes(x=t,y=value,colour=variable)) + 
  geom_line() +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1, alpha=.5)





###############
# Experiment 2
###############

L <- arms <- 3 # TODO - stupid hardcoded L value.
initVals <- data.frame(q=rbeta(arms,shape1=1,shape2=1),p=1) # q is anything from unif(0,1)
campaignLen <- 15
trials <- 10


ts_acqs <- array(0, c(campaignLen, trials))
barl <- array(0, c(campaignLen, trials))
winner <- array(0, c(2, trials))

for(v in 1:trials){  

  ts_acqs[,v] <- r1 <- adtk.mab(DFUN=adtk.ts_acqs,trueVals=initVals,campaignLen=campaignLen,arms=arms)$metrics$regret 
  barl[,v] <- r2 <- adtk.mab(DFUN=adtk.barl,trueVals=initVals,campaignLen=campaignLen,arms=arms)$metrics$regret
  winner[1,v] <- 1*(r1[campaignLen] < r2[campaignLen]) # Winner has least regret
  winner[2,v] <- 1*(r1[campaignLen] > r2[campaignLen])
  # TODO
  # Winner isn't really meaningful here, because we simulate 2 different games.
  # We would need to draw 10 samples from each arm first, then compare.
}



r <- rbind(t(ts_acqs),t(barl))
df <- as.data.frame(r)

means <- ddply(df, .(g=rep(as.factor(c("acqs","barl")),each=trials)), 
               function(x){apply(X=x,MARGIN=2,FUN=mean)})
m <- as.data.frame(t(means[,-1]))
colnames(m) <- c("acqs","barl")
m$t <- 1:campaignLen



sd <- ddply(df, .(g=rep(as.factor(c("acqs","barl")),each=trials)), 
            function(x){apply(X=x,MARGIN=2,FUN=sd)})
s <- as.data.frame(t(sd[,-1]))
colnames(s) <- c("acqs","barl")
s$t <- 1:campaignLen

all <- melt(m,id=c("t"))
all$se <- melt(s,id=c("t"))$value


tot <- apply(X=winner,MARGIN=1,FUN=sum)

ggplot(data=all, aes(x=t,y=value,colour=variable)) + 
  geom_line() +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1, alpha=.5) +
  ggtitle(paste("acqs",tot[1],"-",tot[2],"barl",": true rate",toString(round(initVals$q,2))))

# Look at the literature on this, there are many 
# approaches (Action elimination, approximations)
# to deal with the size of the knowledge state space.
# Difference between bandits and RL? (RL has states?)

###############
# Experiment 3
###############


arms <- 2
initVals <- data.frame(q=rbeta(arms,shape1=10,shape2=20),p=rbeta(arms,shape1=1,shape2=1))
campaignLen <- 15
trials <- 100

ts_both <- array(0, c(campaignLen, trials))
barl_both <- array(0, c(campaignLen, trials))

for(v in 1:trials){  
  ts_both[,v] <- adtk.mab(DFUN=adtk.ts_both,trueVals=initVals,campaignLen=campaignLen,arms=arms)$metrics$regret
  barl_both[,v] <- adtk.mab(DFUN=adtk.barl_both,trueVals=initVals,campaignLen=campaignLen,arms=arms)$metrics$regret
  # TODO - too noisy. Include the winner by comparing the same game
}

r <- rbind(t(ts_both),t(barl_both))
df <- as.data.frame(r)

means <- ddply(df, .(g=rep(as.factor(c("ts","barl")),each=trials)), 
               function(x){apply(X=x,MARGIN=2,FUN=mean)})
m <- as.data.frame(t(means[,-1]))
colnames(m) <- c("ts","barl")
m$t <- 1:campaignLen



sd <- ddply(df, .(g=rep(as.factor(c("ts","barl")),each=trials)), 
            function(x){apply(X=x,MARGIN=2,FUN=sd)})
s <- as.data.frame(t(sd[,-1]))
colnames(s) <- c("ts","barl")
s$t <- 1:campaignLen

all <- melt(m,id=c("t"))
all$se <- melt(s,id=c("t"))$value


# tot <- apply(X=winner,MARGIN=1,FUN=sum)

ggplot(data=all, aes(x=t,y=value,colour=variable)) + 
  geom_line() +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1, alpha=.5) 
# TODO - add winner count




