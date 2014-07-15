
library(ggplot2)
library(plyr)
library(VGAM)
library(inline) 
library(Rcpp)
library(rstan) 


# We give 3 example datasets to illustrate these techniques.

###################
## Independent sim
###################

# How to get the distribution for n?
# Since this is not part of our model, lets just sample it from the data.
setwd("~")
site_domain <- read.csv(file="site_domain_113803.csv")

rindie <- function(num,p,q){
  
  n <- sample(site_domain$imps,size=num,replace=TRUE)
  c <- rbinom(num,size=n,prob=p)
  a <- rbinom(num,size=c,prob=q)
  data.frame(a,c,n)
}

d <- rindie(num=100000,p=0.001,q=0.6) # d <- rindie(num=100000,p=0.0001,q=0.6)
d$p <- d$c/d$n
d$q <- d$a/d$c
d$g <- 1 # Create group to summarize
dindie <- d

# 100,000 is similar to the data
ggplot(data=d) + 
  aes(x = log(a/c), y = log(c/n)) +
  geom_point(aes(size = log(n),alpha = 0.3))

##################################
## Testing - expect no correlation



ddply(d[d$c>0,],~g,summarize,
      Cor=cor(a/c,c/n),ImpWeigthedCor=wcor(a/c,c/n,n),liwc=wcor(a/c,log(c/n),n),RankCor=cor(a/c,c/n,method="spearman"),items=length(n),imps=sum(n))
# g         Cor ImpWeigthedCor        liwc     RankCor items     imps
# 1 1 -0.08153287    -0.01410091 -0.05835615 -0.01166252   369 10869547

ddply(d[d$c>1,],~g,summarize, # No correlation in this group.
      Cor=cor(a/c,c/n),ImpWeigthedCor=wcor(a/c,c/n,n),liwc=wcor(a/c,log(c/n),n),RankCor=cor(a/c,c/n,method="spearman"),items=length(n),imps=sum(n))
# g        Cor ImpWeigthedCor        liwc    RankCor items     imps
# 1 1 -0.2923502    -0.05013228 -0.08579775 -0.2093781   134 10050497

ddply(d[d$a>0,],~g,summarize, # Created correlation!
      Cor=cor(a/c,c/n),ImpWeigthedCor=wcor(a/c,c/n,n),liwc=wcor(a/c,log(c/n),n),RankCor=cor(a/c,c/n,method="spearman"),items=length(n),imps=sum(n))
# g       Cor ImpWeigthedCor      liwc   RankCor items     imps
# 1 1 0.1479859      0.1549478 0.5607463 0.4642423   260 10485552


summary(glm(cbind(a,c-a)~log(p),data=d[d$c>0,],family=binomial)) # No evidence of correlation.
summary(glm(q~log(p),data=d[d$a>0,])) # Strong evidence. We must be creating some relationship with this selection.

#############################
## Dispersed data
#############################


rbb <- function(num){
  
  n <- sample(site_domain$imps,size=num,replace=TRUE)
  c <- rbetabinom.ab(n=num,size=n,shape1=1,shape2=10000)
  a <- rbetabinom.ab(num,size=c,shape1=6,shape2=4)
  data.frame(a,c,n)
}

d <- rbb(num=100000)
d$p <- d$c/d$n
d$q <- d$a/d$c
d$g <- 1 # Create group to summarize
dbb <- d

ggplot(data=d) + 
  aes(x = log(a/c), y = log(c/n)) +
  geom_point(aes(size = log(n),alpha = 0.3))


#############################
## Cluster based correlation
#############################

rclust<- function(num,p,q,theta){
  
  n <- sample(site_domain$imps,size=num,replace=TRUE)
  clust <- rbinom(n=num,size=1,prob=theta)
  c <- rbinom(num,size=n,prob=p[clust+1]) # Vectorization might perform better to group each cluster
  a <- rbinom(num,size=c,prob=q[clust+1])
  data.frame(a,c,n,clust)
}

d <- rclust(num=100000,p=c(0.0005,0.00005),q=c(0.8,0.2),theta=0.9)
d$p <- d$c/d$n
d$q <- d$a/d$c
d$g <- 1 # Create group to summarize
dclust <- d

ggplot(data=d) + 
  aes(x = log(a/c), y = log(c/n)) +
  geom_point(aes(size = log(n),alpha = 0.3, color=as.factor(c("good","bad"))[clust+1]))

#############################
## Testing - expect no correlation



ddply(d[d$c>0,],~g,summarize, # Weak evidence for correlation
      Cor=cor(a/c,c/n),ImpWeigthedCor=wcor(a/c,c/n,n),liwc=wcor(a/c,log(c/n),n),RankCor=cor(a/c,c/n,method="spearman"),items=length(n),imps=sum(n))
# g       Cor ImpWeigthedCor   RankCor items     imps
# 1 1 0.0678018     0.02881429 0.1999053   328 11878585

ddply(d[d$c>1,],~g,summarize,
      Cor=cor(a/c,c/n),ImpWeigthedCor=wcor(a/c,c/n,n),liwc=wcor(a/c,log(c/n),n),RankCor=cor(a/c,c/n,method="spearman"),items=length(n),imps=sum(n))
# g       Cor ImpWeigthedCor  RankCor items     imps
# 1 1 0.4297301       0.187535 0.589616   135 12965971

ddply(d[d$a>0,],~g,summarize,
      Cor=cor(a/c,c/n),ImpWeigthedCor=wcor(a/c,c/n,n),liwc=wcor(a/c,log(c/n),n),RankCor=cor(a/c,c/n,method="spearman"),items=length(n),imps=sum(n))
# g       Cor ImpWeigthedCor   RankCor items    imps
# 1 1 0.1775996      0.2135723 0.6933963   165 9564968

summary(glm(cbind(a,c-a)~p,data=d[d$c>0,],family=binomial)) # no evidence - need to use log values
summary(glm(cbind(a,c-a)~log(p),data=d[d$c>0,],family=binomial)) # Very strong evidence
summary(glm(q~log(p),data=d[d$a>0,])) # Even stronger for a > 0

#########################################
## Linear correlation is more complicated
#########################################

# TODO

###########
## Summary
###########

# 0) Need to use log values of p.
# 1) Clustering based on a_i creates false correlations.
# 2) Clustering based on c does not create false correlations.
# 3) Using c>0, correlations do not show up. But become stronger
#    for c>1, c>2. This is because it reduces the noise of all the small-n zero and one values.
# 4) Filter based on realized CTR is no use as it does not help with (3). 
#     Needs a correlation measure which properly accounts for likelihood weight of response - closest thing is glm.
# 5) Binomial GLM does the best job of (4), though does not consider error in dependent variable.


###################################################
###################################################
## Hypothesis testing 
###################################################
###################################################


# Approaches to Bayesian model diagnostics:
#  Posterior predictive checks
#  Info criteria
#  Bayes factors (marginal likelihood)
#  Sensitivity analysis (change priors/model see difference)

# Read part II BDA3.

# CHECKING MODELS - many uses of posterior predictive replicated data (can use stan for this):
#   - how well does posterior predictive predict (you need to create some 
#         loss function / discrepency measure appropriate to your use). Generate a
#         distribution for the test statistic under the 'model is correct' hypothesis.
#         Then calculate p-value for the observed data. See pg 146. You can do all
#         sorts of stuff like this. Let's you turn this in to hypothesis test.
#   - sample from posterior predictive - visually compare to real data.

# MODEL COMPARISON - which model is better, but not is this model good.
#   - Information criteria - estimates/approximations to cross/externally validated data.
#   - Bayes factors (use model evidence aka marginal likelihood)

# ANOVA - a series of submodels, where you compare the ratio of errors.
#       - assumes normally distributed errors so you can do an F-test.
#       - what if errors not N, what if only 1 sample per group.
# Since ANOVA is comparing errors in NLM, we can use GLM to generate errors.
# we would need to define a dummy variable for every lever 
# (since we are giving each one it's own group). Alternatively, we could model the
# true values as random variables using faraway glme? (with beta?)
# We would compare cbind(a,c)~1 and cbind(a,c-a) ~ factor(1:length(a))
# ANOVA on glm gives you the deviances, but no f-test (since f-test is only for Normal residuals).


# Are p_i all the same, are q_i all the same.
#   - LR test - but how to stop it preferring the more complex? And how to not overfit?
#   - Marginal likelihood for model selection between the 2 models.
#   - AIC / BIC
#   - Bayesian hypothesis testing
# Are true p_i and q_i Beta distributed
#   - Look at tests for beta-binomial, UMP goodness of fit?
#   - LR test
# Are p_i and q_i correlated in some way.
#   - LR test
#   - AIC
#   - 'Bacon With Your Eggs? Applications of a New Bivariate Beta-Binomial Distribution'
#   - 'The Use of a Correlated Binomial Model for the Analysis of Certain Toxicological Experiments'
# Are p_i and q_i clustered.
#   - test/s discussed with AC
# Lee - 'Properties and Applications of the Sarmanov Family of Bivariate Distributions'
# g(p1, p2) = f(p1|??1, ??1)f(p2|??2, ??2) ?[1 +??(p1 ????1)(p2 ????2)].

# Now put it all together in a nice library!

# GEE ?
# VGLM ?
# Log Linear ?

#####################
# Test - are binomials drawn from same distribution?
# Research 'Testing Equality of Binomial Proportions'

# 0) t-statistic for testing true means when variance is unknown.

# If all our samples had the same n, we would aggregate the successes and see a binomial distribution,
# which in large n is approximated by Normal. We would then perform chi-squared to test goodness of fit.
# The residuals are then chisqr. This is same procedure as Pearsons chisqr.
# BUT - our n are not the same, and many are small!

# -> Idea, group the small n together. We can use chisq
# -> If we choose a p, we can use pearson chisq.

# Pearson chi squared test for goodness-of-fit of categorical data
# Contents of each bucket minus the expected value is sum of iid r.v.
# By CLT, this is N(0,Var). This is the same proof for binomial approx to Normal.
# We then add sum the errors for each bucket 
# to get a chisq with #bucket-1 degrees of freedom.

# 1) Assume normal approximation to binomial ( fishy for small n ):
# The chisquared statistic is the function meeting the Monotone Likelihood Ratio
# requirement. Therefore test is UMP.
# For normally distributed variables, under the null hypothesis,
# the SSE follows a chi-squared distribution.

# 1.1) If we group everything, this *is* Bin by definition!


# This is not good for q, as it's very close to zero.

# Pearson chisq (relies on normal approximation? should be able to handle different n)
# Yates chisq (corrects pearson for approximation error)

# Fisher's exact test compares 2 (not UMP, also see Barnard exact test).
# 'An exact method of testing equality of several binomial proportions to a specified standard'
# 'Testing the equality of several binomial proportions to a prespecified standard' - EVEN BETTER!

# Attempt at an 'exact test':
# It would be nice to get the UMP benefits of (1) without the approximation
# assumption. But can't get max invariant statitic to reduce to point-vs-interval
# and then MLR to get UMP. Best can do is make up a stat.
# 2) H0: All samples are drawn from common p - is this an MLR family?
# (p is a nuisance parameter as we don't care what value it actually takes). 
# It's not clear what statistic we can use here to get UMP test.
# If we knew the p to test, we could define some test statistic.
# We might then compute the distribution of this statistic trhough simulation 
# and decide whether the observed statistic was sufficiently extreme to reject.

# Formed as a model comparison:

# 3) Form LR test
# If we can specify alternative model we can use LR test. As ML solution for p is used. 
# ML solution for p is a/c.
# Compare common p vs independent p's.
# Calculate the LR statistic (aka Deviance). 
# This is approx chisq distributied with df1 - df2 degrees of freedom.
# We would be forced to use AIC instead if the models were not nested, but no need here.
# Validate your results against glm (though glm should be slightly less accurate because 
# it uses computational solution, where you use analytic solution).


# 4) Bayes factor 
# We may question whether ML point estimates are representative and robust.
# If instead of using ML values for the test, we defined some prior over p values,
# we could use then marginalize out the parameter and the LR test becomes a Bayes factor test.
# Stan can do this for you.


##########################
# Test single p and q
##########################

# Also,
# 'Goodness-of-Fit Issues in Toxicological Experiments Involving Litters of Varying Size'

# 1) Group the low scoring data and test - passes very strongly.
# Benefit is that we don't have to specify an alternative.

pearson <- function(k,n,grpt=10){ # grpt - threshold for grouping small samples
  
  grp <- ifelse(n > grpt, 1:length(n), rep(-1, length(n)) )
  df <- data.frame(grp,k,n)
  D <- ddply(df, ~grp, summarize, k=sum(k), n=sum(n) )
  p <- sum(k)/sum(n) # ML estimation of rate
  r  <- (D$k - D$n*p)/(D$k*p*(1-p))
  1 - pchisq(q=sum(r^2),df=length(r))
}

pearson(dindie$a,dindie$c) # 1
pearson(dclust$a,dclust$c) # 0
# May be further improved by Yates correction.

# Fails for p even if we aggregate a lot.
# The true q is too small.
pearson(dindie$c,dindie$n,100000) 


# 3) LR test if we mind specifying a (ML) point hypothesis alternative.
# By Neyman-Pearson, is Most powerful for 2 point hypotheses.

lrTest <- function (k,n) {
  p0 <- sum(k)/sum(n)
  L0 <- sum( dbinom(k,size=n,prob=p0,log=TRUE) ) #  -629.711
  p1 <- k/n
  L1 <- sum( dbinom(k,size=n,prob=p1,log=TRUE) ) # -647.4449
  deviance <- -2*L0 + 2*L1
  1 - pchisq(q=deviance,df=length(n)-1)
}

# Too easily accepts the null as many points increase df
lrTest(dindie[dindie$c > 0,"a"],dindie[dindie$c > 0,"c"]) # 1.5e-13

d <- dclust
d$grp <- ifelse(d$c > 10, 1:length(d$c), rep(-1, length(d$c)) )
Dclust <- ddply(d, ~grp, summarize, a=sum(a), c=sum(c) )

# FAILS to reject the cluster model - alternative model is too complex
# Try beta-binomial and cluster as alternative
lrTest(Dclust$a,Dclust$c)


# 4) Bayes - we leave for later, but has benefit of not assuming MLE.
# This might be relevant if likelihood is sharply peaked for one hypothesis, 
# but not the other.


####################
# Test beta-binomial 
####################

# We can not accept the total independence assumption in the previous LR test.
# But nor can we accept the single p either.

# Testing null hypothesis against all others is difficult. This paper has a great
# discussion about the various approaches chosen:
# 'Bootstrap goodness-of-fit test for the beta-binomial model'

# See fitting dirichlet-multinomial aka multinomial-polya.

# We just go for LR:
# ML estimate for DirMult not available in closed form - use numerical optimization (see T.Minka 2003)
# R VGAM package does glm with beta binomial emissions. 
# Also RStan - try both to make sure results agree.
# library(VGAM)


dindiec <- dindie[dindie$c > 0,]
dclustc <- dclust[dclust$c > 0,]

lrTestBB <- function(k,n){
  
  df <- data.frame(k,n)
  p0 <- sum(k)/sum(n)
  L0 <- sum( dbinom(k,size=n,prob=p0,log=TRUE) )
  
  fit <- vglm(cbind(k,n-k) ~ 1, betabinomial, data = df, trace = TRUE)
  #logit(coef(fit,matrix=TRUE)[1],inverse=TRUE) # Correctly finds true p value.
  L1 <- logLik(fit)
  deviance <- -2*L0 + 2*L1
  1 - pchisq(q=deviance,df=2-1)
  }

dbbc <- dbb[dbb$c>0,]

lrTestBB(dindiec$a,dindiec$c) # Correctly keeps null
lrTestBB(dclustc$a,dclustc$c) # Rejects null, but with warnings?
lrTestBB(dbbc$a,dbbc$c)       # Correctly rejects null

# Very slow convergence for these values :(
dindien <- dindie[dindie$n > 0,]
lrTestBB(dindien$c,dindien$n) # 0.3042346 - accepts null, would like strong given amount of data.

dclustn <- dclust[dclust$n > 0,]
lrTestBB(dclustn$c,dclustn$n) # 0.875667 - Accepts (just) the cluster model

# CAREFUL - The sampling statement drops all the terms in the log probability function that are
# constant, whereas the explicit call to normal_log adds all of the terms in the defini-
# tion of the log normal probability function, including all of the constant normalizing
# terms. Therefore, the explicit increment form can be used to recreate the exact log
# probability values for the model.
bb_code <-  'data {
    int<lower=1> N;          // Datapoints
    int<lower=0> a[N];       // aqs
    int<lower=1> c[N];       // views
  }
  parameters {
    real<lower=0> alpha;
    real<lower=0> beta;  
  }
  model {
    // a ~ beta_binomial(c,alpha,beta);
    // Need to keep the explicit addition so that we include all
    // the additive constants and get the real value for lp__
    increment_log_prob(beta_binomial_log(a,c,alpha,beta));
  }'

bb_model <- stan_model(model_code=bb_code)

# We see the optimisation takes alpha as high as possible to minimize variance
# then beta is chosen to get the correct mean -> as values go larger,
# it converges in law/distribution to binomial
lrTestBB_stan <- function(k,n){

  p0 <- sum(k)/sum(n)
  L0 <- sum( dbinom(k,size=n,prob=p0,log=TRUE) )
  
  bb_dat <- list(N=length(k), a=k, c=n)
  
  opt <- optimizing(object=bb_model, data=bb_dat, algorithm="BFGS")
  L1 <- opt$value

  deviance <- -2*L0 + 2*L1

  print(opt)
  1 - pchisq(q=deviance,df=2-1)
}

# It's failing to find the ML estimate!!!
lrTestBB_stan(dindiec$a,dindiec$c) # 1
lrTestBB_stan(dclustc$a,dclustc$c) # 0 - reject the null
lrTestBB_stan(dbbc$a,dbbc$c)       # 0 - reject the null


###################
# Test clustering
###################
# Simple - split in to groups, based on CTR - regress and compare residuals.
# Can do this simple split with both binomial and betabinomial emissions.
# LR test - 1 cluste vs 2 clusters.

# A sample point (and therefore an ML estimate) requires a sample from Z
# which is a categorical vector. 
# If we are only interested in model comparison, we could sum out the 
# p and q values (so both models would be beta-binomial) and compare.
# But this still leaves problem of z values. Can we, for model comparison
# somehow sum them out too?

# 1) The z parameter is, a-priori, distributed according to Categorical(theta),
# so, a priori, the likelihood of a point is SUM that*BB. But in a true ML or MAP
# setting, we would weight the more likely cluster at 1 (not at .6 say, according to
# the posterior estimate of the distribution of z). We could instead consider each
# z to be sampled from an independent prior, but 

# 2) There is a sort of 2 step algo in the stan cluster code - it samples from posterior
# of the emission densities, given a distribution for responsibility param. It then 
# takes sets the responsibility param to be the Expectation given emission values.
# It's similar to EM, except that the M step is an HMC step, not a maximization.

# Either way, given the large number of indices, we would need EM or Gibbs (VB if we
# require speed) to efficiently derive posterior distributions of p and q.

# It's unclear what each of these models does for us - we are neither taking a 
# point hypothesis for LR test, nor being fully Bayesian and summing out all 
# parameters to produce a Bayes factor - write this up!

# 'Naive Bayes classifier' is a term applied to any clustering model where we assume 
# independence of the features. The features are the emissions (3D Gaussian for images,
# multinomial for word frequencies which is the generalization of our work). The 
# independence means we don't need to worry about covariances and so the (approximation) 
# algorithms are quite quick. An SVM, for example, would be able to consider more complex
# relationships in it's classification.



LRCluster <- function(d){

  # Assume median values for cutoff
  t <- median(d$p[d$p>0])
  d$cluster <- (d$p>t)*1
  
  # Single cluster
  m0 <- glm(cbind(a,c-a)~1,data=d,family=binomial)
  # 2 clusters
  m1 <- glm(cbind(a,c-a)~cluster,data=d,family=binomial)
  
  logLik(m1)[1]/logLik(m0)[1]
}


ggplot(data=d) + 
  aes(x = log(a/c), y = log(c/n)) +
  geom_point(aes(size = log(n),alpha = 0.3, color=as.factor(c("good","bad"))[cluster+1]))


LRCluster(dreal[dreal$a<=dreal$c,])

LRCluster(dclust[ (dclust$c>0),])
# 0.7130438




