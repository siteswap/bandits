
# Produce various correlation measures
.wcor <- function(x,y,w){
  x <- x - mean(x)
  y <- y - mean(y)
  corr <- (x*y)/(sd(x)*sd(y))
  sum(corr*w)/sum(w)
}

adtk.corr <- function(d){

	results <- ddply(d[d$c>0,],~g,summarize,
      		Cor=cor(a/c,c/n),
		IWC=.wcor(a/c,c/n,n),
		logIWCc=.wcor(a/c,log(c/n),n),
		RankCor=cor(a/c,c/n,method="spearman"),
		items=length(n),imps=sum(n))
	
  coeffs <- c()
	signifs <- c()
  for(g in unique(d$g)){
    fit <- glm(cbind(a,c-a)~log(p),data=d[(d$c>0) & (d$g==g) ,],family=binomial)
    sc <- summary(fit)$coefficients
    coeffs <- append(coeffs, sc["log(p)","Estimate"])
    signifs <- append(signifs, sc["log(p)","Pr(>|z|)"])
  }
  
	results$glmcoef <- coeffs
	results$glmsignif <- signifs
	results
}

############
# Test 1 - Frequentist goodness-of-fit test for simple binom model
############

adtk.test1 <- function(k,n,grpt=10){ # grpt - threshold for grouping small samples
  
  grp <- ifelse(n > grpt, 1:length(n), rep(-1, length(n)) )
  df <- data.frame(grp,k,n)
  D <- ddply(df, ~grp, summarize, k=sum(k), n=sum(n) )
  p <- sum(k)/sum(n) # ML estimation of rate
  r  <- (D$k - D$n*p)/(D$k*p*(1-p))
  1 - pchisq(q=sum(r^2),df=length(r))
}

############
# Test m1 vs m2 (not documented)
############

.lrTest <- function (k,n) {
  p0 <- sum(k)/sum(n)
  L0 <- sum( dbinom(k,size=n,prob=p0,log=TRUE) ) #  -629.711
  p1 <- k/n
  L1 <- sum( dbinom(k,size=n,prob=p1,log=TRUE) ) # -647.4449
  deviance <- -2*L0 + 2*L1
  1 - pchisq(q=deviance,df=length(n)-1)
}

############
# Test 2 - LR test of binomial vs betabinomial
############

library(VGAM)
library(inline) 
library(Rcpp)
library(rstan) 

adtk.test2 <- function(k,n,method="vglm"){
		
	if(method == "vglm"){
		return(.test2_vglm(k,n))
	}else if(method == "stan"){
		return(.test2_stan(k,n))
	}else{}
  }

.test2_vglm <- function(k,n){
  
  df <- data.frame(k,n)
  p0 <- sum(k)/sum(n)
  L0 <- sum( dbinom(k,size=n,prob=p0,log=TRUE) )
  
  fit <- vglm(cbind(k,n-k) ~ 1, betabinomial, data = df, trace = TRUE)
  #logit(coef(fit,matrix=TRUE)[1],inverse=TRUE) # Correctly finds true p value.
  L1 <- logLik(fit)
  deviance <- -2*L0 + 2*L1
  1 - pchisq(q=deviance,df=2-1)
  }

# TODO - treat with lazy singleton evaluation.
.bb_code='data {
    int<lower=1> N;          // Datapoints
    int<lower=0> a[N];       // aqs
    int<lower=1> c[N];       // views
  }
  parameters {
    real<lower=0> alpha;
    real<lower=0> beta;  
  }
  model {
    // Need to keep the explicit addition to get the real value of lp__
    increment_log_prob(beta_binomial_log(a,c,alpha,beta));
  }'

# We see the optimisation takes alpha as high as possible to minimize variance
# then beta is chosen to get the correct mean -> as values go larger,
# it converges in law/distribution to binomial.
.test2_stan <- function(k,n){

  p0 <- sum(k)/sum(n)
  L0 <- sum( dbinom(k,size=n,prob=p0,log=TRUE) )
  
  bb_dat <- list(N=length(k), a=k, c=n)
 
  bb_model <- stan_model(model_code=.bb_code)
  opt <- optimizing(object=.bb_model, data=bb_dat, algorithm="BFGS")
  L1 <- opt$value

  deviance <- -2*L0 + 2*L1

  print(opt)
  1 - pchisq(q=deviance,df=2-1)
}

##############
# test 3 - LR test 1 binom cluster vs 2 binom clusters.
##############
adtk.test3 <- function(k,n){

  rate <- k/n
  df <- data.frame(k,n)
  
  # Assume median values for cutoff
  t <- median(rate[rate>0])
  df$cluster <- (rate>t)*1
  
  # Single cluster
  m0 <- glm(cbind(k,n-k)~1,data=df,family=binomial)
  # 2 clusters
  m1 <- glm(cbind(k,n-k)~cluster,data=df,family=binomial)
  
  deviance <- -2*logLik(m0)[1] + 2*logLik(m1)[1]
  1 - pchisq(q=deviance,df=4-2)
}

##############
# Posterior test 1
##############

# Get posterior/s, perform CV.

library(plyr)
library(ggplot2)

# # Note, campaign_id is always 3597062
# d <- read.csv(file="~/site_domain_113803.csv")
# d <- d[ d$line_item_id==1023111, c("site_domain","imps","clicks","post_click_convs") ]
#       
# # Use only levers with atleast 10 appearances
# freq <- ddply(d,"site_domain",summarise,l=length(imps) )
# freqDomains <- freq[freq$l >= 10,1]
# d <- d[d$site_domain %in% freqDomains,]
# 
# # Randomize the order of the rows, and order by site_domain
# set.seed(seed=2105) # Set seed so sample is consistent
# d <- d[sample(1:dim(d)[1]),]
# 
# # Take half the data for training, half for testing
# halfOne <- function(l){ sum(l[1:ceiling(length(l)/2)]) }
# halfTwo <- function(l){ sum(l[(ceiling(length(l)/2)+1):length(l)]) }
# 
# train <- ddply(d,"site_domain",summarise,v=halfOne(imps),c=halfOne(clicks),a=halfOne(post_click_convs) )
# test <- ddply(d,"site_domain",summarise,v=halfTwo(imps),c=halfTwo(clicks),a=halfTwo(post_click_convs) )
# dirty <- function(t){ ((t$c < t$a) + (t$v < t$c)) != 0 }
# dirtyData <- union(which(dirty(train)),which(dirty(test)))
# train <- train[-dirtyData,]
# test <- test[-dirtyData,]

# # Test away !
# gmp.test(gmp(train),test)
# # 5.50897e-05
# 
# gmie.test(gmie(train),test)
# # 2.59007e-05
# 
# gm1.test(gm1(train),test)
# # 2.839369e-05 # warnings
# 
# gm2.test(gm2(train),test)
# # 2.790645e-05 # warnings
# 
# # Model 1 with mixture.
# gm1c.test(gm1c(train),test)
# 
# # Model 3 - something wrong.
# gm2c.test(gm2c(train),test)
# 


