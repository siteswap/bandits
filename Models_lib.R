
library(inline) 
library(Rcpp)
library(rstan) 

###############
# Baseline non-hierarchical models:
# 1) Complete pooling
# 2) Individual estimation
###############

# Complete pooling
gmp <- function(d){
  predictedCVR <- sum(d$a) /sum(d$v)
  predictedCVR
}

gmp.test <- function(gmpfit,t){
  realizedCVR <- t$a / t$v
  risk <- abs(realizedCVR - gmpfit)
  # risk[is.nan(risk)] <- 0
  sum(risk * t$v) / sum(t$v)
}

# Individual estimation
gmie <- function(d){
  predictedCVR <- d$a / d$v
  predictedCVR
}

gmie.test <- function(gmiefit,t){
  realizedCVR <- t$a / t$v
  risk <- abs(realizedCVR - gmiefit)
  sum(risk * t$v) / sum(t$v)
}

##########
# Model 1
##########

# Expects v, c, a
gm1 <- function(d,prior=c(-1,-1)){

  # Take Empirical Bayes approach to priors.
  if( all( prior == c(-1,-1) ) ){  
    # Estimate data driven parameters:
    # Can not use method of moments as var < mean, use maximum likelihood
    prior <- optim(par=c(1,1), fn=gm1.logLike, vi=d$v, ai=d$a, method = "Nelder-Mead")$par
    }
  
  d$alpha <- prior[1]
  d$beta <- prior[2]
  
  d$alpha_p <- d$alpha + d$a
  d$beta_p  <- d$v - d$a + d$beta
  d
}

gm1.logLike <- function(x,vi,ai){
  
  a <- x[1]
  b <- x[2]
  # Beta-Binomial model - intermediate p parameter integrates to a Beta.
  ll <- sum( lchoose(vi,ai)
             + lbeta( ai+a, vi-ai+b )
             - lbeta(a,b)
  )
  -ll # optim minimizes
}

# Calculate volume weighted mean loss - vml1
# Expects test data (t) to have same ordering as training data.
gm1.test <- function(gm1fit,t){
  
  # Taking the Bayesian view, we calculate the expected loss.
  # This works out simply as the true value minus the expectation ( which is alpha/(alpha + beta) )
  realizedCVR <- t$a / t$v
  # Loss is difference in values, risk is expected loss
  risk <- abs(realizedCVR - gm1fit$alpha_p/(gm1fit$alpha_p + gm1fit$beta_p))
  
  # Volume weighted mean loss:
  risk[is.nan(risk)] <- 0
  sum(risk * t$v) / sum(t$v)
}


##########
# Model 2
##########

# The 2 unknown variables are independent, given clicks
gm2 <- function(d,prior_q=c(-1,-1),prior_p=c(-1,-1)){
  
  # Take Empirical Bayes approach to priors.
  if( all( prior_q == c(-1,-1) ) ){  
    prior_q <- optim(par=c(1,1), fn=gm1.logLike, vi=d$c, ai=d$a, method = "Nelder-Mead")$par
  }
  d$alpha_q <- prior_q[1]
  d$beta_q <- prior_q[2]
  
  if( all( prior_p == c(-1,-1) ) ){  
    prior_p <- optim(par=c(1,1), fn=gm1.logLike, vi=d$v, ai=d$c, method = "Nelder-Mead")$par
  }
  d$alpha_p <- prior_p[1]
  d$beta_p <- prior_p[2]
  
  d$alpha_p1 <- d$alpha_p + d$c
  d$beta_p1  <- d$v - d$c + d$beta_p
  d$alpha_q1 <- d$alpha_q + d$a
  d$beta_q1  <- d$c - d$a + d$beta_q
  d
}

# Calculate volume weighted mean loss
gm2.test <- function(gm2fit,t){
  
  # Taking the Bayesian view, we calculate the expected loss.
  # This works out simply as the true value minus the expectation ( which is alpha/(alpha + beta) )
  realizedCVR <- t$a / t$v
  realizedCVR[is.infinite(realizedCVR)] <- 0
  
  # Predicted CVR is q*p
  # Loss is difference in values, risk is expected loss
  # This works out as CVR - E[p]*E[q]
  risk <- abs(realizedCVR - (gm2fit$alpha_p1/(gm2fit$alpha_p1 + gm2fit$beta_p1))*(gm2fit$alpha_q1/(gm2fit$alpha_q1 + gm2fit$beta_q1)))
  
  # Volume weighted mean loss:
  risk[is.nan(risk)] <- 0
  sum(risk * t$v) / sum(t$v)
}



#################
# Model 1 Cluster
#################


# Note - with multiple chains, results suffer from index switching.
# Also need to be wary if index switching occurring within single chain.
gm1c <- function(d){
  
  gm1c_dat <- list(N = dim(d)[1], 
                   K = 2,
                   a = d$a,
                   v = d$v
                   )
  
  model1_clust_code <-  'data {
    int<lower=0> N;  // number of data points
    int<lower=1> K;  // number of clusters
    int<lower=0> a[N];       // aqs
    int<lower=1> v[N];       // views
  }
  transformed data {
    real<upper=0> neg_log_K;
    neg_log_K <- -log(K); // Equal prior weights
  }
  parameters {
    real<lower=0> alpha[K]; // cluster means
    real<lower=0> beta[K];  
  }
  transformed parameters {
    real<upper=0> soft_z[N,K]; // log unnormalized cluster assigns
    for (n in 1:N)
      for (k in 1:K)
        soft_z[n,k] <- neg_log_K + beta_binomial_log(a[n],v[n],alpha[k],beta[k]);
  }
  model {
    for (k in 1:K){
      alpha[k] ~ uniform(0,10000);  // prior??? 
      beta[k] ~ uniform(0,10000);  // prior???
    }
    for (n in 1:N)
      increment_log_prob(log_sum_exp(soft_z[n])); // likelihood
  }'
  
  # TODO - chains don't run in parrallel?!?!?
  # https://groups.google.com/forum/#!msg/stan-users/3goteHAsJGs/UmSDRo77Vn4J
  stan(model_code = model1_clust_code, data = gm1c_dat, iter = 1000, chains = 4)
  
}


#################
# Model 2 Cluster
#################


gm2c <- function(d){
  
  m_dat <- list(N = dim(d)[1], 
                   K = 2,
                   a = d$a,
                   c = d$c,
                   v = d$v
  )
  
  model_code <-  'data {
    int<lower=0> N;  // number of data points
    int<lower=1> K;  // number of clusters
    int<lower=0> a[N];       // acquisitions
    int<lower=0> c[N];       // clicks
    int<lower=1> v[N];       // views
  }
  transformed data {
    real<upper=0> neg_log_K;
    neg_log_K <- -log(K); // Equal prior weights
  }
  parameters {
    real<lower=0> alphaq[K]; // cluster means
    real<lower=0> betaq[K];  
    real<lower=0> alphap[K]; // cluster means
    real<lower=0> betap[K];  
  }
  transformed parameters {
    real<upper=0> soft_z[N,K]; // log unnormalized cluster assigns
    for (n in 1:N)
      for (k in 1:K)
        soft_z[n,k] <- neg_log_K + beta_binomial_log(a[n],c[n],alphaq[k],betaq[k])
                                 + beta_binomial_log(c[n],v[n],alphap[k],betap[k]);
  }
  model {
    for (k in 1:K){
      alphaq[k] ~ uniform(0,10000);  // prior??? 
      betaq[k] ~ uniform(0,10000);  // prior???
      alphap[k] ~ uniform(0,10000);  // prior??? 
      betap[k] ~ uniform(0,10000);  // prior???
    }
    for (n in 1:N)
      increment_log_prob(log_sum_exp(soft_z[n])); // likelihood
  }'
  
  stan(model_code = model_code, data = m_dat, iter = 1000, chains = 1)
}


gm2c.test <- function(gm2cfit,t){
  
  m <- extract(gm2cfit)
  mi <- which(max(m$lp__)==m$lp__)[1]
  soft_z <- m$soft_z[mi,,]  
  hard_z <- apply(X=soft_z,MARGIN=1,FUN=function(x){which(x==min(x))})
  
  aq <- m$alphaq[mi,hard_z] + t$a 
  bq <- m$betaq[mi,hard_z] + t$c - t$a
  ap <- m$alphap[mi,hard_z] + t$c
  bp <- m$betap[mi,hard_z] + t$v - t$c
  
  expectedCVR <- (aq/aq+bq)*(ap/ap+bp)
  
  realizedCVR <- t$a / t$v
  risk <- abs(realizedCVR - expectedCVR)
  risk[is.nan(risk)] <- 0
  sum(risk * t$v) / sum(t$v)
}

# TODO consider a zero process as the second process.