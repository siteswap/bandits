
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
  risk <- realizedCVR - gmpfit
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
  risk <- realizedCVR - gmiefit
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
  risk <- realizedCVR - gm1fit$alpha_p/(gm1fit$alpha_p + gm1fit$beta_p)
  
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
  risk <- realizedCVR - (gm2fit$alpha_p1/(gm2fit$alpha_p1 + gm2fit$beta_p1))*(gm2fit$alpha_q1/(gm2fit$alpha_q1 + gm2fit$beta_q1))
  
  # Volume weighted mean loss:
  risk[is.nan(risk)] <- 0
  sum(risk * t$v) / sum(t$v)
}


################
# Model 1 Mixture
################

# No option to override prior
gm1m <- function(d){
  
  gm1m_dat <- list(J = dim(d)[1], 
                     a = d$a,
                     v = d$v,
                     K = 2)

  model1_mix_code <- '
    data {
      int<lower=0> J; // Levers
      int<lower=0> a[J]; // acquisitions
      int<lower=0> v[J]; // views
      int<lower=1> K; // groups
    }
    parameters {
      // Can we do something with beta-binomial to save computation?
      simplex[K] theta;
      real<lower=0> alpha[K]; 
      real<lower=0> beta[K];
      real<lower=0,upper=1> q[J];
    }
    model {
      real ps[K]; // hold the k values for summation
      for(j in 1:J){
        for(k in 1:K){
          // the k need to be summed together so we must use log_sum_exp
          ps[k] <- binomial_log(a[j],v[j],q[j]) + 
                        beta_log(q[j],alpha[k],beta[k]) + log(theta[k]);
        }
        increment_log_prob(log_sum_exp(ps));
      }    
    }
    '
  
  # Error : Error in function boost::math::lgamma<d>(d): numeric overflow
  # error occurred during calling the sampler; sampling not done
  # TODO - chains don't run in parrallel?!?!?
  # https://groups.google.com/forum/#!msg/stan-users/3goteHAsJGs/UmSDRo77Vn4J
  # Could use beta-binomial if you don't need q estimates at first.
  stan(model_code = model1_mix_code, data = gm1m_dat, iter = 2500, chains = 1)
}

gm1m.test <- function(gm1mfit,t){
  
  q <- extract(gm1mfit)$q
  realizedCVR <- t$a / t$v
  # Loss is difference in values, risk is expected loss
  risk <- c()
  for(j in 1:length(realizedCVR)){
    risk[j] <- realizedCVR[j] - mean(q[,j])
  }
  
  # Volume weighted mean loss:
  risk[is.nan(risk)] <- 0
  sum(risk * t$v) / sum(t$v)
}


#################
# Model 1 Cluster
#################


# Test data
# k1 <- rbinom(n=120,size=100,prob=rbeta(n=12,shape1=1,shape2=10))
# k2 <- rbinom(n=80,size=100,prob=rbeta(n=8,shape1=8,shape2=3))
# 
# gm1c_dat <- list(N = 200, 
#                  K = 2,
#                  a = c(k1,k2),
#                  v = rep(100,200)
#                  )


# No option to override prior
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
  
  stan(model_code = model1_clust_code, data = gm1c_dat, iter = 1000, chains = 1)
  
}

# Running this on the training data, we see 2 
# groupings emerge.

#                   mean se_mean     sd    2.5%     25%     50%     75%   97.5% n_eff Rhat
# alpha[1]           1.3     0.0    0.2     0.9     1.1     1.3     1.4     1.7    40  1.1
# alpha[2]           0.1     0.0    0.0     0.1     0.1     0.1     0.1     0.2    17  1.1
# beta[1]         8542.6   209.1 1212.0  5547.1  7841.8  8877.7  9518.3  9974.7    34  1.0
# beta[2]          154.4    11.0   73.7    55.3    98.4   131.7   213.4   306.6    45  1.1

# TODO - plot the sites' log(v) vs log(a) and colour by mean mixing weight.


