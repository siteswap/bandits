
# Example datasets to illustrate these techniques.

###################
## Model 1 - pooled
###################

# Since this is not part of our model, we sample.
setwd("~")
site_domain <- read.csv(file="site_domain_113803.csv")

rm1 <- function(num,p,q){
  
  n <- sample(site_domain$imps,size=num,replace=TRUE)
  c <- rbinom(num,size=n,prob=p)
  a <- rbinom(num,size=c,prob=q)
  d <- data.frame(a,c,n)
  d$p <- d$c/d$n
  d$q <- d$a/d$c 
  d
}


#############################
## Model 3 - hierarchical 
#############################


rm3 <- function(num){
  
  n <- sample(site_domain$imps,size=num,replace=TRUE)
  c <- rbetabinom.ab(n=num,size=n,shape1=1,shape2=10000)
  a <- rbetabinom.ab(num,size=c,shape1=6,shape2=4)
  d <- data.frame(a,c,n)
  d$p <- d$c/d$n
  d$q <- d$a/d$c
  d
}

##############################
## Model 5 - binomial clusters 
##############################

rm5 <- function(num,p,q,theta){
  
  n <- sample(site_domain$imps,size=num,replace=TRUE)
  clust <- rbinom(n=num,size=1,prob=theta)
  c <- rbinom(num,size=n,prob=p[clust+1]) 
  a <- rbinom(num,size=c,prob=q[clust+1])
  d <- data.frame(a,c,n,clust)
  d$p <- d$c/d$n
  d$q <- d$a/d$c
  d
}

