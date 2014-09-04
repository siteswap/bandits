
library(hash)
library(VGAM) 

# t = N-a-b
# Value of arm
v <- function(a,b,t){
  if(t==1){ 
    return( R(a,b) )
  }else{
    
    k <- toString(c(a,b,t))
    if(has.key(k,gi.h)){
      return( unname(values(gi.h,keys=k)) )
    }
    
    q <- c()
    # Arrange the outer loop
    for(m in 1:(t-1)){
      q <- rbind(q,Q(a,b,m,t-m))
    }
    
    r <- sup(a,b,q)
    gi.h[k] <- r
    return(r)   
  }
}

# Q is the weighting we give to each state - weighting
# would be 0 if it's value is s.t. we will never choose
# it. Otherwise it's the probability we will reach it.

# Since we don't know 

# For given prior values, a and b,
# enumerate all possible states, m steps in to future 
# and give the probability of reaching that state (p).
# R is the value for playing that arm once more.
# v is the value for playing with option to keep 
# playing (up to t rounds total).


# TODO - contains a bug - if you manually calculate
# the index for a=3,b=3,t=4 you will see the code
# gives a higher value. This is because we are counting
# path realizations which would have been terminated
# earlier - we can not use betabinomial for value p
# without considering the path taken to reach the state.
Q <- function(a,b,m,t){
  # How likely you are to get this result
  
  # TODO - bug here (see comments).
  # need to remove cases where you would have stopped by now.
  p <- dbetabinom.ab(x=0:m,size=m,shape1=a,shape2=b)
  # Remove cases where you would have stopped by now.
  # Return a pair - the value and the limit for u
  data.frame(p,
             v=mapply(FUN=v,a=a+0:m,b=b+m:0,t=t), 
             R=R(a+0:m,b+m:0), m=m)
}

# Immediate value of the arm
R <- function(a,b){
  # (a+1)/(a+b+2)
  # Gittins (1979) assumes we have flat priors
  # and alpha, beta are wins, losses. Hence 
  # a 'game' would always start off with alpha==0
  # beta==0.
  # But here alpha, beta are the current priors 
  # so we do not add 1.
  a/(a+b)
}

# Supremum function
# Assume a second arm with fixed payoff u,
# we look for the supremum of u values s.t.
# we prefer to take the unknown arm.
# We only need consider 'len' critical 
# values for u.
sup <- function(a,b,q){
  
  len <- length(q$v)
  res <- rep(0,len)
  
  # For each possible value of u (i.e. each v),
  for(i in 1:len){
    # qofu - like Q function defined in Gittins 1979,
    # its the probability of getting to this state and
    # that all previous steps to get here had value greater
    # than u.
    qofu <- q$p*(q$v>=q$v[i])
    # Expected rate of payoff with this u
    res[i] <- (R(a,b) + sum(qofu*q$R)) / (1 + sum(qofu))
  }
  
  max(res)
}

# Let's take the Q function from Gittin's paper and elide the 
# r and u parameters:
# Q(r,a=1,b=1,trials=1,u)

# So, we get a result back for each possible value of r. Which is 
# the probability of the event. This value
# is then 'switched' to zero based on the value function of the 
# state it leads in to.

gi.v <- v

gi.h <- hash()
# clear(gi.h)



#########################
### Profiling GI code ###
#########################

## There may be faster implementations already out there:
#https://sites.google.com/site/lorenzodigregorio/gittins-index
## Nino-Mora:
#http://halweb.uc3m.es/jnino/eng/pubs/2011-ijoc.pdf
#
#clear(gi.h)
#
#> system.time( v(1,1,10) )
#user  system elapsed 
#0.2     0.0     0.2 
#> system.time( v(1,1,20) )
#user  system elapsed 
#2.39    0.00    2.39 
#> system.time( v(1,1,30) )
#user  system elapsed 
#10.48    0.00   10.48 
#> system.time( v(1,1,40) )
#user  system elapsed 
#30.64    0.00   30.66 
#> system.time( v(1,1,50) )
#user  system elapsed 
#71.77    0.32   72.15
#> system.time( v(1,1,60) )
#user  system elapsed 
#145.76    0.03  145.94 
#> system.time( v(1,1,100) )
#user  system elapsed 
#1141.94    0.33 1200.94 
#
## Only computing 10k states, quite slow?
#length(gi.h)
#[1] 9395
#
## Even with GI, it gets expensive:
#plot(c(10*1:6,100),log(c(0.2,2.39,10.48,30.66,72.15,145.94,1200)),
#     ylab="seconds",xlab="horizon",type='l')
#


# 3 issues:
# 1) Slow, single thread implementation in R (linear performance improvement)
# 2) Recomputing for small differences in arm value.
# 3) Exponential increase in compute time.

# Goal 1:
# Get to 100 rounds with 10 arms to compare with ts algos.

# Options:
# 1) Only small improvements can be made in R. Would need to move to java
#    to really improve this (parrallel hashmaps, compiled code etc)
# 2) Look for nearest index value (an unprincipled sort of approximation)


