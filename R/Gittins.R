
# GI has (t^2)/2 order.
# Hash values.

Q <- function(a,b,m,t){
  # print(paste("Q",a,b,m,t))
  # How likely you are to get this result
  p <- dbinom(x=0:m,size=m,prob=R(a,b))
  # Remove cases where you would have stopped by now.
  # Return a pair - the value and the limit for u
  data.frame(p,
             v=mapply(FUN=v,a=a+0:m,b=b+m:0,t=t), 
             R=R(a+0:m,b+m:0), m=m)
}

R <- function(a,b){
  (a+1)/(a+b+2)
}

# Supremum function
sup <- function(a,b,q){
  
  len <- length(q$v)
  res <- rep(0,len)
  
  for(i in 1:len){
    qofu <- q$p*(q$v>=q$v[i])
    res[i] <- (R(a,b) + sum(qofu*q$R)) / (1 + sum(qofu))
  }
  
  max(res)
}

# t = N-a-b
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

gi.h <- hash()
# clear(gi.h)


