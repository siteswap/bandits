

L <- 3 # levers
# matrix of params
mp <- data.frame(a=rep(1,L),b=rep(1,L))

# Probability of getting result on lever l, given params
p <- function(l,mp){
  mp[l,"a"]/(mp[l,"a"]+mp[l,"b"])
}

# Quality of lever l, with t trials ahead.
q <- function(l,mp,t){
  
  q1 <- p(l,mp) # 1 step quality
  
  if(t==0){
    
    return( q1 ) # Immediate value of lever l
    
  }else{
    
    # Quality under success
    mp$a[l]  <- mp$a[l] + 1
    # sq <- max(q(1,mp,t-1),q(2,mp,t-1))
    sq <- max(sapply(1:L,FUN=q,mp=mp,t=t-1))
    
    # Quality under failure
    mp$a[l]  <- mp$a[l] - 1
    mp$b[l]  <- mp$b[l] + 1
    fq <- max(sapply(1:L,FUN=q,mp=mp,t=t-1))
    
    return( q1 + q1*sq + (1 - q1)*fq )
  }
}

q.all <- function(mp,t){
  sapply(1:L,FUN=q,mp=mp,t=t)/(t+1)
}


### Experiments ###

# When both arms have same mean, chooses higher variance.
> mp
a b
1 2 2
2 1 1
> q.all(mp,1)
[1] 0.5250000 0.5416667


# As time horizon goes up, net benefit goes up.
> diff(q.all(mp,5))*6
[1] -0.02452547


# It prefers exploitation (arm 2) when rounds are low
# but as future rounds increase it prefers exploration (arm 1)
> mp
a b
1 2 2
2 6 5
> q.all(mp,0)
[1] 0.5000000 0.5454545 # Exploit arm 2
> q.all(mp,1)
[1] 0.5363636 0.5454545 # Exploit arm 2
> q.all(mp,2)
[1] 0.5515152 0.5530303 # Exploit arm 2
> q.all(mp,3)
[1] 0.5590909 0.5568182 # Explore arm 1
> q.all(mp,4)
[1] 0.5646753 0.5611389 # Explore arm 1


# 1) Compare this against other algos in simulation.
# 2) See if we might approximate it for large values of t and L.

# Note in a single call, we have played out every possible game
# So we would only ever need to call the function once, if we can
# utilise a nice cache of pre-aggregated values.

# A bit of plotting - could we estimate these curves more cheaply?

getVals <- function (mp,maxt) {
  vals <- matrix(data=rep(0,maxt*3),nrow=maxt)
  for(r in 0:(maxt-1)){
    qvals <- q.all(mp,t=r)
    vals[r+1,] <- qvals / sum(qvals) # sum(qvals)
  }
  vals
}

# Note the outcome of Arm 1 *DOES* change the order of arm2 and arm3
# This is not true for Gittins indexes.

par(mfrow=c(2,3))

mp <- data.frame(a=c(1,4,7),b=c(1,3,5))
vals <- getVals(mp,maxt=8)
plot(vals[,1],type='l',ylim=c(min(vals),max(vals)),xlim=c(1,8),
     main="Index value of each arm \n Prior(a=c(1,4,7),b=c(1,3,5))")
lines(vals[,2],col="red")
lines(vals[,3],col="blue")
abline(v=1)

# Assume we play black, under success:
mp <- data.frame(a=c(2,4,7),b=c(1,3,5))
valss1 <- getVals(mp,maxt=7)
plot(x=2:8,valss1[,1],type='l',ylim=c(min(valss1),max(valss1)),xlim=c(1,8),
     main="Under success \n Prior(a=c(2,4,7),b=c(1,3,5))")
lines(x=2:8,valss1[,2],col="red")
lines(x=2:8,valss1[,3],col="blue")
abline(v=2)

# Under failure:
mp <- data.frame(a=c(1,4,7),b=c(2,3,5))
valsf1 <- getVals(mp,maxt=7)
plot(x=2:8,valsf1[,1],type='l',ylim=c(min(valsf1),max(valsf1)),xlim=c(1,8),
     main="Under failure \n Prior(a=c(1,4,7),b=c(2,3,5))")
lines(x=2:8,valsf1[,2],col="red")
lines(x=2:8,valsf1[,3],col="blue")
abline(v=2)

# Play red, under success:
mp <- data.frame(a=c(1,5,7),b=c(2,3,5))
valsf1s2 <- getVals(mp,maxt=6)
plot(x=3:8,valsf1s2[,1],type='l',ylim=c(min(valsf1s2),max(valsf1s2)),xlim=c(1,8),
     main="Under success \n Prior(a=c(1,5,7),b=c(2,3,5))")
lines(x=3:8,valsf1s2[,2],col="red")
lines(x=3:8,valsf1s2[,3],col="blue")
abline(v=3)

# Play red, under failure:
mp <- data.frame(a=c(1,4,7),b=c(2,4,5))
valsf1f2 <- getVals(mp,maxt=6)
plot(x=3:8,valsf1f2[,1],type='l',ylim=c(min(valsf1f2),max(valsf1f2)),xlim=c(1,8),
     main="Under failure \n Prior(a=c(1,4,7),b=c(2,4,5))")
lines(x=3:8,valsf1f2[,2],col="red")
lines(x=3:8,valsf1f2[,3],col="blue")
abline(v=3)



