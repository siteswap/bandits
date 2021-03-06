
library(xtable)

source('~/git/bandits/R/Models_gen.R')
source('~/git/bandits/R/Models_test.R')

# Plot example data ( 100,000 points is similar to the data )

dmadv <- adtk.adv()
dm1 <- adtk.m1(num=100000,p=0.001,q=0.6) 
dm3 <- adtk.m3(num=100000,s1p=1,s2p=10000,s1q=6,s2q=4)
dm5 <- adtk.m5(num=100000,p=c(0.0005,0.00005),q=c(0.8,0.2),theta=0.9)
dm6 <- adtk.m6(num=100000,s1p=c(5,5),s2p=c(10000,100000),
                          s1q=c(8,2),s2q=c(10,10),
                          theta=0.9)
dmall <- rbind(dmadv[,-1],dm1,dm3,dm5[,-4],dm6[,-4])

adtk.plot(dmadv,"Advanced data")
p1 <- adtk.plot(dm1,"Model 1 \n p=0.001,q=0.6")
p3 <- adtk.plot(dm3,"Model 3 \n s1p=1,s2p=10000,\ns1q=6,s2q=4")
p5 <- adtk.plot(dm5,"Model 5 \n p=c(0.0005,0.00005),\nq=c(0.8,0.2),\ntheta=0.9")
p6 <- adtk.plot(dm6,"Model 6 \n s1p=c(5,5),s2p=c(10000,100000),\n s1q=c(8,2),s2q=c(10,10),\ntheta=0.9")
# grid.arrange(p1,p3,p5,p6, ncol = 2)

# Correlation tests
results_corr <- adtk.corr(dmall)

# Diagnostic testsk
results_diag <- ddply(dmall[dmall$c>0,],~g,summarize,
      test1=adtk.test1(k=a,n=c,grpt=10),# Binom
      test2=adtk.test2(k=a,n=c), # BetaBinom
      test3_30=adtk.test3(k=a,n=c,q=.3),
      test3_45=adtk.test3(k=a,n=c,q=.45),         # 2 Cluster Binom
      test3_60=adtk.test3(k=a,n=c,q=.6)           # 2 Cluster Binom
      #,test4_2ClustBetaBinom=adtk.test4(k=a,n=c) # 2 Cluster BetaBinom
      )

# Residual plots:
grp <- ifelse(n > grpt, 1:length(n), rep(-1, length(n)) )
df <- data.frame(grp,k,n)
D <- ddply(df, ~grp, summarize, k=sum(k), n=sum(n) )
p <- sum(k)/sum(n) # ML estimation of rate
r  <- (D$k - D$n*p)/(D$k*p*(1-p))
plot(sort(r),type='l',ylim=c(-2.5,2.5))
par(new=TRUE)
plot(qnorm(p=(1:159)/160,mean=0,sd=1),type='l',ylim=c(-2.5,2.5))

# Print results
xtable(results_corr)
xtable(results_diag)

