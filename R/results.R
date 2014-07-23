
source('~/git/bandits/R/Models_gen.R')
source('~/git/bandits/R/Models_test.R')

# Plot example data ( 100,000 points is similar to the data )

dmadv <- adtk.adv()
dm1 <- adtk.m1(num=100000,p=0.001,q=0.6) 
dm3 <- adtk.m3(num=100000)
dm5 <- adtk.m5(num=100000,p=c(0.0005,0.00005),q=c(0.8,0.2),theta=0.9)
dmall <- rbind(dmadv[,-1],dm1,dm3,dm5[,-4])

adtk.plot(dmadv)
adtk.plot(dm1)
adtk.plot(dm3)
adtk.plot(dm5)

# Correlation tests
adtk.corr(dmall)

# Frequentist style tests
ddply(dmall[dmall$c>0,],~g,summarize,
      test1_BinomGoF=adtk.test1(k=a,n=c),         # Binom
      test2_BinomVsBetaBinom=adtk.test2(k=a,n=c), # BetaBinom
      test3_1clustVs2Clust=adtk.test3(k=a,n=c)    # 2 Cluster
      )

# Bayesian style posterior checks

