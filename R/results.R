
# Plot example data ( 100,000 points is similar to the data )

dm1 <- rm1(num=100000,p=0.001,q=0.6) 
dm1$g <- "dm1"
dm3 <- rm3(num=100000)
dm3$g <- "dm3"
dm5 <- rm5(num=100000,p=c(0.0005,0.00005),q=c(0.8,0.2),theta=0.9)
dm5$g <- "dm5"

ggplot(data=dm1) +
  aes(x = log(a/c), y = log(c/n)) +
  geom_point(aes(size = log(n),alpha = 0.3))

ggplot(data=dm3) +
  aes(x = log(a/c), y = log(c/n)) +
  geom_point(aes(size = log(n),alpha = 0.3))

ggplot(data=dm5) +
  aes(x = log(a/c), y = log(c/n)) +
  geom_point(aes(size = log(n),alpha = 0.3, 
	color=as.factor(c("good","bad"))[clust+1]))

# Correlation tests
corr(rbind(dm1,dm3,dm5))

ggplot(data=d) + 
  aes(x = log(a/c), y = log(c/n)) +
  geom_point(aes(size = log(n),alpha = 0.3, color=as.factor(c("good","bad"))[cluster+1]))



