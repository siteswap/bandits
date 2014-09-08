
library(plyr)
library(ggplot2)
library(xtable

# Weighted correlation function
wcor <- function(x,y,w){
  x <- x - mean(x)
  y <- y - mean(y)
  corr <- (x*y)/(sd(x)*sd(y))
  sum(corr*w)/sum(w)
}


# We take the top 3 lineitems by total acquisitions and look only at where we have clicks.
# We also drop lineitem 767753 as it only has acquisitions through facebook (no good for correlations analysis).
# setwd("./ProjectBandits")
d <- data <- read.csv(file="site_domain_113803.csv")
acqs <- ddply(d,~line_item_id,summarize,C=sum(post_click_convs))
acqs <- acqs[order(acqs$C,decreasing=TRUE)[1:3],1]
acqs <- acqs[-which(acqs == 767753)] 
d <- d[d$line_item_id %in% acqs,]
d <- ddply(d,c("line_item_id","site_domain"),summarise,n=sum(imps),c=sum(clicks),a=sum(post_click_convs) )
d$q <- d$a/d$c
d$p <- d$c/d$n
d$r <- d$a/d$n
d$abin <- 1*(d$a==0) # for logistic regression later
d$line_item_id <- as.factor(d$line_item_id)
d <- d[d$c>0,]

       
ggplot(d, aes(x = log(q), y = log(p))) +
  facet_grid(line_item_id ~ .) +
  geom_point() + ggtitle("CTR vs CVR by line_item_id") +
  theme(plot.title = element_text(lineheight=.8, face="bold"))
  
# We consider both correlation and weighted correlation where the weight is the number of impressions. The idea here is that we want the many low impression sites not to obscure some relationship in the more important sites. Also, the lower the clicks, the higher the noise will be.
        
tmp <- ddply(d,~line_item_id,summarize,Cor=cor(a/c,c/n),WeigthedCor=wcor(a/c,c/n,n),items=length(n),imps=sum(n))
xtable(tmp)
        
# We see no correlation, so we strip out the groups with no acquisitions to see if we can find a relationship in this subgroup. Hopefully this can also strip out some noise.

        
tmp <- ddply(d[d$a>0,],~line_item_id,summarize,Cor=cor(a/c,c/n),WeigthedCor=wcor(a/c,c/n,n),items=length(n),imps=sum(n))
xtable(tmp)        
        
# And also strip out the group with q is exactly 1.

tmp <- ddply(d[(d$a>0) & (d$q!=1),],   
              ~line_item_id,summarize,Cor=cor(a/c,c/n),WeigthedCor=wcor(a/c,c/n,n),items=length(n),imps=sum(n))
xtable(tmp)      
        
# When we break the data down by campaign, we see that much of the apparent correlation is due clusters within lineitems:

        # Look at the data by campaign_id and line_item_id
        # We can see that both are important
        d <- data <- read.csv(file="site_domain_113803.csv")
        acqs <- ddply(d,~line_item_id,summarize,C=sum(post_click_convs))
        acqs <- acqs[order(acqs$C,decreasing=TRUE)[1:3],1]
        acqs <- acqs[-which(acqs == 767753)] 
        d <- d[d$line_item_id %in% acqs,]
        d <- ddply(d,c("line_item_id","campaign_id","site_domain"),summarise,n=sum(imps),c=sum(clicks),a=sum(post_click_convs) )
        d$q <- d$a/d$c
        d$p <- d$c/d$n
        d$r <- d$a/d$n
        d$abin <- 1*(d$a>0) # for logistic regression later
        d$line_item_id <- as.factor(d$line_item_id)
        d$campaign_id <- as.factor(d$campaign_id)
        d <- d[d$c>0,]
        
        
        ggplot(d, 
          aes(x = log(q), y = log(p))) +
          facet_grid(campaign_id ~ .) +
          geom_point() + ggtitle("CTR vs CVR by campaign") +
          theme(plot.title = element_text(lineheight=.8, face="bold"))
        
tmp <- ddply(d[(d$a>0) ,],c("line_item_id","campaign_id"),
             summarize,Cor=cor(q,p),WeigthedCor=wcor(q,p,n),items=length(n),imps=sum(n))
xtable(tmp)  

        
# The previous section discussed the clustering of q around 0 and 1. Here we look at whether we could predict which cluster it will enter. 
# This could lead to a compound model.      
        
# We remove sites with 10 or fewer impressions to reduce noise.         
        
xtable(
  summary( glm(formula=abin~log(p)+campaign_id, data=d[d$n > 100,],family=binomial) )
  )
        
# This result makes sense - we see consistently lower CTR for the no acqs sites:
t <- ddply(d, c("line_item_id","campaign_id","abin"),summarize,CTR=100*sum(c)/sum(n))
xtable(t)
        
# Plot this - suggests the opposite?!!!
        d$abinf <- as.factor(d$abin)
        ggplot(d[(d$c > 0) & (d$n > 0),], aes(x=log(p), fill=abinf)) + 
          geom_density(alpha=.3) +
          facet_grid(campaign_id ~ .)
        
        
        d$abinf <- as.factor(d$abin)
        ggplot(d[(d$c > 0),], aes(x=log(c), fill=abinf)) + 
          geom_density(alpha=.3) +
          facet_grid(campaign_id ~ .)
        
        
        d[ (d$c>0)  & (d$campaign_id==3683130 ),]

        
