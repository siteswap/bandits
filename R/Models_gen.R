
# Example datasets to illustrate these techniques.

########################
## Advance campaign data
########################

# Since this is not part of our model, we sample.
site_domain <- read.csv(file="~/site_domain_113803.csv") # TODO

adtk.adv <- function(lid=1023111) {
  
  d <- ddply( site_domain[ site_domain$line_item_id==lid, ],"site_domain",
         summarise,n=sum(imps),c=sum(clicks),a=sum(post_click_convs) )
  d$p <- d$c/d$n
  d$q <- d$a/d$c 
  d$g <- paste("m_",as.character(lid),sep="")
  d$a <- ifelse(d$a > d$c,d$c,d$a) # Set c to be max a.
  d
}


###################
## Model 1 - pooled
###################

adtk.m1 <- function(num,p,q){
  
  n <- sample(site_domain$imps,size=num,replace=TRUE)
  c <- rbinom(num,size=n,prob=p)
  a <- rbinom(num,size=c,prob=q)
  d <- data.frame(a,c,n)
  d$p <- d$c/d$n
  d$q <- d$a/d$c 
  d$g <- "m1"
  d
}


#############################
## Model 3 - hierarchical 
#############################


adtk.m3 <- function(num,s1p,s2p,s1q,s2q){
  
  n <- sample(site_domain$imps,size=num,replace=TRUE)
  c <- rbetabinom.ab(n=num,size=n,shape1=s1p,shape2=s2p)
  a <- rbetabinom.ab(num,size=c,shape1=s1q,shape2=s2q)
  d <- data.frame(a,c,n)
  d$p <- d$c/d$n
  d$q <- d$a/d$c
  d$g <- "m3"
  d
}

##############################
## Model 5 - binomial clusters 
##############################

adtk.m5 <- function(num,p,q,theta){
  
  n <- sample(site_domain$imps,size=num,replace=TRUE)
  clust <- rbinom(n=num,size=1,prob=theta)
  c <- rbinom(num,size=n,prob=p[clust+1]) 
  a <- rbinom(num,size=c,prob=q[clust+1])
  d <- data.frame(a,c,n,clust)
  d$p <- d$c/d$n
  d$q <- d$a/d$c
  d$g <- "m5"
  d
}

##############################
## Model 6 - beta binomial clusters 
##############################

adtk.m6 <- function(num,s1p,s2p,s1q,s2q,theta){
  
  n <- sample(site_domain$imps,size=num,replace=TRUE)
  clust <- rbinom(n=num,size=1,prob=theta)
  c <- rbetabinom.ab(n=num,size=n,shape1=s1p[clust+1],shape2=s2p[clust+1]) 
  a <- rbetabinom.ab(n=num,size=c,shape1=s1q[clust+1],shape2=s2q[clust+1])
  d <- data.frame(a,c,n,clust)
  d$p <- d$c/d$n
  d$q <- d$a/d$c
  d$g <- "m6"
  d
}



##############################
## Util Code
##############################


adtk.plot <- function(d,title){

	gp <- ggplot(data=d) + aes(x = log(a/c), y = log(c/n))
	if(is.null(d$clust)){
	  gp <- gp + geom_point(aes(size = log(n),alpha = 0.3)) + ggtitle(title)
	}else{
	  d$group <- as.factor(c("good","bad"))[d$clust+1]
	  gp <- ggplot(data=d) + aes(x = log(a/c), y = log(c/n))
	  gp <- gp + geom_point(aes(size = log(n),alpha = 0.3, color=group)) +
      ggtitle(title)
	}
  gp
}

