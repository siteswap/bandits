
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


adtk.m3 <- function(num){
  
  n <- sample(site_domain$imps,size=num,replace=TRUE)
  c <- rbetabinom.ab(n=num,size=n,shape1=1,shape2=10000)
  a <- rbetabinom.ab(num,size=c,shape1=6,shape2=4)
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
## Util Code
##############################


adtk.plot <- function(d){

	gp <- ggplot(data=d) + 
	  aes(x = log(a/c), y = log(c/n))
	if(is.null(d$clust)){
	  gp + geom_point(aes(size = log(n),alpha = 0.3))
	}else{
	  gp + geom_point(aes(size = log(n),alpha = 0.3, color=as.factor(c("good","bad"))[clust+1]))
	}
}

