library(plyr)
library(ggplot2)

# Note, campaign_id is always 3597062
d <- read.csv(file="~/site_domain_113803.csv")
d <- d[ d$line_item_id==1023111, c("site_domain","imps","clicks","post_click_convs") ]
      
# Use only levers with atleast 10 appearances
freq <- ddply(d,"site_domain",summarise,l=length(imps) )
freqDomains <- freq[freq$l >= 10,1]
d <- d[d$site_domain %in% freqDomains,]

# Randomize the order of the rows, and order by site_domain
set.seed(seed=2105) # Set seed so sample is consistent
d <- d[sample(1:dim(d)[1]),]

# Take half the data for training, half for testing
halfOne <- function(l){ sum(l[1:ceiling(length(l)/2)]) }
halfTwo <- function(l){ sum(l[(ceiling(length(l)/2)+1):length(l)]) }

train <- ddply(d,"site_domain",summarise,v=halfOne(imps),c=halfOne(clicks),a=halfOne(post_click_convs) )
test <- ddply(d,"site_domain",summarise,v=halfTwo(imps),c=halfTwo(clicks),a=halfTwo(post_click_convs) )
dirty <- function(t){ ((t$c < t$a) + (t$v < t$c)) != 0 }
dirtyData <- union(which(dirty(train)),which(dirty(test)))
train <- train[-dirtyData,]
test <- test[-dirtyData,]

# Test away !

gmp.test(gmp(train),test)
# -1.332663e-05

gmie.test(gmie(train),test)
# -8.668857e-06

gm1.test(gm1(train),test)
# -9.292773e-06 # warnings

gm2.test(gm2(train),test)
# -8.773567e-06 # warnings

# Model 1 with mixture.
gm1c.test(gm1c(train),test)

# Model 3 - something wrong.
gm2c.test(gm2c(train),test)

# Looks bimodal:
plot(density(log(train$a/train$v)))





