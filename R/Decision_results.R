
initVals <- data.frame(q=rbeta(arms,shape1=5,shape2=10),p=rbeta(arms,shape1=5,shape2=5))
ts_clicks <- as.list(rep(0,4))
ts_acqs <- as.list(rep(0,4))

for(v in 1:4){  
  # 1. one only tries to maximise number of clicks, ignoring conversions
  ts_clicks[[v]] <- adtk.mab(DFUN=adtk.ts_clicks,trueVals=initVals)
  # 2. one only tries to maximise number of conversions, ignoring clicks
  ts_acqs[[v]] <- adtk.mab(DFUN=adtk.ts_acqs,trueVals=initVals)
}

par(mfrow=c(4,2),oma=c(0,0,2,0))
for(v in 1:4){
  adtk.mabplot(ts_clicks[[v]])
  adtk.mabplot(ts_acqs[[v]])
}
title("UCB max clicks - UCB max acqs",outer=TRUE)  

# 3. one tries to maximise number of conversions, considering clicks as well, assuming however, 
#     that the two rates are independent (formulate the model carefully and extract a confidence interval for UCB). 