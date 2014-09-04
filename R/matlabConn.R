
# Start and connect to local matlab server
library(R.matlab)
Matlab$startServer()
matlab <- Matlab()
isOpen <- open(matlab)
evaluate(matlab, "addpath '/home/alex/git/bandits/matlab' ")
gi <- evaluate(matlab, "GIBoth(50,100)")
close(matlab)
