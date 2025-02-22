# Code to construct CI for clustering coefficient for Facebook Network

library(igraph)

################################################################################
# Function performs jackknife on the sample
# sam.sub: sample graph of size n
# f: sampling fraction

jack.fun.cl <- function(sam.sub, f){
  
  jack.sam <- unlist(lapply(V(sam.sub), function(i){
    
    # create the sampled network deleting vertex i
    ints <- setdiff(V(sam.sub),i)
    sub.graph <- subgraph(sam.sub, ints)
    
    # Compute clustering coefficient
    ct <- (transitivity(sub.graph))
    return(ct)
  }))
  # jackknife variance
  jack.var <- (1-f)*sum((jack.sam - mean(jack.sam))^2)
  # return clustering coefficient of the sample and jackknife variance
  return(c( transitivity(sam.sub), jack.var))
}



################################################################################
# Load the data

unzip("fb-pages-tvshow.zip", list = TRUE)
unz <- unzip("fb-pages-tvshow.zip", "fb-pages-tvshow.edges")
# readLines(unz, n=10)
dat <- read.table(unz, skip=1, sep=",") # skip first line to avoid % bipartite unweighted" 

# # look
# head(dat)
# str(dat)

# Load as a graph
pop.g <- graph_from_data_frame(dat, directed = FALSE)

################################################################################
# Implementation

pop.val <- transitivity(pop.g) # 0.5906502

N <- length(V(pop.g)) # population size
f <- c(0.05, 0.1, 0.2, 0.3) # sampling fraction
n <- round(N*f) # sample size
alpha <- 0.1 # 1 - confidence level

res.jack <- matrix(NA, nrow = 3, ncol = length(f)) # matrix to store lower CI, mid point, Upper CI
colnames(res.boot) <- f

for (ii in 1:length(f)) {
  set.seed(1234)
  # select a sample of size n
  sam.v <- sort(sample(V(pop.g), n[ii] , replace = F))
  sam.sub <- subgraph(pop.g, sam.v) # sample subgraph
  res.jack.val <- jack.fun.cl(sam.sub, f[ii])
  res.jack[2,ii] <- res.jack.val[1] # point estimate
  res.jack[c(1,3),ii] <- c((res.jack.val[1] - qnorm(1-alpha/2)*sqrt(res.jack.val[2])),
                           (res.jack.val[1] + qnorm(1-alpha/2)*sqrt(res.jack.val[2]))) # CI
}


# write.csv(res.jack, "jack_cl_ci.csv")
