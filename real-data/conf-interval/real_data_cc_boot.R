# Code to construct CI for clustering coefficient for Facebook Network

library(igraph)
library(parallel)

ncores <- 10 # number of cores for parallel computing

################################################################################
# Function performs bootstrap on the sample
# sam.sub: sample graph of Size n
# f: sampling fraction
# n: sample size
# N: population size
# B: Bootstrap iterations

boot.fun.cl <- function(sam.sub, f, n, N, B){
  
  # construction of the pseudo population
  K <- floor(1/f)
  
  # Adjacency matrix of the sample graph
  An <- as_adjacency_matrix(sam.sub)
  
  boot.rep <- unlist(mclapply(1:B, function(boot_itr){
    if(n*K == N){
      new_pop <<- rep(1:n, K) # N is a multiple of n
    }
    else{
        new_pop <<- c( rep(1:n, K), c(1:n)[1:(N-n*K)] ) # N is not a multiple of n
    }
    # Take a bootstrap sample of size n from N
    boot_sam <- sort(new_pop[sample(1:N, n, replace = F)])
    # Take a bootstrap sample of size n from N
    boot.graph <- graph_from_adjacency_matrix(An[boot_sam, boot_sam])
    # Compute clustering coefficient
    return(transitivity(boot.graph))
  }, mc.cores = ncores))
  # return clustering coefficient of the sample and bootstrap variance
  return(c(transitivity(sam.sub), var(boot.rep)))
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
n <- floor(N*f) # sample size
alpha <- 0.1 # 1 - confidence level
B <- 10 # Bootstrap iteration

res.boot <- matrix(NA, nrow = 3, ncol = length(f)) # matrix to store lower CI, mid point, Upper CI
colnames(res.boot) <- f

for (ii in 1:length(f)) {
  set.seed(1234)
  # select a sample of size n
  sam.v <- sort(sample(V(pop.g), n[ii] , replace = F))
  sam.sub <- subgraph(pop.g, sam.v) # sample subgraph
  res.boot.val <- boot.fun.cl(sam.sub, f[ii], n[ii], N, B)
  res.boot[2,ii] <- res.boot.val[1] # point estimate
  res.boot[c(1,3),ii] <- c((res.boot.val[1] - qnorm(1-alpha/2)*sqrt(res.boot.val[2])),
                           (res.boot.val[1] + qnorm(1-alpha/2)*sqrt(res.boot.val[2]))) # CI
}


# write.csv(res.boot, "boot_cl_ci.csv", row.names = F)

