# Code to construct CI of ratio of triangle density for two networks

library(parallel)
library(igraph)

ncores <- 10 # number of cores for parallel computing

################################################################################
# Function performs bootstrap on the sample
# An: sample network of size n
# f: sampling fraction
# n: sample size
# N: population size

resam.boot.fun.tri <- function(An, n, N, f){
  
  # construction of the pseudo population
  K <- floor(1/f)
  sam.v_boot <- as.numeric(V(An))
  
  if(n*K == N){
    new_pop <<- rep(sam.v_boot, K) # N is a multiple of n
  }
  else{
      new_pop <<- c( rep(sam.v_boot, K), sam.v_boot[1:(N-n*K)] ) # N is not a multiple of n
  }
  # select a sample of size n From N
  boot_sam <- sort(new_pop[sample(1:N, n, replace = F)])
  # Construct bootstrap graph
  boot.graph <- graph_from_adjacency_matrix(get.adjacency(An)[boot_sam, boot_sam])
  # return triangle density for bootstrap sample
  ct <- (length(triangles(boot.graph))/3)/(choose(n,3))
}

################################################################################
# Main function to perform bootstrap
# An: sample network of size n
# f: sampling fraction
# n: sample size
# N: population size
# B: number of bootstrap samples

boot.fun.tri <- function(An, n, N, f, B){
  # compute bootstrap values
  boot_val <- unlist(mclapply(1:B, function(boot_rep){
    resam.boot.fun.tri(An, n, N, f)
  }, mc.cores = ncores))
  # bootstrap estimate of variance (multipled with n)
  boot_var <- n*var(boot_val)
  # return triangle density of the sample and jackknife variance
  return(c((length(triangles(An))/3)/(choose(n,3)), boot_var))
}

################################################################################
# Load the data

unz1 <- unzip("fb-pages-politician.zip", "fb-pages-politician.edges")
dat1 <- read.table(unz1, skip = 1, sep = ",")
pop.g1 <- graph_from_data_frame(dat1, directed = FALSE)


unz2 <- unzip("fb-pages-public-figure.zip", "fb-pages-public-figure.edges")
dat2 <- read.table(unz2, skip = 1, sep = ",")
pop.g2 <- graph_from_data_frame(dat2, directed = FALSE)

################################################################################
# Implementation

n <- c(2000,2500,3000) # sample size
alpha <- 0.1 # 1 - confidence level
B <- 10 # number of boostrap samples

N1 <- length(V(pop.g1)) # population size of network 1
f1 <- n/N1 # sampling fraction for network 1

N2 <- length(V(pop.g2)) # population size of network 2
f2 <- n/N2 # sampling fraction for network 2

test_res_boot_mat <- matrix(NA, nrow = 3, ncol = length(n)) # matrix to store lower CI, mid point, Upper CI
colnames(test_res_boot_mat) <- n

for (ii in 1:length(n)){
  set.seed(2345)
  # select sample of n from each network
  sam.v1 <- sort(sample(V(pop.g1), n[ii], replace = F))
  sam.v2 <- sort(sample(V(pop.g2), n[ii], replace = F))

  sam.sub1 <- subgraph(pop.g1, sam.v1)
  sam.sub2 <- subgraph(pop.g2, sam.v2)

  res.tri1_boot <- boot.fun.tri(sam.sub1, n[ii], N1, f1[ii], B)
  res.tri2_boot <- boot.fun.tri(sam.sub2, n[ii], N2, f2[ii], B)

  v_boot <- sqrt(c(res.tri1_boot[2], res.tri2_boot[2])) # sd estimate vector
  u_boot <- c(res.tri1_boot[1], res.tri2_boot[1]) # point estimate vector
  
  mu_boot <- u_boot/v_boot # Standardized vector: estimator/sd
  
  test_res_boot_mat[2,ii] <-  u_boot[1]/u_boot[2] # ratio of estimate
  
  denom.g_boot <- sqrt( (1/(mu_boot[2]^2)) + ((mu_boot[1])^2/(mu_boot[2]^4))) # denominator of the delta method
  
  test_res_boot_mat[c(1,3),ii] <- c( u_boot[1]/u_boot[2] - qnorm(1-alpha/2)*(v_boot[1]*denom.g_boot)/(sqrt(n[ii])*v_boot[2]), 
                                     u_boot[1]/u_boot[2] + qnorm(1-alpha/2)*(v_boot[1]*denom.g_boot)/(sqrt(n[ii])*v_boot[2])) # CI
}

# write.csv(test_res_boot_mat, "test_res_boot_mat.csv", row.names = F)
