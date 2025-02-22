# Code to compute the ratio of bootstrap variance to the true variance

library(parallel)
library(igraph)
library(graphon)

detectCores()
n.core <- 10 # number of codes to be used in parallel


################################################################################
# Function computes the estimate of true population variance using MC
# AN: population graph of Size N
# N: population size
# n: sample size
# B: MC replication

true.fun <- function(AN,N, n, B){
  return(unlist(lapply(1:B, function(i){
    # Take sample of size n from N
    b.sam <- sort(sample(1:N, n, replace = F))
    # Compute triangle density
    return((length(triangles(subgraph(AN, b.sam)))/3)/(choose(n,3)))
  })))
}

################################################################################
# Function performs boostrap on the sample
# AN: population graph of Size N
# n: sample size
# N: population size
# sam: sampled vertices

boot.fun <- function(AN, sam, n, N){
  
  # construction of the pseudo population
  K <- floor(N/n)
  
  if(n*K == N){
    new_pop <<- rep(sam, K) # N is a multiple of n
  }
  else{
      new_pop <<- c( rep(sam, K), sam[1:(N-n*K)] ) # N is not a multiple of n
  }
  # Take a boostrap sample of size n from N
  boot_sam <- sort(new_pop[sample(1:N, n, replace = F)])
  # Graph from boostrap samples
  boot.graph <- graph_from_adjacency_matrix(get.adjacency(AN)[boot_sam, boot_sam])
  # Compute triangle density 
  ct <- (length(triangles(boot.graph))/3)/(choose(n,3))
  return(ct)
}

################################################################################
# Function computes boostrap varince for the sample
# AN: population graph of Size N
# n: sample size
# sam: sampled vertices
# f: sampling fraction
# B: number of boostrap replication

resam.fun <- function(AN, sam, n, f, B){
  
  st <- Sys.time()
  # subgraph density from B boostrap samples
  boot_val <- unlist(lapply(1:B, function(boot_rep){
    boot.fun(AN, sam, n, N)
  }))
  # Boostrap variance
  boot_var <- var(boot_val)
  boot_time <- as.numeric(Sys.time()-st)/60
  
  # return boostrap estimate of variance and time taken
  return(c(boot_var, boot_time))
}

################################################################################
# Main function for jackknife
# AN: population graph of Size N
# n: sample size
# f : sampling fraction
# B : number of boostrap replication

get.val <- function(AN, N, f, B){
  n <- round(N*f)
  sam <- sort(sample(1:N, n, replace = FALSE))
  val <- resam.fun(AN, sam, n, f, B)
  return(c(val))
}

################################################################################
# Simulation parameters

N <- 100 # population size
f <- c(0.05, 0.1, 0.3) # sampling fractions
n <- round(N*f) # sample size
B <- 10 # Bootstrap replicate to estimate true variance
M <- 10 # MC replicate for bootstrap


################################################################################
# Random Graph Generating Models

#########################################
# W random Graph
# Kernel: |eps_i - eps_j|

Xi <- matrix(rep(runif(N), times = N), nrow = N, byrow = F) # Latent positions
H <- (abs(Xi - t(Xi))) # Matrix of connection probabilities
AN <- graph_from_adjacency_matrix(gmodel.P(H, rep = 1, noloop = TRUE, symmetric.out = TRUE)) # population graph

#########################################
# SBM

th.mat <- matrix(c(0.4, 0.1, 0.2, 0.1, 0.5, 0.3, 0.2, 0.3, 0.7), ncol = 3) # Block probability matrix
lambda <- c(0.3, 0.3, 0.4) # Probability of group membership
b <- N*lambda
AN <- sample_sbm(N, th.mat, b, directed = F, loops = F) # population graph

#########################################
# RDPG

norm_vec <- function(x) x/sqrt(sum(x^2)) # Norm of a vector

d <- 6 # Dimension of latent positions
X <- matrix(runif(N*d), nrow = N, byrow = F) # Latent positions
nX <- apply(X, 1, norm_vec)
nX <- t(nX) %*% nX 
diag(nX) <- 0 # Matrix of connection probabilities
AN <- graph_from_adjacency_matrix(gmodel.P(nX, rep = 1, noloop = TRUE, symmetric.out = TRUE))

#########################################
# Sparse Graphon Model
# Kernel: rho*|eps_i - eps_j|

rho <- 1/N^(1/12) # Sparsity parameter
Xi <- matrix(rep(runif(N), times = N), nrow = N, byrow = F) # Latent positions
H <- rho*(abs(Xi - t(Xi))) # Matrix of connection probabilities
AN <- graph_from_adjacency_matrix(gmodel.P(H, rep = 1, noloop = TRUE, symmetric.out = TRUE)) # population graph

#########################################
# Sparse SBM

rho <- 1/N^(1/12) # Sparsity parameter
th.mat <- rho*matrix(c(0.4, 0.1, 0.2, 0.1, 0.5, 0.3, 0.2, 0.3, 0.7), ncol = 3) # Block probability matrix
lambda <- c(0.3, 0.3, 0.4) # Probability of group membership
b <- N*lambda
AN <- sample_sbm(N, th.mat, b, directed = F, loops = F) # population graph

################################################################################
# Implementation 

val.df.list <- list() # Data frame to store bootstrap values 
val.df.true <- numeric() # store true variance

for (f_i in 1:length(f)) {
  
  n_i <- round(N*f[f_i])
  # Store values for each sampling fraction
  
  val.df.list[[f_i]] <- matrix( unlist( mclapply(1:M, function(j) {
    get.val(AN, N, f[f_i], B)}, mc.cores = n.core) ) , byrow = T, ncol = 2) # bootstrap values
  
  true.sam <- true.fun(AN, N, n_i, B)
  val.df.true[f_i] <- var(true.sam) # population variance (not multiplied with n) 
}

################################################################################
# Visualization

res.df <- matrix(unlist(lapply(1:length(f), function(f_i){
  val.df.list[[f_i]][,1]/val.df.true[f_i] # ratio: bootstrap var/ true var
})), byrow = F, ncol = length(f))
colnames(res.df) <- f

boxplot(res.df, ylim = c( min(1, min(res.df)), max(res.df) ))
abline(h = 1, col = 'red', lwd = 2)

