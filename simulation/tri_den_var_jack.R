# Code to compute the ratio of jackknife variance to the true variance

library(parallel)
library(igraph)
library(graphon)

detectCores()
n.core <- 10 # number of codes to be used in parallel

################################################################################
# Function computes the estimate of true population variance using MC
# AN: population graph of Size N
# n: sample size
# B: MC replication

true.fun <- function(AN, n, B){
  return(unlist(lapply(1:B, function(i){
    # Take sample of size n from N
    b.sam <- sort(sample(1:N, n, replace = F))
    # Compute triangle density
    return((length(triangles(subgraph(AN, b.sam)))/3)/(choose(n,3)))
  })))
}

################################################################################
# Function performs jackknife on the sample
# AN: population graph of Size N
# n: sample size
# sam: sampled vertices

jack.fun <- function(AN, sam, n){
  unlist(lapply(sam, function(i){
    # create the sampled network deleting vertex i
    ints <- setdiff(sam,i)
    sub.graph <- subgraph(AN, ints)
    # Compute triangle density 
    ct <- (length(triangles(sub.graph))/3)/(choose(n-1,3))
    return(ct)
  }))
}

################################################################################
# Function computes jackknife varince for the sample
# AN: population graph of Size N
# n: sample size
# sam: sampled vertices
# f : sampling fraction

resam.fun <- function(AN, sam, n, f){
  
  st <- Sys.time()
  # compute the 1-deleted subgraph denisty counts
  jack.sam <- jack.fun(AN, sam ,n)
  # jackknife varaince
  jack.var <- (1-f)*(n-1)*sum((jack.sam - mean(jack.sam))^2)
  jack.time <- as.numeric(Sys.time()-st)/60
  
  # resturn jackknife estimate of varaince and time taken
  return(c(jack.var, jack.time))
}

################################################################################
# Main function for jackknife
# AN: population graph of Size N
# n: sample size
# f : sampling fraction
# B : MC replicate for jackknife

get.val <- function(AN, N, f, B){
  n <- round(N*f)
  sam <- sort(sample(1:N, n, replace = FALSE))
  val <- resam.fun(AN, sam, n, f)
  return(c(val))
}


################################################################################
# Simulation paraemters

N <- 100 # population size
f <- c(0.05, 0.1, 0.3) # sampling fractions
n <- round(N*f) # sample size
B <- 10 # MC replicate to eatimate true variance
M <- 10 # MC replicate for jackknife

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

val.df.list <- list() # Data frame to store jackknife values 
val.df.true <- numeric() # store true variance

for (f_i in 1:length(f)) {
  # Store values for each sampling fraction
  
  val.df.list[[f_i]] <- matrix( unlist( mclapply(1:M, function(j) {
    get.val(AN, N, f[f_i], B)}, mc.cores = n.core) ) , byrow = T, ncol = 2) # jackknife values
  
  n_i <- round(N*f[f_i])
  true.sam <- true.fun(AN, n_i, B)
  val.df.true[f_i] <- n_i*var(true.sam) # population variance
}

################################################################################
# Visualization

res.df <- matrix(unlist(lapply(1:length(f), function(f_i){
  val.df.list[[f_i]][,1]/val.df.true[f_i] # ratio: jackknife var/ true var
})), byrow = F, ncol = length(f))
colnames(res.df) <- f

boxplot(res.df, ylim = c( min(1, min(res.df)), max(res.df) ))
abline(h = 1, col = 'red', lwd = 2)

