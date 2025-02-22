# Code to construct CI of ratio of triangle density for two networks

library(igraph)

################################################################################
# Function performs jackknife on the sample
# sam.sub: sample graph of size n
# f: sampling fraction

jack.fun.tri <- function(sam.sub, f){
  n <- length(V(sam.sub))
  
  jack.sam <- unlist(lapply(V(sam.sub), function(i){
    # create the sampled network deleting vertex i
    ints <- setdiff(V(sam.sub),i)
    sub.graph <- subgraph(sam.sub, ints)
    
    # Compute triangle density
    ct <- (length(triangles(sub.graph))/3)/choose(n-1,3)
    return(ct)
  }))
  # jackknife variance
  jack.var <- (n-1)*(1-f)*sum((jack.sam - mean(jack.sam))^2)
  # return triangle density of the sample and jackknife variance
  return(c( (length(triangles(sam.sub))/3)/choose(n,3), jack.var))
}

################################################################################
# load the data

unz1 <- unzip("fb-pages-politician.zip", "fb-pages-politician.edges")
dat1 <- read.table(unz1, skip = 1, sep = ",") # skip first line to avoid % bipartite unweighted" 
pop.g1 <- graph_from_data_frame(dat1, directed = FALSE)


unz2 <- unzip("fb-pages-public-figure.zip", "fb-pages-public-figure.edges")
dat2 <- read.table(unz2, skip = 1, sep = ",") # skip first line to avoid % bipartite unweighted" 
pop.g2 <- graph_from_data_frame(dat2, directed = FALSE)


################################################################################
# Implementation

n <- c(2000,2500,3000) # sample size
alpha <- 0.1 # 1 - confidence level

N1 <- length(V(pop.g1)) # population size of network 1
f1 <- n/N1 # sampling fraction for network 1

N2 <- length(V(pop.g2)) # population size of network 2
f2 <- n/N2 # sampling fraction for network 2

test_res_jack_mat <- matrix(NA, nrow = 3, ncol = length(n)) # matrix to store lower CI, mid point, Upper CI
colnames(test_res_jack_mat) <- n

for (ii in 1:length(n)){
  set.seed(2345)
  # select sample of n from each network
  sam.v1 <- sort(sample(V(pop.g1), n[ii] , replace = F))
  sam.v2 <- sort(sample(V(pop.g2), n[ii] , replace = F))
  
  sam.sub1 <- subgraph(pop.g1, sam.v1)
  sam.sub2 <- subgraph(pop.g2, sam.v2)
  
  res.tri1_jack <- jack.fun.tri(sam.sub1, f1[ii])
  res.tri2_jack <- jack.fun.tri(sam.sub2, f2[ii])
  
  v_jack <- sqrt(c(res.tri1_jack[2], res.tri2_jack[2])) # sd estimate vector
  u_jack <- c(res.tri1_jack[1], res.tri2_jack[1]) # point estimate vector
  
  mu_jack <- u_jack/v_jack # Standardized vector: estimator/sd
  
  test_res_jack_mat[2,ii] <- u_jack[1]/u_jack[2] # ratio of estimate
  
  denom.g_jack <- sqrt( (1/(mu_jack[2]^2)) + ((mu_jack[1])^2/(mu_jack[2]^4))) # denominator of the delta method
  
  test_res_jack_mat[c(1,3),ii] <- c( u_jack[1]/u_jack[2] - qnorm(1-alpha/2)*(v_jack[1]*denom.g_jack)/(sqrt(n[ii])*v_jack[2]), 
                                     u_jack[1]/u_jack[2] + qnorm(1-alpha/2)*(v_jack[1]*denom.g_jack)/(sqrt(n[ii])*v_jack[2])) # CI 
}

# write.csv(test_res_jack_mat, "test_res_jack_mat.csv", row.names = F)




