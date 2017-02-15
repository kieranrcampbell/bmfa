library(rhdf5)



# MH - gibbs --------------------------------------------------------------

source("scripts/mh-gibbs.R")

fname <- "data/synthetic.h5"
X <- h5read(fname, "basic_branching/X")

delta_true <- h5read(fname, "basic_branching/delta")

pc1 <- prcomp(X)$x[,1]
pst_init <- (pc1 - min(pc1) + 0.01) / (max(pc1) - min(pc1) + 0.02)

# fixed <- list(pst = true_t, delta0 = delta_true[,1], delta1 = delta_true[,2])
fixed <- list(delta0 = delta_true[,1], delta1 = delta_true[,2])#, gamma = branch)#,

proposals = list(kappa_t = 0.1, kappa_k = 0.8, kappa_delta = 0.1); set.seed(123)
g <- mh_gibbs(X, iter = 20000, thin = 10, proposals = proposals, 
              fixed = fixed, ptype = "individual", pst_init = pst_init)

save(g, file = "data/clusterg.Rdata")
