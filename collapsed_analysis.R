library(cowplot)
library(ggplot2)
library(rhdf5)
library(coda)
library(MCMCglmm)
library(viridis)
library(ggmcmc)
library(dplyr)
library(readr)
library(scater)
library(matrixStats)
library(tidyr)


# Synthetic ---------------------------------------------------------------

fname <- "data/synthetic.h5"
X <- h5read(fname, "basic_branching/X")
branch <- h5read(fname, "basic_branching/branch_assignment")
true_t <- h5read(fname, "basic_branching/pseudotime")

dp <- data.frame(prcomp(X)$x[,1:2], branch = as.factor(branch), pseudotime = true_t)

ggplot(dp, aes(x = PC1, y = PC2, colour = pseudotime)) + geom_point() + scale_color_viridis()
ggplot(dp, aes(x = PC1, y = PC2, colour = branch)) + geom_point()

source("scripts/gibbs_constr.R")

g <- mfa_gibbs_constr(t(X))

mc <- to_ggmcmc(g)

ggs_traceplot(filter(mc, Parameter == "lp__")) + stat_smooth()

ggs_traceplot(filter(mc, Parameter == "theta[1]")) + stat_smooth()
ggs_traceplot(filter(mc, Parameter == "theta[2]")) + stat_smooth()


tmap <- posterior.mode(mcmc(g$pst_trace))
gamma_mean <- colMeans(g$gamma_trace)

dp2 <- data.frame(prcomp(X)$x[,1:2], branch = gamma_mean, pseudotime = tmap)

ggplot(dp2, aes(x = PC1, y = PC2, colour = pseudotime)) + geom_point() +
  scale_color_viridis(name = "MAP\npseudotime")
ggplot(dp2, aes(x = PC1, y = PC2, colour = branch)) + geom_point() +
  scale_color_viridis(name = "MAP\ngamma")


data_frame(true_t, tmap) %>% 
  ggplot(aes(x = true_t, y = tmap)) +
  geom_point(shape = 21, fill = "grey", color = "black", size = 2) +
  ylab("MAP pseudotime") + xlab("True pseudotime") + geom_rug(alpha = 0.5)




# Wishbone ----------------------------------------------------------------




raw <- read_csv("data/wishbone_mouse_marrow_scrnaseq.csv")

set.seed(1L)
cells_to_sample <- sample(seq_len(nrow(raw)), 500)

tpm <- as.matrix(raw[cells_to_sample,-1])

nc <- ncol(tpm)
wishbone_data <- data.frame(tpm[, (nc - 3):nc])
tpm <- tpm[, 1:(nc - 4)]

cellnames <- raw[[1]][cells_to_sample]
pd <- data.frame(sample = cellnames, Trajectory = wishbone_data$Trajectory, Branch = wishbone_data$Branch)
rownames(pd) <- cellnames
tpm <- t(tpm)
colnames(tpm) <- cellnames

sce <- newSCESet(exprsData = tpm, phenoData = new("AnnotatedDataFrame", pd))
sce <- calculateQCMetrics(sce)

plotPCA(sce, colour_by = "Branch")
plotPCA(sce, colour_by = "Trajectory")

cd34 <- exprs(sce)["CD34", ]
gata1 <- exprs(sce)["GATA1", ]
gata2 <- exprs(sce)["GATA2", ]
mpo <- exprs(sce)["MPO", ]

save(sce, file = "data/wishbone.Rdata")

markers <- c("CD34", "GATA1", "GATA2", "MPO")
sm <- sce[match(markers, featureNames(sce)), ]

sce <- plotPCA(sce, colour_by = "MPO", return_SCESet = TRUE)
plotPCA(sce, colour_by = "CD34")
plotPCA(sce, colour_by = "GATA1")

# X <- prcomp(t(exprs(sce)))$x[,1:20] # redDim(sce)
X <- redDim(sce)


sc <- sce[rowMeans(exprs(sce) > 0) > 0.2, ] # expressed in at least 20% of cells
Cv2 <- matrixStats::rowVars(exprs(sc)) / rowMeans(exprs(sc))^2

genes <- order(Cv2, decreasing = TRUE)

# X <- exprs(sm)
source("scripts/gibbs_constr.R")
g <- mfa_gibbs_constr(t(X[,1:10]), iter = 30000, thin = 15, tau_k = 5, tau_c = 5)

mc <- to_ggmcmc(g)

ggs_traceplot(filter(mc, Parameter == "lp__")) + stat_smooth()

# ggs_traceplot(mc, "k0")
# ggs_histogram(mc, "k0")
# ggs_histogram(mc, "k1")

tmap <- posterior.mode(mcmc(g$pst_trace))
gamma_mean <- colMeans(g$gamma_trace)

tsne <- wishbone_data[,1:2]

dp2 <- data.frame(redDim(sce)[,1:2], tsne, branch = gamma_mean, pseudotime = tmap,
                  wishbone_branch = sce$Branch, wishbone_pseudotime = sce$Trajectory) %>% tbl_df()

t1 <- plot_grid(ggplot(dp2, aes(x = tSNE1, y = tSNE2, colour = pseudotime)) + geom_point() +
            scale_color_viridis(name = "mfa"),
          ggplot(dp2, aes(x = tSNE1, y = tSNE2, colour = wishbone_pseudotime)) + geom_point() +
            scale_color_viridis(name = "wishbone"))


t2 <- plot_grid(ggplot(dp2, aes(x = tSNE1, y = tSNE2, colour = branch)) + geom_point() +
            scale_color_viridis(name = "mfa"),
          ggplot(dp2, aes(x = tSNE1, y = tSNE2, colour = wishbone_branch)) + geom_point() +
            scale_color_viridis(name = "wishbone"))

plot_grid(t1, t2, labels = c("pseudotime", "branching"), nrow = 2)
ggsave("~/Desktop/mfa.png")

# plot_grid(ggplot(dp2, aes(x = PC1, y = PC2, colour = pseudotime)) + geom_point() +
#             scale_color_viridis(name = "mfa"),
#           ggplot(dp2, aes(x = PC1, y = PC2, colour = wishbone_pseudotime)) + geom_point() +
#             scale_color_viridis(name = "wishbone"))
          

ggsave("~/Desktop/mfa_pst.png",width=8,height=4)

plot_grid(ggplot(dp2, aes(x = PC1, y = PC2, colour = branch)) + geom_point() +
            scale_color_viridis(name = "mfa"),
          ggplot(dp2, aes(x = PC1, y = PC2, colour = wishbone_branch)) + geom_point() +
            scale_color_viridis(name = "wishbone"))


ggsave("~/Desktop/mfa_branching.png",width=8,height=4)


ggplot(dp2, aes(x = wishbone_pseudotime, y = pseudotime)) + geom_point()
ggplot(dp2, aes(x = as.factor(wishbone_branch), y = branch)) + geom_boxplot()

df <- data_frame(tmap, cd34, gata1, gata2, mpo, branch = as.factor(round(gamma_mean)))
df <- gather(df, Gene, Expression, -tmap, -branch)

ggplot(df, aes(x = tmap, y = Expression, color = branch)) + geom_point(alpha = 0.2) +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_color_brewer(palette = "Set1") + 
  stat_smooth(se = FALSE, method = "loess") 


