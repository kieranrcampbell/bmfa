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
library(magrittr)


# Synthetic ---------------------------------------------------------------

fname <- "data/synthetic.h5"
X <- h5read(fname, "basic_branching/X")
branch <- h5read(fname, "basic_branching/branch_assignment")
true_t <- h5read(fname, "basic_branching/pseudotime")

dp <- data.frame(prcomp(X)$x[,1:2], branch = as.factor(branch), pseudotime = true_t)

ggplot(dp, aes(x = PC1, y = PC2, colour = pseudotime)) + geom_point() + scale_color_viridis()
ggplot(dp, aes(x = PC1, y = PC2, colour = branch)) + geom_point()

source("scripts/gibbs_ard.R")

g <- mfa_gibbs_ard(t(X), iter = 40000, thin = 20)

mc <- to_ggmcmc(g)

ggs_traceplot(filter(mc, Parameter == "lp__")) + stat_smooth()

ggs_autocorrelation(filter(mc, Parameter == "lp__"))


branching_features <- paste0("tau_k[", 1:20, "]")
tau_means <- filter(mc, grepl("tau_k", Parameter)) %>% 
  group_by(Parameter) %>% summarise(map = posterior.mode(mcmc(value))) %>% 
  mutate(is_branching = !(Parameter %in% branching_features))

ggplot(tau_means, aes(x = Parameter, y = map, fill = is_branching)) + 
  geom_bar(stat = "identity") + coord_flip()



tmap <- posterior.mode(mcmc(g$pst_trace))
gamma_mean <- colMeans(g$gamma_trace)

dp2 <- data.frame(prcomp(X)$x[,1:2], branch = gamma_mean, pseudotime = tmap)

ggplot(dp2, aes(x = PC1, y = PC2, colour = pseudotime)) + geom_point() +
  scale_color_viridis(name = "MAP\npseudotime")
ggplot(dp2, aes(x = PC1, y = PC2, colour = branch)) + geom_point() +
  scale_color_viridis(name = "MAP\ngamma")

ggplot(data_frame(branch = as.factor(branch), gamma_mean), 
       aes(x = branch, y = gamma_mean)) + geom_boxplot()


data_frame(true_t, tmap) %>% 
  ggplot(aes(x = true_t, y = tmap)) +
  geom_point(shape = 21, fill = "grey", color = "black", size = 2) +
  ylab("MAP pseudotime") + xlab("True pseudotime") + geom_rug(alpha = 0.5)

filter(mc, Parameter %in% c("k0[2]", "k1[2]", "theta[2]")) %>% 
  ggs_histogram()

filter(mc, Parameter %in% c("k0[26]", "k1[26]", "theta[26]")) %>% 
  ggs_histogram()

filter(mc, Parameter %in% c("c0[1]", "c1[1]", "eta[1]")) %>% 
  ggs_histogram()

filter(mc, Parameter %in% c("c0[26]", "c1[26]", "theta[26]")) %>% 
  ggs_histogram()

branching_features <- paste0("tau_c[", 1:20, "]")
tau_means <- filter(mc, grepl("tau_c", Parameter)) %>% 
  group_by(Parameter) %>% summarise(map = posterior.mode(mcmc(value))) %>% 
  mutate(is_branching = Parameter %in% branching_features)

ggplot(tau_means, aes(x = Parameter, y = map, fill = is_branching)) + 
  geom_bar(stat = "identity") + coord_flip()


tau_map <- posterior.mode(mcmc(g$tau_trace))
tau_k_map <- posterior.mode(mcmc(g$tau_k_trace))

k0 <- posterior.mode(mcmc(g$k0_trace))
k1 <- posterior.mode(mcmc(g$k1_trace))
c0 <- posterior.mode(mcmc(g$c0_trace))
c1 <- posterior.mode(mcmc(g$c1_trace))

plot_expected_behaviour <- function(X, branch, tmap, gamma_mean, k0, k1, c0, c1) {
  d1 <- as_data_frame(X) %>% 
    mutate(pseudotime = tmap, true_branch = as.factor(branch)) %>% 
    gather(feature, expression, -pseudotime, -true_branch)
  
  branch0_mean <- sapply(seq_along(k0), function(g) c0[g] + k0[g] * tmap) %>% as_data_frame()
  branch1_mean <- sapply(seq_along(k1), function(g) c1[g] + k1[g] * tmap) %>% as_data_frame()
  
  branch0_mean %<>% mutate(true_branch = 1, pseudotime = tmap) %>% 
    gather(feature, expression, -true_branch, -pseudotime)
  branch1_mean %<>% mutate(true_branch = 0, pseudotime = tmap) %>% 
    gather(feature, expression, -true_branch, -pseudotime)
  d2 <- bind_rows(branch0_mean, branch1_mean) %>% 
    mutate(true_branch = as.factor(true_branch))
  
  ggplot(d1, aes(x = pseudotime, y = expression, color = true_branch)) +
    geom_point(alpha = 0.5) + scale_color_brewer(palette = "Set1") +
    facet_wrap(~ feature, scales = "free_y") +
    geom_line(data = d2)
  

}



























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

## Visualise what these actually look like across pseudotime
d <- as_data_frame(X) %>% 
  mutate(pseudotime = wishbone_data$Trajectory,
         branch = as.factor(wishbone_data$Branch)) %>% 
  gather(feature, expression, -pseudotime, -branch) %>% 
  mutate(feature = factor(feature, levels = paste0("PC", 1:ncol(X))))

df_branch2 <- df_branch3 <- filter(d, branch == 1)
df_branch2 %<>% mutate(branch = replace(branch, branch == 1, 2))
df_branch3 %<>% mutate(branch = replace(branch, branch == 1, 3))

d %<>% filter(branch != 1)
d <- bind_rows(d, df_branch2, df_branch3)

ggplot(d, aes(x = pseudotime, y = expression, color = branch)) +
  #geom_point(alpha = 0.3) + 
  stat_smooth(se = TRUE, method = "lm") + facet_wrap(~ feature, scales = "free_y") +
  scale_color_brewer(palette = "Set1")

ggplot(d, aes(x = pseudotime, y = expression, color = branch)) +
  geom_point(alpha = 0.1) + 
  stat_smooth(se = FALSE, method = "lm") + facet_wrap(~ feature, scales = "free_y") +
  scale_color_brewer(palette = "Set1")


sc <- sce[rowMeans(exprs(sce) > 0) > 0.2, ] # expressed in at least 20% of cells
Cv2 <- matrixStats::rowVars(exprs(sc)) / rowMeans(exprs(sc))^2

genes <- order(Cv2, decreasing = TRUE)

X <- t(exprs(sc)[genes[1:40], ])

# X <- exprs(sm)
source("scripts/gibbs_constr.R")

whiten <- function(x) (x - mean(x)) / sd(x)
X <- apply(X, 2, whiten)

g <- mfa_gibbs_constr(t(X[,1:2]), iter = 40000, 
                      thin = 5, tau_k = 1, tau_c = 1,
                      pc_initialise = 2, collapse = FALSE)

mc <- to_ggmcmc(g)

ggs_traceplot(filter(mc, Parameter == "lp__")) + stat_smooth()
ggs_autocorrelation(filter(mc, Parameter == "lp__"))


ggs_density(mc, family = "theta")
ggs_density(mc, family = "eta")

# ggs_traceplot(mc, "k0")
# ggs_histogram(mc, "k0")
# ggs_histogram(mc, "k1")

tmap <- posterior.mode(mcmc(g$pst_trace))
gamma_mean <- colMeans(g$gamma_trace)

tsne <- wishbone_data[,1:2]

dp2 <- data.frame(redDim(sce)[,1:2], tsne, branch = gamma_mean, pseudotime = tmap,
                  wishbone_branch = sce$Branch, wishbone_pseudotime = sce$Trajectory) %>% tbl_df()

plot_grid(ggplot(dp2, aes(x = tSNE1, y = tSNE2, colour = pseudotime)) + geom_point() +
            scale_color_viridis(name = "pseudotime"),
          ggplot(dp2, aes(x = tSNE1, y = tSNE2, colour = branch)) + geom_point() +
            scale_color_viridis(name = "branching"))

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

data_frame(gamma_mean, tmap) %>% 
  ggplot(aes(x = tmap, y = gamma_mean)) + geom_point()

df <- data_frame(tmap, cd34, gata1, gata2, mpo, branch = as.factor(round(gamma_mean)))
df <- gather(df, Gene, Expression, -tmap, -branch)

ggplot(df, aes(x = tmap, y = Expression, color = branch)) + geom_point(alpha = 0.2) +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_color_brewer(palette = "Set1") + 
  stat_smooth(se = FALSE, method = "lm") 


