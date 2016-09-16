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

Xz <- X
Xz[Xz < 0] <- 0

dp <- data.frame(prcomp(Xz)$x[,1:2], branch = as.factor(branch), pseudotime = true_t)

ggplot(dp, aes(x = PC1, y = PC2, colour = pseudotime)) + geom_point() +
  scale_color_viridis(name = "Pseudotime")

# ggsave("figs/for_presentation/synthetic_pst.png", width = 5, height = 3)

ggplot(dp, aes(x = PC1, y = PC2, colour = branch)) + geom_point() +
  scale_color_brewer(palette = "Set1", name = "Branch")

ggsave("figs/for_presentation/synthetic_branching.png", width = 4.5, height = 3)


source("houdini/R/zi.R")
sourceCpp("houdini/src/gibbs.cpp")

set.seed(1L)
g <- mfa_zi(Xz, iter = 2e5, lambda = 1,
             thin = 1e2, collapse = T, eta_tilde = mean(Xz))

#source("scripts/gibbs_semi_ard.R")
#g <- mfa_gibbs_semi_ard(t(X), iter = 200000, thin = 100, collapse = TRUE, eta_tilde = 5)

mc <- to_ggmcmc(g)

qplot(as.vector(Xz), as.vector(imputed)) + xlab("Measured") + ylab("Imputed")

ggs_traceplot(filter(mc, Parameter == "lp__")) + stat_smooth()

ggs_autocorrelation(filter(mc, Parameter == "lp__"))





tmap <- posterior.mode(mcmc(g$pst_trace))
gamma_mean <- colMeans(g$gamma_trace)

dp2 <- data.frame(prcomp(Xz)$x[,1:2], branch = gamma_mean, pseudotime = tmap)

plt_map_pseudotime <- ggplot(dp2, aes(x = PC1, y = PC2, colour = pseudotime)) + 
  geom_point(size = 2) +
  scale_color_viridis(name = expression(paste("MAP ", t))) + theme_classic()
plt_map_pseudotime

plt_map_branch <- ggplot(dp2, aes(x = PC1, y = PC2, colour = branch)) + 
  geom_point(size = 2) +
  scale_color_viridis(name = expression(paste("MAP ", gamma))) + theme_classic()
plt_map_branch

plt_true_map <- data_frame(true_t, tmap) %>% 
  ggplot(aes(x = true_t, y = tmap)) +
  geom_point(shape = 21, fill = "grey", color = "black", size = 2) +
  ylab("MAP pseudotime") + xlab("True pseudotime") + geom_rug(alpha = 0.5) +
  theme_classic()
plt_true_map

print(cor(true_t, tmap, method = "spearman"))

tau_map <- filter(mc, grepl("tau[", Parameter, fixed =T)) %>% 
  group_by(Parameter) %>% summarise(map = mean(value)) %>% 
  mutate(is_branching = !(Parameter %in% branching_features))

branching_features <- paste0("tau_k[", 1:20, "]")
chi_means <- filter(mc, grepl("tau_k", Parameter)) %>% 
  group_by(Parameter) %>% summarise(map = 1 / mean(value)) %>% # 1 / posterior.mode(mcmc(value))) %>%  #
  mutate(is_branching = !(Parameter %in% branching_features)) %>% 
  mutate(map = map * tau_map$map)

chi_means %<>% arrange(desc(map)) %>% 
  mutate(Parameter = as.character(Parameter)) %>% 
  mutate(Parameter = factor(Parameter, levels = Parameter))

plt_chi_branch <- ggplot(chi_means, aes(x = Parameter, y = map, fill = is_branching)) + 
  geom_bar(stat = "identity") + coord_flip() + #  scale_y_sqrt() +
  scale_fill_brewer(name = "gene branches?", palette = "Set1") +
  ylab(expression(paste("[MAP ", chi[g] ,"]" ^ "-1"))) + xlab("Gene") + theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
plt_tau_branch
# ggsave("figs/ard.png", width = 5, height = 3)

plot_grid(plt_map_pseudotime, plt_map_branch,
          plt_true_map, plt_tau_branch, labels = "AUTO")

ggsave("~/Google Drive/campbell/mfa/manuscript/figs/toy.png", width = 8, height = 5)
ggsave("~/Google Drive/campbell/mfa/nips/figs/nips_fig_1.png", width = 8, height = 5)

ggsave("figs/for_presentation/synthetic_results.png", width = 11, height = 7)












# Wishbone ----------------------------------------------------------------




raw <- read_csv("data/wishbone_mouse_marrow_scrnaseq.csv")

set.seed(1L)
cells_to_sample <- seq_len(nrow(raw)) # sample(seq_len(nrow(raw)), 2000)

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
is_exprs(sce) <- exprs(sce) > 0
sce <- calculateQCMetrics(sce)

sc <- sce[rowMeans(exprs(sce) > 0) > 0.2, ] # expressed in at least 20% of cells

## find genes correlated with trajectory

means <- rowMeans(exprs(sc))
vars <- matrixStats::rowVars(exprs(sc))
to_use <- vars > 5 # can change to 5 too

s <- sc[to_use, ]
plotPCA(s, colour_by = "Branch")
plotPCA(s, colour_by = "Trajectory")

Y <- t(exprs(s))


set.seed(123L)
lambda <- empirical_lambda(t(exprs(sce)))

grna <- mfa_zi(Y, iter = 5e4, thin = 25, collapse = F,
                pc_initialise = 2, lambda = lambda)

mc <- to_ggmcmc(grna)

qplot(as.vector(Y), as.vector(grna$x_mean_trace)) + xlab("Measured") + ylab("Imputed")


ggs_traceplot(filter(mc, Parameter == "lp__")) + stat_smooth()
ggs_autocorrelation(filter(mc, Parameter == "lp__"))


# MAP pseudotime plots ---

tmap <- posterior.mode(mcmc(grna$pst_trace))
gamma_mean <- colMeans(grna$gamma_trace)
branch <- gamma_mean # round(gamma_mean)

tsne <- wishbone_data[,1:2]

dp2 <- data.frame(tsne, branch = gamma_mean, pseudotime = tmap,
                  wishbone_branch = sce$Branch, wishbone_pseudotime = sce$Trajectory) %>% tbl_df()
  
plot_grid(ggplot(dp2, aes(x = tSNE1, y = tSNE2, colour = pseudotime)) + geom_point() +
            scale_color_viridis(name = expression(paste("MAP ", t))),
          ggplot(dp2, aes(x = tSNE1, y = tSNE2, colour = branch)) + geom_point() +
            scale_color_viridis(name = expression(paste("MAP ", gamma))),
          nrow = 2, labels = "AUTO")

map_plt <- last_plot()

# ggsave("figs/for_presentation/wishbone_tsne.png", width = 10, height = 4)
# ggsave("~/Google Drive/campbell/pseudogp/mix-fa/manuscript/figs/wishbone_tsne.png",
#        width = 10, height = 4)

data_frame(tmap, trajectory = sc$Trajectory) %>% 
  ggplot(aes(x = trajectory, y = tmap)) + geom_point(size = 2, alpha = 0.5)


# Gene plots --- 

tau_means <- filter(mc, grepl("tau_k", Parameter)) %>% 
  group_by(Parameter) %>% summarise(map = 1 / mean(value)) %>% 
  mutate(Gene = featureNames(s)) %>% 
  arrange(desc(map)) %>% 
  mutate(Gene = as.character(Gene)) %>% 
  mutate(Gene = factor(Gene, levels = Gene))

ggplot(tau_means, aes(x = Gene, y = map)) + 
  geom_bar(stat = "identity") + coord_flip() +
  ylab(expression(paste("[MAP ", chi[g] ,"]" ^ "-1"))) + xlab("Gene") +
  theme(axis.text.y = element_text(size = 6))

chi_plot <- last_plot()

plot_tsne_expression <- function(gene) {
  d <- as_data_frame(tsne) %>% mutate(expr = exprs(sce)[gene, ])
  
  ggplot(d, aes(x = tSNE1, y = tSNE2, color = log2(expr + 1))) +
    geom_point() + scale_color_viridis(name = gene)
}

chi_plot

gplt <- plot_grid(plot_tsne_expression("ELANE"),
          plot_tsne_expression("CAR2"),
          plot_tsne_expression("RPL26"), 
          nrow = 3, labels = c("B", "C", "D"))

# plot_grid(chi_plot, gplt, nrow = 1, labels = c("A", ""))
# 
# ggsave("figs/for_presentation/wishbone_ard.png", width = 8, height = 8)
# ggsave("~/Google Drive/campbell/pseudogp/mix-fa/manuscript/figs/wishbone_ard.png",
#        width = 8, height = 8)


##### NIPS plot #######

gplt <- plot_grid(plot_tsne_expression("ELANE"),
                  plot_tsne_expression("CAR2"),
                  plot_tsne_expression("RPL26"), 
                  nrow = 3, labels = c("D", "", ""))

# lower_grid <- plot_grid(chi_plot, gplt, ncol = 2, labels = c("C", ""))

plot_grid(map_plt, chi_plot, gplt, ncol = 3, labels = c("", "C", ""))
ggsave("figs/nips_fig_2.png", width = 11, height = 7)
ggsave("~/Google Drive/campbell/mfa/nips/figs/nips_fig_2.png", width = 11, height = 7)

###########

dx <- as_data_frame(prcomp(X)$x[,1:2]) %>% mutate(tmap, gamma_mean) 

ggplot(dx, aes(x = PC1, y = PC2, color = tmap)) + geom_point() +
  scale_color_viridis()

ggplot(dx, aes(x = PC1, y = PC2, color = gamma_mean)) + geom_point() +
  scale_color_viridis()

# plot the expected behaviour

k0 <- posterior.mode(mcmc(g$k0_trace))
k1 <- posterior.mode(mcmc(g$k1_trace))
c0 <- posterior.mode(mcmc(g$c0_trace))
c1 <- posterior.mode(mcmc(g$c1_trace))

#plot_expected_behaviour <- function(X, branch, tmap, gamma_mean, k0, k1, c0, c1) {
d1 <- as_data_frame(X) %>% 
  mutate(pseudotime = tmap, true_branch = as.factor(branch)) %>% 
  gather(feature, expression, -pseudotime, -true_branch)

branch0_mean <- sapply(seq_along(k0), function(g) c0[g] + k0[g] * tmap) %>% as_data_frame()
branch1_mean <- sapply(seq_along(k1), function(g) c1[g] + k1[g] * tmap) %>% as_data_frame()

names(branch0_mean) <- names(branch1_mean) <- colnames(X)

branch0_mean %<>% mutate(true_branch = 0, pseudotime = tmap) %>% 
  gather(feature, expression, -true_branch, -pseudotime)
branch1_mean %<>% mutate(true_branch = 1, pseudotime = tmap) %>% 
  gather(feature, expression, -true_branch, -pseudotime)
d2 <- bind_rows(branch0_mean, branch1_mean) %>% 
  mutate(true_branch = as.factor(true_branch))

ggplot(d1, aes(x = pseudotime, y = expression, color = true_branch)) +
  geom_point(alpha = 0.5) + scale_color_brewer(palette = "Set1") +
  facet_wrap(~ feature, scales = "free_y") +
  geom_line(data = d2)



# What genes are we actually using ----------------------------------------


