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

ggplot(dp, aes(x = PC1, y = PC2, colour = pseudotime)) + geom_point() +
  scale_color_viridis(name = "Pseudotime")

ggsave("figs/for_presentation/synthetic_pst.png", width = 5, height = 3)

ggplot(dp, aes(x = PC1, y = PC2, colour = branch)) + geom_point() +
  scale_color_brewer(palette = "Set1", name = "Branch")

ggsave("figs/for_presentation/synthetic_branching.png", width = 4.5, height = 3)


source("houdini/R/houdini.R")
sourceCpp("houdini/src/gibbs.cpp")

set.seed(1L)
g <- mfa_cpp(t(scale(X)), iter = 5e5, 
             thin = 5e2, collapse = T, eta_tilde = 0)

#source("scripts/gibbs_semi_ard.R")
#g <- mfa_gibbs_semi_ard(t(X), iter = 200000, thin = 100, collapse = TRUE, eta_tilde = 5)

mc <- to_ggmcmc(g)


ggs_traceplot(filter(mc, Parameter == "lp__")) + stat_smooth()

ggs_autocorrelation(filter(mc, Parameter == "lp__"))


branching_features <- paste0("tau_k[", 1:20, "]")
tau_means <- filter(mc, grepl("tau_k", Parameter)) %>% 
  group_by(Parameter) %>% summarise(map = 1 / mean(value)) %>% #posterior.mode(mcmc(value))) %>% 
  mutate(is_branching = !(Parameter %in% branching_features))
tau_means %<>% arrange(desc(map)) %>% 
  mutate(Parameter = as.character(Parameter)) %>% 
  mutate(Parameter = factor(Parameter, levels = Parameter))

plt_tau_branch <- ggplot(tau_means, aes(x = Parameter, y = map, fill = is_branching)) + 
  geom_bar(stat = "identity") + coord_flip() + #  scale_y_sqrt() +
  scale_fill_brewer(name = "gene branches?", palette = "Set1") +
  ylab(expression(paste("[MAP ", chi[g] ,"]" ^ "-1"))) + xlab("Gene") + theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
# ggsave("figs/ard.png", width = 5, height = 3)


tmap <- posterior.mode(mcmc(g$pst_trace))
gamma_mean <- colMeans(g$gamma_trace)

dp2 <- data.frame(prcomp(X)$x[,1:2], branch = gamma_mean, pseudotime = tmap)

plt_map_pseudotime <- ggplot(dp2, aes(x = PC1, y = PC2, colour = pseudotime)) + geom_point() +
  scale_color_viridis(name = expression(paste("MAP ", t))) + theme_classic()

plt_map_branch <- ggplot(dp2, aes(x = PC1, y = PC2, colour = branch)) + geom_point() +
  scale_color_viridis(name = expression(paste("MAP ", gamma))) + theme_classic()


plt_true_map <- data_frame(true_t, tmap) %>% 
  ggplot(aes(x = true_t, y = tmap)) +
  geom_point(shape = 21, fill = "grey", color = "black", size = 2) +
  ylab("MAP pseudotime") + xlab("True pseudotime") + geom_rug(alpha = 0.5) +
  theme_classic()

print(cor(true_t, tmap, method = "spearman"))

plot_grid(plt_map_pseudotime, plt_map_branch,
          plt_true_map, plt_tau_branch, labels = "AUTO")
ggsave("~/Google Drive/campbell/pseudogp/mix-fa/manuscript/figs/toy.png", width = 11, height = 7)

ggsave("figs/for_presentation/synthetic_results.png", width = 11, height = 7)

filter(mc, Parameter %in% c("k0[29]", "k1[29]", "theta[29]", "eta[1]",
                            "tau_k[29]", "c0[29]", "c1[29]")) %>% 
  ggs_traceplot()




filter(mc, Parameter %in% c("eta[1]", "eta[2]")) %>% ggs_histogram()

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
  

#}



























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

sce <- plotPCA(sce, colour_by = "Branch", return_SCESet = TRUE)
ggsave("figs/for_presentation/hpscs.png", width = 6, height = 6)


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
  gather(feature, expression, -pseudotime, -branch) #%>% 
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
# source("scripts/gibbs_semi_ard.R")

g <- mfa_gibbs_semi_ard(t(X[,1:2]), iter = 40000, thin = 20, collapse = TRUE,
                        pc_initialise = 2, eta_tilde = 0)


mc <- to_ggmcmc(g)

ggs_traceplot(filter(mc, Parameter == "lp__")) + stat_smooth()
ggs_autocorrelation(filter(mc, Parameter == "lp__"))


tau_means <- filter(mc, grepl("tau_k", Parameter)) %>% 
  group_by(Parameter) %>% summarise(map = mean(value)) %>% 
  arrange(desc(map)) %>% 
  mutate(Parameter = as.character(Parameter)) %>% 
  mutate(Parameter = factor(Parameter, levels = Parameter))

ggplot(tau_means, aes(x = Parameter, y = map)) + 
  geom_bar(stat = "identity") + coord_flip() + scale_y_sqrt() +
  scale_fill_brewer(name = "gene branches?", palette = "Set1") +
  ylab(expression(paste(MAP, sqrt(tau)))) + xlab("Gene") 

tmap <- posterior.mode(mcmc(g$pst_trace))
gamma_mean <- colMeans(g$gamma_trace)

tsne <- wishbone_data[,1:2]

dp2 <- data.frame(redDim(sce)[,1:2], tsne, branch = gamma_mean, pseudotime = tmap,
                  wishbone_branch = sce$Branch, wishbone_pseudotime = sce$Trajectory) %>% tbl_df()

plot_grid(ggplot(dp2, aes(x = tSNE1, y = tSNE2, colour = pseudotime)) + geom_point() +
            scale_color_viridis(name = "pseudotime"),
          ggplot(dp2, aes(x = tSNE1, y = tSNE2, colour = branch)) + geom_point() +
            scale_color_viridis(name = "branching"))



# 
# t1 <- plot_grid(ggplot(dp2, aes(x = tSNE1, y = tSNE2, colour = pseudotime)) + geom_point() +
#             scale_color_viridis(name = "mfa"),
#           ggplot(dp2, aes(x = tSNE1, y = tSNE2, colour = wishbone_pseudotime)) + geom_point() +
#             scale_color_viridis(name = "wishbone"))


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


ggplot(dp2, aes(x = wishbone_pseudotime, y = pseudotime, color = as.factor(round(branch)))) + 
  geom_point() + scale_color_brewer(name = "Branch", palette = "Set1") +
  xlab("Wishbone pseudotime") + ylab("Our pseudotime") +
  geom_rug()

ggplot(dp2, aes(x = as.factor(wishbone_branch), y = branch)) + geom_boxplot()

data_frame(gamma_mean, tmap) %>% 
  ggplot(aes(x = tmap, y = gamma_mean)) + geom_point()

df <- data_frame(tmap, cd34, gata1, gata2, mpo, branch = as.factor(round(gamma_mean)))
df <- gather(df, Gene, Expression, -tmap, -branch)

ggplot(df, aes(x = tmap, y = Expression, color = branch)) + geom_point(alpha = 0.2) +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_color_brewer(palette = "Set1") + 
  stat_smooth(se = FALSE, method = "lm") 



# 

# Wishbone analysis here --------------------------------------------------



sc <- sce[rowMeans(exprs(sce) > 0) > 0.2, ] # expressed in at least 20% of cells

## find genes correlated with trajectory

gcors <- apply(exprs(sc), 1, function(x) cor(x, sc$Trajectory))
bcors <- apply(exprs(sc), 1, function(x) {
  a <- aov(x ~ as.factor(sc$Branch))
  summary(a)[[1]][["Pr(>F)"]][1]
})


means <- rowMeans(exprs(sc))
vars <- matrixStats::rowVars(exprs(sc))
to_use <- vars > 8 # can change to 5 too

# to_use_g <- abs(gcors) > 0.35
# 
# to_use <- to_use_g | to_use_b

s <- sc[to_use, ]
plotPCA(s, colour_by = "Branch")
plotPCA(s, colour_by = "Trajectory")

Y <- exprs(s)
X <- t(Y) 

X <- scale(X) # , scale = FALSE)


set.seed(123L)
grna <- mfa_cpp(t(X), iter = 50000, thin = 25, collapse = FALSE,
                      pc_initialise = 2, eta_tilde = 0)


mc <- to_ggmcmc(grna)

ggs_traceplot(filter(mc, Parameter == "lp__")) + stat_smooth()
ggs_autocorrelation(filter(mc, Parameter == "lp__"))


# MAP pseudotime plots ---

tmap <- posterior.mode(mcmc(grna$pst_trace))
gamma_mean <- colMeans(grna$gamma_trace)
branch <- gamma_mean # round(gamma_mean)

tsne <- wishbone_data[,1:2]

dp2 <- data.frame(redDim(sce)[,1:2], tsne, branch = gamma_mean, pseudotime = tmap,
                  wishbone_branch = sce$Branch, wishbone_pseudotime = sce$Trajectory) %>% tbl_df()
  
plot_grid(ggplot(dp2, aes(x = tSNE1, y = tSNE2, colour = pseudotime)) + geom_point() +
            scale_color_viridis(name = expression(paste("MAP ", t))),
          ggplot(dp2, aes(x = tSNE1, y = tSNE2, colour = branch)) + geom_point() +
            scale_color_viridis(name = expression(paste("MAP ", gamma))),
          labels = "AUTO")

map_plt <- last_plot()

ggsave("figs/for_presentation/wishbone_tsne.png", width = 10, height = 4)
ggsave("~/Google Drive/campbell/pseudogp/mix-fa/manuscript/figs/wishbone_tsne.png",
       width = 10, height = 4)

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

plot_grid(chi_plot, gplt, nrow = 1, labels = c("A", ""))

ggsave("figs/for_presentation/wishbone_ard.png", width = 8, height = 8)
ggsave("~/Google Drive/campbell/pseudogp/mix-fa/manuscript/figs/wishbone_ard.png",
       width = 8, height = 8)


##### NIPS plot #######

gplt <- plot_grid(plot_tsne_expression("ELANE"),
                  plot_tsne_expression("CAR2"),
                  plot_tsne_expression("RPL26"), 
                  nrow = 3, labels = c("D", "", ""))

lower_grid <- plot_grid(chi_plot, gplt, ncol = 2, labels = c("C", ""))

plot_grid(map_plt, lower_grid, nrow = 2, rel_heights = c(1, 2))
ggsave("figs/nips_fig_2.png", width = 8, height = 11)
ggsave("~/Google Drive/campbell/mfa/nips/figs/nips_fig_2.png", width = 8, height = 11)

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


