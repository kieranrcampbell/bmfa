library(scater)
library(dplyr)
library(readxl)
library(tidyr)
library(magrittr)
library(Rcpp)
library(coda)
library(MCMCglmm)
library(viridis)
library(ggmcmc)
library(ggplot2)
library(preprocessCore)

raw <- read_excel("data/nbt.3154-S3.xlsx")

gex <- as.matrix(raw[,-1])

cell_names <- raw$Cell
cell_type <- sapply(strsplit(cell_names, "_"), `[`, 1)

cts <- c("4SFG", "4SG", "HF", "NP", "PS")

for(ct in cts) cell_type[grep(ct, cell_type)] <- ct

pd <- AnnotatedDataFrame(data.frame(cell_name = cell_names, cell_type = cell_type))
rownames(pd) <- rownames(gex) <- cell_names

sce <- newSCESet(exprsDat = t(gex), phenoData = pd)
is_exprs(sce) <- exprs(sce) > 20
sce <- calculateQCMetrics(sce)


plotPCA(sce, colour_by = "cell_type", ncomponents = 3)
plotPCA(sce, colour_by = "pct_dropout", ncomponents = 3)



## tidy data analysis

raw_tidy <- mutate(raw, cell_type) %>% 
  gather(Gene, Expression, -cell_type, -Cell)

ggplot(raw_tidy, aes(x = cell_type, y = Expression)) + geom_violin()

ggplot(raw_tidy, aes(x = Gene, y = Expression)) + geom_boxplot()


saturation_gene <- group_by(raw_tidy, Gene) %>% 
  summarise(median_expression = median(Expression))

non_saturated_genes <- filter(saturation_gene, median_expression < 25) %>% extract2("Gene")

plotPCA(sce[non_saturated_genes, ], colour_by = "cell_type", ncomponents = 3)
plotPCA(sce[non_saturated_genes, ], colour_by = "pct_dropout", ncomponents = 3)

sce <- sce[non_saturated_genes, ]
## Let's check the distribution of cells

by_cell <- filter(raw_tidy, Gene %in% featureNames(sce)) %>% 
  group_by(cell_type, Cell) %>% 
  summarise(pct_60 = quantile(Expression, 0.6))

ggplot(by_cell, aes(x = cell_type, y = pct_60)) + geom_violin()

non_saturated_cells <- filter(by_cell, pct_60 < 22.5) %>% extract2("Cell")


plotPCA(sce[, non_saturated_cells], colour_by = "cell_type", ncomponents = 3)
plotPCA(sce[, non_saturated_cells], colour_by = "pct_dropout", ncomponents = 3)

ns_cell_type = as.character(sce[, non_saturated_cells]$cell_type)
ns_cell_names <- sampleNames(sce[, non_saturated_cells])

X <- exprs(sce[, non_saturated_cells])
X <- scale(t(X))


## Intensity density plots

set.seed(123L)
X %>%  as_data_frame() %>% mutate(ns_cell_names) %>%  sample_n(10) %>% 
  gather(Gene, Expression, -ns_cell_names) %>% 
  ggplot(aes(x = Expression, color = ns_cell_names)) + geom_density()

dd <- X %>%  as_data_frame() %>% mutate(ns_cell_type) %>% #  sample_n(200) %>% 
  gather(Gene, Expression, -ns_cell_type) 

ggplot(dd, aes(x = Expression, color = ns_cell_type)) + geom_density()
ggplot(dd, aes(x = Expression, color = Gene)) + geom_density()
ggplot(dd, aes(x = Gene, y = Expression)) + geom_boxplot()



Y <- t(normalize.quantiles(t(X)))
Y %>%  as_data_frame() %>% mutate(ns_cell_names) %>%  sample_n(10) %>% 
  gather(Gene, Expression, -ns_cell_names) %>% 
  ggplot(aes(x = Expression, color = ns_cell_names)) + geom_density()





sce2 <- sce[, non_saturated_cells]
exprs(sce2) <- t(Y)
plotPCA(sce2, colour_by = "cell_type", ncomponents = 3)
plotPCA(sce2, colour_by = "pct_dropout", ncomponents = 3)

source("houdini/R/houdini.R")
sourceCpp("houdini/src/gibbs.cpp")

xpca <- prcomp(Y)$x[,1:3] %>% as_data_frame() %>% mutate(ns_cell_type) 

xpca %>% 
  ggplot(aes(x = PC3, y = PC1, color = ns_cell_type)) + geom_point() +
  scale_colour_brewer(palette = "Set1")
xpca %>% 
  ggplot(aes(x = PC2, y = PC1, color = ns_cell_type)) + geom_point() +
  scale_colour_brewer(palette = "Set1")


g <- mfa_cpp(t(Y), iter = 60000, thin = 30, collapse = FALSE,
             pc_initialise = 3, eta_tilde = 0)


mc <- to_ggmcmc(g)

ggs_traceplot(filter(mc, Parameter == "lp__")) + stat_smooth()
ggs_autocorrelation(filter(mc, Parameter == "lp__"))

tmap <- posterior.mode(mcmc(g$pst_trace))
gamma_mean <- colMeans(g$gamma_trace)

xpca$tmap <- tmap
xpca$gamma_mean <- gamma_mean

ggplot(xpca, aes(x = PC3, y = PC1, color = ns_cell_type)) + geom_point() + scale_color_brewer(palette = "Set1")
ggplot(xpca, aes(x = PC2, y = PC1, color = ns_cell_type)) + geom_point()
ggplot(xpca, aes(x = PC3, y = PC1, color = tmap)) + geom_point() + scale_color_viridis()
ggplot(xpca, aes(x = PC2, y = PC1, color = tmap)) + geom_point() + scale_color_viridis()
ggplot(xpca, aes(x = PC3, y = PC1, color = gamma_mean)) + geom_point() + scale_color_viridis()

ggplot(xpca, aes(x = ns_cell_type, y = tmap)) + geom_boxplot() 


# Gene expression along pseudotime ---------------------------------------

dy <- Y %>% as_data_frame()
names(dy) <- featureNames(sce)
dy %<>% mutate(tmap, cell_type = ns_cell_type, gamma_mean) %>% 
  gather(gene, expression, -tmap, -cell_type, -gamma_mean)

ggplot(dy, aes(x = tmap, y = expression, color = round(gamma_mean))) + geom_point() +
  facet_wrap(~ gene)

ggplot(dy, aes(x = tmap, y = expression, color = cell_type)) + geom_point() +
  facet_wrap(~ gene) + scale_color_brewer(palette = "Set1")
