
library(ggplot2)
library(viridis)
library(dplyr)
library(RColorBrewer)
library(magrittr)
library(tidyr)


load("data/wishbone.Rdata")

mean(exprs(sce) == 0)

make_dropout_plot <- function(sce, title) {
  df <- data_frame(
    mean_expr = rowMeans(exprs(sce)),
    prop_drop = rowMeans(exprs(sce) == 0),
    var_expr = matrixStats::rowVars(exprs(sce)))
  
  fit <- nls(prop_drop ~ exp(-lambda * mean_expr), data = df)
  fit_d <- nls(prop_drop ~ exp(-lambda * mean_expr^2), data = df)
  
  f_dropout <- function(x) exp(-coef(fit)[1] * x)
  fd_dropout <- function(x) exp(-coef(fit)[1] * x^2)
  
  
  df %<>% mutate(exp = f_dropout(mean_expr),
                 dexp = fd_dropout(mean_expr))
  
  dfm <- gather(df, model, value, -mean_expr, -prop_drop, -var_expr)
  
  ggplot(arrange(dfm, mean_expr), aes(x = mean_expr, y = prop_drop)) +
    geom_point(shape = 21, alpha = 0.5) + scale_fill_viridis() +
    geom_line(aes(y = value, color = model)) + ggtitle(title) +
    scale_color_brewer(palette = "Set1")
}

paul_plot <- make_dropout_plot(sce, "Paul")

library(HSMMSingleCell)
library(monocle)
hsmm <- HSMM; rm(HSMM)
exprs(hsmm) <- log2(exprs(hsmm) + 1)
trapnell_plot <- make_dropout_plot(hsmm, "Trapnell")

cowplot::plot_grid(paul_plot, trapnell_plot)


## experimental


