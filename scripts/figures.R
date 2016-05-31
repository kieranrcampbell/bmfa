library(cowplot)
library(ggplot2)


# Multiple trajectory figures ---------------------------------------------
ksigmoid <- function(t, phi, delta, k) {
  2 * phi / (1 + exp(-k * (t - delta)))
}

plt1 <- ggplot(data.frame(x = c(0, 1)), aes(x = x)) +
  stat_function(fun = ksigmoid, args = list(phi = 1, delta = 0.75, k = 20), 
                colour = "darkred", size = 1.5) + 
  stat_function(fun = ksigmoid, args = list(phi = 0, delta = 0.75, k = 0), 
                colour = "darkblue", size = 1.5) +
  xlab("Pseudotime") + ylab("Expression")

plt2 <- ggplot(data.frame(x = c(0, 1)), aes(x = x)) +
  stat_function(fun = ksigmoid, args = list(phi = 1, delta = 0.75, k = -20), 
                colour = "darkred", size = 1.5) + 
  stat_function(fun = ksigmoid, args = list(phi = 2, delta = 0.75, k = 0), 
                colour = "darkblue", size = 1.5) +
  xlab("Pseudotime") + ylab("Expression")

plot_grid(plt1, plt2, nrow = 1)
ggsave("figs/branch_expression.png", width = 6, height = 2.5)
