
source('../../../scripts/ref_free.R')

library(ggplot2)

set_model(list('../../../data/Simulated/10_4_only_3rd_nonspec/genotypes.txt',
               '../../../data/Simulated/10_4_only_3rd_nonspec/phenotypes.txt', '1', 2))

s <- unique(q_list) # Number of states must be identical across sites.

# True effects.
e_true <- as.numeric(readLines('../../../data/simulated/10_4_only_3rd_nonspec/e.txt'))


# Complete sampling -------------------------------------------------------

lambda_list <- c(0, 10 ^ seq(-5, 1, by = 1))
k <- 10L # Five-fold cross-validation

# First-order.

cv_F1 <- cross_validation(genotypes, y, generate_site_combn(1L), k, lambda_list)
lambda_F1 <- identify_lambda(lambda_list, cv_F1, show_plot = TRUE)
F1 <- infer_model(genotypes, y, generate_site_combn(1L), lambda_F1, norm_e = TRUE)
colnames(cv_F1) <- lambda_list
write.csv(round(cv_F1, 3L), 'cv_F1.txt', row.names = FALSE)

# Second-order.

cv_F2 <- cross_validation(genotypes, y, generate_site_combn(2L), k, lambda_list)
lambda_F2 <- identify_lambda(lambda_list, cv_F2, show_plot = TRUE)
F2 <- infer_model(genotypes, y, generate_site_combn(2L), lambda_F2, norm_e = TRUE)
colnames(cv_F2) <- lambda_list
write.csv(round(cv_F2, 3L), 'cv_F2.txt', row.names = FALSE)

# Third-order.

cv_F3 <- cross_validation(genotypes, y, generate_site_combn(3L), k, lambda_list)
lambda_F3 <- identify_lambda(lambda_list, cv_F3, show_plot = TRUE)
F3 <- infer_model(genotypes, y, generate_site_combn(3L), lambda_F3, norm_e = TRUE)
colnames(cv_F3) <- lambda_list
write.csv(round(cv_F3, 3L), 'cv_F3.txt', row.names = FALSE)

# Partial sampling --------------------------------------------------------

set.seed(10L)
subset <- sample(nrow(genotypes), size = nrow(genotypes) * 0.4)

# First-order.

cv_P1 <- cross_validation(genotypes[subset, ], y[subset], generate_site_combn(1L), k, lambda_list)
lambda_P1 <- identify_lambda(lambda_list, cv_P1, show_plot = TRUE)
P1 <- infer_model(genotypes[subset, ], y[subset], generate_site_combn(1L), lambda_P1, norm_e = TRUE)
colnames(cv_P1) <- lambda_list
write.csv(round(cv_P1, 3L), 'cv_P1.txt', row.names = FALSE)

# Second-order.

cv_P2 <- cross_validation(genotypes[subset, ], y[subset], generate_site_combn(2L), k, lambda_list)
lambda_P2 <- identify_lambda(lambda_list, cv_P2, show_plot = TRUE)
P2 <- infer_model(genotypes[subset, ], y[subset], generate_site_combn(2L), lambda_P2, norm_e = TRUE)
colnames(cv_P2) <- lambda_list
write.csv(round(cv_P2, 3L), 'cv_P2.txt', row.names = FALSE)

# Third-order.

cv_P3 <- cross_validation(genotypes[subset, ], y[subset], generate_site_combn(3L), 10, lambda_list)
lambda_P3 <- identify_lambda(lambda_list, cv_P3, show_plot = TRUE)
P3 <- infer_model(genotypes[subset, ], y[subset], generate_site_combn(3L), lambda_P3, norm_e = TRUE)
colnames(cv_P3) <- lambda_list
write.csv(round(cv_P3, 3L), 'cv_P3.txt', row.names = FALSE)


# Plotting ----------------------------------------------------------------

# Genetic score.

s <- -log(1 / y - 1)
h <- hist(s, breaks = seq(-10, 6, by = 1), plot = FALSE)
h$counts <- h$counts / sum(h$counts) * 100
plot(h, xlim = c(-10, 5), ylim = c(0, 20))
points(seq(-10, 5, by = 0.01), 20 / (1 + exp(-seq(-10, 5, by = 0.01))), type = 'l', col = 'red')

# Model fit.

model <- P3
model_ord <- 3L
ord <- unlist(lapply(1L:model_ord, function(ord) rep(ord, choose(n, ord) * s ^ ord)))

data <- data.frame(e = model$e[-1L], ord = as.factor(ord))

ggplot(data[ord == 1L, ], aes(x = ord, y = e)) +
  geom_violin() +
  ylim(-4.5, 4.5)


