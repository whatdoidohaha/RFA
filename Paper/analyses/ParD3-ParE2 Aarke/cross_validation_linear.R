
# Reference-free analysis without modeling nonspecific epistasis; cross-validation.
# Non-binary state space; using the script ref_free.R.


# Data import -------------------------------------------------------------

source('../../scripts/ref_free.R')
dataset <- 'ParD3-ParE2 Aarke'

# Control file:
#   1. Directory to the genotype matrix.
#   2. Directory to the phenotype vector.
#   3. Model order (or directory to a site-combination list).
#   4. Type of nonspecific epistasis.
control <- list(paste0('../../data/', dataset, '/genotypes.txt'), paste0('../../data/', dataset, '/phenotypes.txt'), '1', 1)
set_model(control)


# Run ---------------------------------------------------------------------

# First-order.

lambda_list <- 10 ^ seq(-3, 1, by = 0.5) # Regularization strengths.

# Ten-fold cross-validation repeated three times.
cv_R2 <- cross_validation(genotypes, y, generate_site_combn(1L), k = 10L, lambda_list, niter = 3L, epsilon = 5e-4)
colnames(cv_R2) <- log(lambda_list, 10)

write.csv(round(cv_R2, 3L), 'cv_R2_1st_linear.txt', row.names = FALSE)

# Visualization.
# dotplot_with_sd(log(lambda_list, 10), t(cv_R2), plot_SEM = TRUE, ylim = c(0, 1))


# Second-order.

lambda_list <- 10 ^ seq(-3, 1, by = 0.5) # Regularization strengths.

# Ten-fold cross-validation repeated three times.
cv_R2 <- cross_validation(genotypes, y, generate_site_combn(2L), k = 10L, lambda_list, niter = 3L, epsilon = 5e-4)
colnames(cv_R2) <- log(lambda_list, 10)

write.csv(round(cv_R2, 3L), 'cv_R2_2nd_linear.txt', row.names = FALSE)


# Third-order.

lambda_list <- 10 ^ seq(-3, 1, by = 1) # Regularization strengths.

# Ten-fold cross-validation, repeated three times.
cv_R2 <- cross_validation(genotypes, y, generate_site_combn(3L), k = 10L, lambda_list, niter = 3L, epsilon = 5e-4)
colnames(cv_R2) <- log(lambda_list, 10)

write.csv(round(cv_R2, 3L), 'cv_R2_3rd_linear.txt', row.names = FALSE)
