
# Reference-free analysis modeling nonspecific epistasis; cross-validation.
# Binary state space; using the script ref_free_binary.R.


# Data import -------------------------------------------------------------

source('../../scripts/ref_free_binary.R')
dataset <- 'CR6261-H9'

# Control file:
#   1. Directory to the genotype matrix.
#   2. Directory to the phenotype vector.
#   3. Model order (or directory to a site-combination list).
#   4. Type of nonspecific epistasis.
control <- list(paste0('../../data/', dataset, '/genotypes.txt'), paste0('../../data/', dataset, '/phenotypes.txt'), '1', 2)
set_model(control)


# Run ---------------------------------------------------------------------

# First-order.

lambda_list <- 10 ^ seq(-3, 1, by = 0.5) # Regularization strengths.

# Ten-fold cross-validation repeated three times.
cv_cor <- cross_validation(genotypes, y, generate_site_combn(1L), k = 10L, lambda_list, niter = 3L, epsilon = 5e-4)
colnames(cv_cor) <- log(lambda_list, 10)

write.csv(signif(cv_cor, 3L), 'cv_R2_1st.txt', row.names = FALSE)

# Visualization.
# dotplot_with_sd(log(lambda_list, 10), t(cv_cor), plot_SEM = TRUE, ylim = c(0, 1))


# Second-order.

lambda_list <- 10 ^ seq(-3, 1, by = 0.5) # Regularization strengths.

# Ten-fold cross-validation repeated three times.
cv_cor <- cross_validation(genotypes, y, generate_site_combn(2L), k = 10L, lambda_list, niter = 3L, epsilon = 5e-4)
colnames(cv_cor) <- log(lambda_list, 10)

write.csv(signif(cv_cor, 3L), 'cv_R2_2nd.txt', row.names = FALSE)


# Third-order.

lambda_list <- 10 ^ seq(-3, 1, by = 1) # Regularization strengths.

# Ten-fold cross-validation, repeated three times.
cv_cor <- cross_validation(genotypes, y, generate_site_combn(3L), k = 10L, lambda_list, niter = 3L, epsilon = 5e-4)
colnames(cv_cor) <- log(lambda_list, 10)

write.csv(signif(cv_cor, 3L), 'cv_R2_3rd.txt', row.names = FALSE)

