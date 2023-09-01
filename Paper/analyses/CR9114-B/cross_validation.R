
# Reference-free analysis modeling nonspecific epistasis; cross-validation.
# Cross-validation is performed across measurement replicates.
# Binary state space; using the script ref_free_binary.R.


# Data import -------------------------------------------------------------

source('../../scripts/ref_free_binary.R')
dataset <- 'CR9114-B'

# Control file:
#   1. Directory to the genotype matrix.
#   2. Directory to the phenotype vector.
#   3. Model order (or directory to a site-combination list).
#   4. Type of nonspecific epistasis.
control <- list(paste0('../../data/', dataset, '/genotypes.txt'), paste0('../../data/', dataset, '/phenotypes.txt'), '1', 2)
set_model(control)

# Measurement replicates.
replicates <- as.matrix(read.csv(paste0('../../data/', dataset, '/phenotype_replicates.txt'), header = TRUE))


# Run ---------------------------------------------------------------------

# First-order.

lambda_list <- 10 ^ seq(-3, 1, by = 0.5) # Regularization strengths.

# Optimization in this dataset is prone to local optima (because of strong lower bounding) and requires a close initial value.
init_model <- infer_model(genotypes, y, generate_site_combn(1L))

cv_R2 <- cross_validation_across_replicates(genotypes, replicates, generate_site_combn(1L), lambda_list,
                                            init_model = init_model)
colnames(cv_R2) <- log(lambda_list, 10)

write.csv(round(cv_R2, 3L), 'cv_R2_1st.txt', row.names = FALSE)

# Visualization.
# dotplot_with_sd(log(lambda_list, 10), t(cv_R2), plot_SEM = TRUE, ylim = c(0, 1))


# Second-order.

lambda_list <- 10 ^ seq(-3, 1, by = 0.5) # Regularization strengths.

init_model$e <- c(init_model$e, rep(0, length(generate_site_combn(2L)) - length(init_model$e)))

cv_R2 <- cross_validation_across_replicates(genotypes, replicates, generate_site_combn(2L), lambda_list,
                                            init_model = init_model)
colnames(cv_R2) <- log(lambda_list, 10)

write.csv(round(cv_R2, 3L), 'cv_R2_2nd.txt', row.names = FALSE)


# Third-order.

lambda_list <- 10 ^ seq(-3, 1, by = 1) # Regularization strengths.

init_model$e <- c(init_model$e, rep(0, length(generate_site_combn(3L)) - length(init_model$e)))

cv_R2 <- cross_validation_across_replicates(genotypes, replicates, generate_site_combn(3L), lambda_list,
                                            init_model = init_model)
colnames(cv_R2) <- log(lambda_list, 10)

write.csv(round(cv_R2, 3L), 'cv_R2_3rd.txt', row.names = FALSE)
