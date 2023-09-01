
# Reference-free analysis modeling nonspecific epistasis; sparsity analysis.
# Binary state space; using the script ref_free_binary.R.


# Data import -------------------------------------------------------------

source('../../scripts/ref_free_binary.R')
dataset <- 'CR9114-H1'

# Control file:
#   1. Directory to the genotype matrix.
#   2. Directory to the phenotype vector.
#   3. Model order (or directory to a site-combination list).
#   4. Type of nonspecific epistasis.
control <- list(paste0('../../data/', dataset, '/genotypes.txt'), paste0('../../data/', dataset, '/phenotypes.txt'), '1', 2)
set_model(control)

# Model to determine the order of parameter expansion.
model <- import_model('regularized_fit_3rd.txt')


# Setting and run ---------------------------------------------------------

# Cross-validation is performed to calculate the out-of-sample R2 of each model.
k <- 5L
lambda_list <- 10 ^ seq(-3, 2, by = 1)

analyze_sparsity(genotypes, y, model, site_combn, k, lambda_list)
# Output is written in a .txt file named sparsity.

