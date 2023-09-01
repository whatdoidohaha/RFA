
# Reference-free analysis modeling nonspecific epistasis; sparsity analysis.
# Cross-validation is performed across replicate measurements.
# Non-binary state space; using the script ref_free.R.


# Data import -------------------------------------------------------------

source('../../scripts/ref_free.R')
dataset <- 'DHFR'

# Control file:
#   1. Directory to the genotype matrix.
#   2. Directory to the phenotype vector.
#   3. Model order (or directory to a site-combination list).
#   4. Type of nonspecific epistasis.
control <- list(paste0('../../data/', dataset, '/genotypes.txt'), paste0('../../data/', dataset, '/phenotypes.txt'), '1', 2)
set_model(control)

# Model to determine the order of parameter expansion.
model <- import_model('regularized_fit_2nd.txt')

# Phenotype replicates.
replicates <- as.matrix(read.csv(paste0('../../data/', dataset, '/phenotype_replicates.txt'), header = TRUE))


# Setting and run ---------------------------------------------------------

# Cross-validation is performed across measurement replicates.
lambda_list <- 10 ^ seq(-3, 1, by = 1)

analyze_sparsity_replicates(genotypes, replicates, model, site_combn, lambda_list)
# Output is written in a .txt file named sparsity.
