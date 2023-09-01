
# Reference-free analysis modeling nonspecific epistasis; prediction from subsample.
# Binary state space; using the script ref_free_binary.R.


# Data import -------------------------------------------------------------

source('../../scripts/ref_free_binary.R')
dataset <- 'MPH'

# Control file:
#   1. Directory to the genotype matrix.
#   2. Directory to the phenotype vector.
#   3. Model order (or directory to a site-combination list).
#   4. Type of nonspecific epistasis.
control <- list(paste0('../../data/', dataset, '/genotypes.txt'), paste0('../../data/', dataset, '/phenotypes.txt'), '1', 2)
set_model(control)

# Phenotype replicates.
replicates <- as.matrix(read.csv(paste0('../../data/', dataset, '/phenotype_replicates.txt'), header = TRUE))


# Setting -----------------------------------------------------------------

# Model orders to test.
# The model that best describes the subsample is used for predicting the phenotypes of all other genotypes.
order_list <- c(1L, 2L)

# Sample sizes to test (fraction of all possible genotypes, not that of experimentally characterized genotypes).
frac_list <- seq(0.5, 0.9, by = 0.1)

# Number of samples to try for each sample size.
niter <- 25L

# Regularization strengths.
lambda_list <- 10 ^ seq(-3, 2, by = 1)


# Run ---------------------------------------------------------------------

prediction <- predict_from_subsample_replicates(genotypes, y, replicates, order_list, frac_list, niter, lambda_list)

write.csv(round(t(prediction), 3L), 'prediction.txt', row.names = FALSE)
