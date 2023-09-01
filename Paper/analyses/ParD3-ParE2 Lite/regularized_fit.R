
# Reference-free analysis modeling nonspecific epistasis; regularized model.
# Non-binary state space; using the script ref_free.R.


# Data import -------------------------------------------------------------

source('../../scripts/ref_free.R')
dataset <- 'ParD3-ParE2 Lite'

# Control file:
#   1. Directory to the genotype matrix.
#   2. Directory to the phenotype vector.
#   3. Model order (or directory to a site-combination list).
#   4. Type of nonspecific epistasis.
control <- list(paste0('../../data/', dataset, '/genotypes.txt'), paste0('../../data/', dataset, '/phenotypes.txt'), '1', 2)
set_model(control)


# Run ---------------------------------------------------------------------

# Best-fit first-order model.

lambda <- 10 ^ identify_lambda(dir = 'cv_R2_1st.txt', show_plot = TRUE) # Best regularization strength
site_combn <- generate_site_combn(1L)
model <- infer_model(genotypes, y, site_combn, lambda)
export_model('regularized_fit_1st.txt', n, q_list, site_combn, nonspec_type, model$nonspec_param, model$e)


# Best-fit second-order model.

lambda <- 10 ^ identify_lambda(dir = 'cv_R2_2nd.txt', show_plot = TRUE) # Best regularization strength
site_combn <- generate_site_combn(2L)
model <- infer_model(genotypes, y, site_combn, lambda)
export_model('regularized_fit_2nd.txt', n, q_list, site_combn, nonspec_type, model$nonspec_param, model$e)


# Best-fit third-order model.
# Inapplicable because there are only three sites.
