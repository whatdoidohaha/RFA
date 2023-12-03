
# This is the fast, approximate implementation of reference-free analysis that
# is tailored for binary genotype spaces (combining the worflows RFA-fast and
# RFA-binary).
#
# Written in R version 4.3.1.
# Requires R package glmnet (4.1-7) and lbfgs (1.2.1.2).


# Settings (recommended not to change in a run) ---------------------------

# Set working directory.
setwd('/Users/yeonwoopark/Desktop/RFA/Tutorials/')

# Source the script RFA-fast-binary-functions.R to import function definitions.
source('../Scripts/RFA-fast-binary-functions.R')

# Import genotype and phenotype data.
import_data('Data/avGFP/Genotypes.txt', 'Data/avGFP/Phenotypes.txt')
# Arguments:
#   Path to a csv file listing the genotypes (no header)
#   Path to a text file listing the phenotypes
# The following global variables are created:
#   genotypes : (numeric matrix) Genotype matrix
#   y : (numeric vector) Phenotype vector
#   n : (int) Number of sites
#   q_list : (int vector) Number of states in each site

# Choose a nonspecific epistasis model; two options are currently supported:
# 'none' and 'sigmoid'.
nonspec_type <- 'sigmoid'


# Defining the model ------------------------------------------------------

# When nonspecific epistasis is modeled, the parameters must be specified. For
# the sigmoid, this is the lower bound (L) and range (R = U - L). A reasonable
# initial value is the lower bound and range of the observed phenotype.
nonspec_param <- c(min(y), diff(range(y)))

# To include reference-free effects for all site-combinations up to some order:
model_ord <- 3
site_combn <- generate_site_combn(model_ord)

# Alternatively, to only include specific site-combinations:
# site_combn <- import_site_combn('Data/avGFP/Site-combinations.txt')
# Argument: Path to a text file listing the site-combinations.


# Inference without regularization ----------------------------------------

model <- infer_model(genotypes, y, site_combn)
# The following arguments are optional:
#   lambda : Lasso penalty (default 0, no regularization)
#
# A list with the following elements is returned:
#   nonspec_param : The nonspecific epistasis parameters used
#   e : The inferred reference-free effects
#   R2 : Model fit
#   residual : Residual for each genotype

# Model-predicted phenotypes.
p <- y - model$residual

# To obtain the genetic scores, first construct the phenotype operator,
# which is a mapping from reference-free effects to genetic score.
X <- construct_phenotype_operator(genotypes, site_combn)

# Then perform a matrix multiplication:
s <- as.vector(X %*% model$e)

# Genetic score can be transformed into predicted phenotype/
p <- apply_nonspecific_epistasis(s, model$nonspec_param)

# Check for any systematic misfit.
plot(p, y, cex = 0.5, col = rgb(0, 0, 0, 0.5), pch = 16L)

# Optimize the nonspecific epistasis parameters while fixing the reference-free
# effects. The global variable 'nonspec_param' is updated.
refine_nonspecific_epistasis(genotypes, y, site_combn, model$e, update = TRUE)

# Reinfer the reference-free effects based on the optimized sigmoid.
model <- infer_model(genotypes, y, site_combn)

# This stepwise optimization can be continued until the model fit does not
# increase any more. Update is terminated when R2 does not improve by at lesat
# 'epsilon'.
model <- sequential_update(genotypes, y, site_combn, model, epsilon = 0.001)

# Plot the distributon of genetic score and the sigmoid.
show_latent_space(model, genotypes, y, site_combn)

# Fraction of phenotypic variance attributable to each site-combination. The
# fraction attributable to each effect is this number divided by the number of
# effects in the site-combination.
v_site_combn <- model$e[-1] ^ 2 / sum(model$e[-1] ^ 2)

# The model can be exported.
export_model('model.txt', n, q_list, site_combn, model)

# The exported model can be imported:
# model <- import_model('model.txt')
# This creates the following global variables:
#   n, q_list, site_combn, nonspec_type, nonspec_param
# Make sure these global variables are compatiable with the genotype and
# phenotype data.


# Inference with regularization -------------------------------------------

# Number of genotype partitions for cross-validation.
k <- 5

# Number of times the k-fold cross-validation should be performed.
# Usually once is enough.
niter <- 1

# Performing cross-validation; calcuating the out-of-sample R2.
cvR2 <- cross_validation(genotypes, y, site_combn, k, niter)

# Choosing the penalty that maximizes the out-of-sample R2.
lambda <- cvR2[which.max(apply(cvR2[, -1, drop = FALSE], 1, mean)), 1]

# Diagnose overfit by comparing in- and out-of-sample fit.
diagnose_overfit(genotypes, y, site_combn, lambda, k)

# Cross-validated model.
model <- infer_model(genotypes, y, site_combn, lambda)


# Minimal models with maximal explanatory power ---------------------------

# Perform cross-validation to evaluate a series of nested models constructed by
# including reference-free effects in the order of their phenotypic contribution.

# Maximum number of effects to include.
max_size <- 25

# Cross-validation setting.
k <- 5 # Number of genotype partitions

# Construct the minimal models and evaluate them by cross-validation.
minR2 <- minimal_models(genotypes, y, site_combn, model, max_size, k)
# The argument `model` is the model to use for computing the phenotypic
# contribution of each effect, which determines the order of their inclusion.
