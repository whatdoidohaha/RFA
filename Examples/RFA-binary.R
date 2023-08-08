
# Here we illustrate reference-free analysis of a binary genotype space. When
# there are only two states per site, the zero-mean property causes all effects
# in a site-combination to equal in magnitude and vary only in sign. Therefore,
# only one effect needs to be explicitly modeled and inferred for each site-
# combination. (Simplex encoding and graph Fourier transform exploit the zero-
# mean property for multiple states to always achieve this compact encoding.)
# This compact encoding causes the phenotype of any genotype to depend on every
# model parameter and therefore the inference to be more sensitive to missing
# genotypes (see Fig. 1E in our paper). However, the impact on sensitivity is
# minor when the number of states is small and is outweighed by the benefit of
# compactness.
#
# For each site-combination, this script only models the effect corresponding to
# the second state in every site. As in the example 'RFA', reference-free
# effects and nonspecific epistasis are jointly inferred using an exact but slow
# algorithm. For a fast, approximate inference on binary genotype spaces, check
# the example 'RFA-binary-fast'.
#
# Written in R version 4.3.1.
# Requires R packages lbfgs (1.2.1.2).


# Settings (recommended not to change in a run) ---------------------------

# Set working directory.
setwd('/Users/yeonwoopark/Desktop/RFA/Examples/')

# Source the script RFA-binary-functions.R to import function definitions.
source('../Scripts/RFA-binary-functions.R')

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
# none (1) or sigmoid (2).
nonspec_type <- 2


# Defining the reference-free model ---------------------------------------

# To include effects for all site-combinations up to a given order:
model_ord <- 1
site_combn <- generate_site_combn(model_ord)

# Alternatively, to only include specific site-combinations:
# site_combn <- import_site_combn('Data/avGFP/Site-combinations.txt')
# Argument: Path to a text file listing the site-combinations.


# Inference without regularization ----------------------------------------

model <- infer_model(genotypes, y, site_combn)
# The following arguments are optional:
#   lambda : Lasso penalty (default 0, no regularization)
#   init_model : Model to use as initial value for optimization
#   silent : Show (0) or suppress (1) lbfgs::lbfgs console output (default 0)
#   epsilon : Precision; see the documentation for lbfgs::lbfgs (default 1e-5)
#
# A list with the following elements is returned:
#   nonspec_param : The inferred nonspecific epistasis parameters; for the
#      sigmoid, this is the lower bound (L) and the range (U - L)
#   e : The inferred reference-free effects
#   R2 : Model fit
#   residual : Residual for each genotype

# Model-predicted phenotypes.
p <- y - model$residual

# To obtain the genetic score, first construct the phenotype operator, which is
# a matrix mapping reference-free effects to genetic score.
X <- construct_phenotype_operator(genotypes, site_combn)

# Then perform a matrix multiplication:
s <- as.vector(X %*% model$e)

# Genetic score can be transformed into predicted phenotype by applying the
# inferred nonspecific epistasis.
p <- apply_nonspecific_epistasis(s, model$nonspec_param)

# Plot the distributon of genetic score and the inferred nonspecific epistasis.
show_latent_space(model, genotypes, y, site_combn)

# Fraction of phenotypic variance attributable to each site-combination. The
# fraction attributable to each effect is this number divided by the number of
# effects in the site-combination.
v_site_combn <- model$e[-1] ^ 2 / sum(model$e[-1] ^ 2)

# The model can be exported.
export_model('model.txt', n, q_list, site_combn, nonspec_type, model)

# The exported model can be imported:
# model <- import_model('model.txt')
# This creates the following global variables:
#   n, q_list, site_combn, nonspec_type
# Make sure these global variables are compatiable with the genotype and
# phenotype data.


# Inference with regularization -------------------------------------------

# Lasso penalties to evaluate by cross-validation.
lambda_list <- c(0, 1e-2, 1e-1, 1, 10)

# Number of genotype partitions for cross-validation.
k <- 5

# Number of times the k-fold cross-validation should be performed.
# Usually once is enough.
niter <- 1

# Performing cross-validation; calcuating the out-of-sample R2.
cvR2 <- cross_validation(genotypes, y, site_combn, k, lambda_list, niter)

# Choosing the penalty that maximizes the out-of-sample R2.
lambda <- lambda_list[which.max(apply(cvR2, 2, mean))]

# Cross-validated model.
model <- infer_model(genotypes, y, site_combn, lambda)

# When performing cross-validation across genotype partitions is not possible
# because there are too few genotypes or too few are in the dynamic range of
# measurement, cross-validation can be performed across measurement replicates.
# Import the measurement replicates into a matrix (`replicates`) where the rows
# correspond to genotypes and columns to replicates. Then run:
# cvR2 <- cv_across_replicates(genotypes, replicates, site_combn, lambda_list)


# Minimal models with maximal explanatory power ---------------------------

# Perform cross-validation to evaluate a series of nested models constructed by
# including reference-free effects in the order of their phenotypic contribution.

# Maximum number of effects to include.
max_size <- 25

# Cross-validation setting.
k <- 3 # Number of genotype partitions
lambda_list <- c(1e-2, 1e-1, 1, 10) # Lasso penalties; do not include 0

# Construct the minimal models and evaluate them by cross-validation.
minR2 <- minimal_models(genotypes, y, site_combn, model, max_size, k, lambda_list)
# The argument `model` is the model to use for computing the phenotypic
# contribution of each effect, which determines the order of their inclusion.

# As before, if cross-validation cannot be performed across genotype partitions,
# it can be performed across measurement replicates.
# minR2 <- minimal_models_replicates(genotypes, replicates, site_combn, model,
#                                    max_size, lambda_list)
