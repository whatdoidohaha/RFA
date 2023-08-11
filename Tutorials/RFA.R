
# Here we illustrate reference-free analysis on a genetic architecture that
# consists of four sites and a variable number of states per site. We begin with
# a simple regression that jointly infers reference-free effects and nonspecific
# epistasis. We then proceed to regularized regression, performing cross-
# validation to identify the optimal lasso penalty. Finally, we build minimal
# reference-free models with maximal explanatory power, which reveal the
# sparsity of genetic architecture.
#
# Jointly inferring reference-free effects and nonspecific epistasis requires
# regularized nonlinear regression, which is implemented here using a generic
# numerical optimizer (R package lbfgs). This makes the script slow; when the
# number of genotypes exceeds 100,000, cross-validating a second- or third-order
# model may take days. In addition, we find that the optimizer may get stuck in
# local optima when the fraction of genotypes in the dynamic range of
# measurement is very small. An approximate optimization that is faster and more
# robust yet largely retains the accuracy is presented in the workflow RFA-fast.
#
# Written in R version 4.3.1.
# Requires R packages abind (1.4-5), lbfgs (1.2.1.2), and Matrix (1.6-0).


# Settings (recommended not to change in a run) ---------------------------

# Set working directory.
setwd('/Users/yeonwoopark/Desktop/RFA/Tutorials/')

# Source the script RFA-functions.R to import function definitions.
source('../Scripts/RFA-functions.R')

# Import genotype and phenotype data.
import_data('Data/ParD3-ParE3/Genotypes.txt', 'Data/ParD3-ParE3/Phenotypes.txt')
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


# Defining the reference-free model ---------------------------------------

# To include effects for all site-combinations up to a given order:
model_ord <- 1
site_combn <- generate_site_combn(model_ord)

# Alternatively, to only include specific site-combinations:
# site_combn <- import_site_combn('Data/ParD3-ParE3/Site-combinations.txt')
# Argument: Path to a text file listing the site-combinations.


# Inference without regularization ----------------------------------------

model <- infer_model(genotypes, y, site_combn)
# The following arguments are optional:
#   lambda : Lasso penalty (default 0, no regularization)
#   init_model : Model to use as initial value for optimization
#   norm_e : Whether to enforce the zero-mean property on the returned effects
#     (default TRUE)
#   silent : Show (0) or suppress (1) lbfgs::lbfgs console output (default 0)
#   epsilon : Precision; see the documentation for lbfgs::lbfgs (default 1e-5)
#
# A list with the following elements is returned:
#   nonspec_param : The inferred nonspecific epistasis parameters; for the
#      sigmoid, this is the lower bound (L) and the range (U - L)
#   e : The inferred reference-free effects
#   R2 : Model fit
#   residual : Residual for each genotype

# Inferred effects; effects in each site-combination are arranged in an array.
e_site_combn <- parse_effects_by_site_combn(site_combn, model$e, TRUE)

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

# Fraction of phenotypic variance attributable to each effect.
v <- partition_variance(site_combn, model$e)

# Fraction of phenotypic variance attributable to each site-combination.
v_site_combn <- sapply(v, sum)

# The model can be exported.
export_model('model.txt', n, q_list, site_combn, model)

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
k <- 3

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
