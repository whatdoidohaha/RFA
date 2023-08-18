
# Here we illustrate a fast, approximate inference of reference-free effects
# using an example genetic architecture consisting of four sites and a variable
# number of states per site. Instead of jointly inferring reference-free
# effects with nonspecific epistasis, this script asks the user to specify
# the nonspecific epistasis parameters. It then performs generalized linear
# regression using the specified function as link function. The highly
# optimized R package glmnet makes this inference fast and robust. Once the
# reference-free effects are inferred, the nonspecific epistasis parameters can
# be optimized while holding the effects fixed. This two-step procedure can be
# iterated to further improve the model fit.
#
# We recommend first performing unregularized regression in the worflow RFA
# to identify suitable nonspecific epistasis parameters. This may be possible
# using a first-order model or may require a second-order model.
#
# Written in R version 4.3.1.
# Requires R packages abind (1.4-5), glmnet (4.1-7), and lbfgs (1.2.1.2).


# Settings (recommended not to change in a run) ---------------------------

# Set working directory.
setwd('/Users/yeonwoopark/Desktop/RFA/Tutorials/')

# Source the script RFA-fast-functions.R to import function definitions.
source('../Scripts/RFA-fast-functions.R')

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


# Defining the model ------------------------------------------------------

# When nonspecific epistasis is modeled, the parameters must be specified. For
# the sigmoid, this is the lower bound (L) and range (R = U - L). A reasonable
# initial value is the lower bound and range of the observed phenotype.
nonspec_param <- c(min(y), diff(range(y)))

# To include reference-free effects for all site-combinations up to some order:
model_ord <- 2
site_combn <- generate_site_combn(model_ord)

# Alternatively, to only include specific site-combinations:
# site_combn <- import_site_combn('Data/ParD3-ParE3/Site-combinations.txt')
# Argument: Path to a text file listing the site-combinations.


# Inference without regularization ----------------------------------------

model <- infer_model(genotypes, y, site_combn)
# The following arguments are optional:
#   lambda : Lasso penalty (default 0, no regularization)
#   norm_e : Whether to enforce the zero-mean property on the returned 
#     effects (default TRUE)
#
# A list with the following elements is returned:
#   nonspec_param : Nonspecific epistasis parameters
#   e : Reference-free effects on genetic score
#   R2 : Model fit
#   residual : Residual for each genotype

# Arranging the effects for each site-combination in an array.
e_site_combn <- parse_effects_by_site_combn(site_combn, model$e, TRUE)

# Model-predicted phenotypes.
p <- y - model$residual

# To obtain the genetic scores, first construct the phenotype operator,
# which is a mapping from reference-free effects to genetic score.
X <- construct_phenotype_operator(genotypes, site_combn)

# Then perform a matrix multiplication.
s <- as.vector(X %*% model$e)

# Genetic score can be transformed into predicted phenotype.
p <- apply_nonspecific_epistasis(s, model$nonspec_param)

# Check for any systematic misfit.
plot(p, y, cex = 0.5, col = rgb(0, 0, 0, 0.5), pch = 16L)

# Optimize the nonspecific epistasis parameters while fixing the reference-free
# effects. The global variable 'nonspec_param' is updated.
refine_nonspecific_epistasis(genotypes, y, site_combn, model$e, update = TRUE)

# Reinfer the reference-free effects based on the optimized nonspecific
# epistasis parameters.
model <- infer_model(genotypes, y, site_combn)

# This stepwise optimization can be continued until the model fit does not
# increase any more. Update is terminated when R2 does not improve by at lesat
# 'epsilon'.
model <- sequential_update(genotypes, y, site_combn, model, epsilon = 0.001)

# Plot the distributon of genetic score and the nonspecific epistasis.
show_latent_space(model, genotypes, y, site_combn)

# Fraction of phenotypic variance due to each effect.
v <- partition_variance(site_combn, model$e)

# Fraction of phenotypic variance due to each site-combination.
v_site_combn <- sapply(v, sum)

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
max_size <- 50

# Cross-validation setting.
k <- 5 # Number of genotype partitions

# Construct the minimal models and evaluate them by cross-validation.
minR2 <- minimal_models(genotypes, y, site_combn, model, max_size, k)
# The argument `model` is the model to use for computing the phenotypic
# contribution of each effect, which determines the order of their inclusion.
