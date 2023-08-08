
# Here we illustrate a fast, approximate inference of reference-free effects
# using an example genetic architecture consisting of four sites and a variable
# number of states per site. Instead of jointly inferring nonspecific epistasis
# with the reference-free effects, this script asks the user to specify a
# sigmoid function for modeling nonspecific epistasis. It then performs a
# generalized linear regression using the sigmoid as link function. The highly
# optimized R package glmnet makes this inference fast and robust. Once the
# reference-free effects are inferred, the sigmoid parameters can be optimized
# with the effects fixed This two-step procedure can be iterated until the
# parameters stabilize.
#
# We recommend first performing unregularized regression to identify the optimal
# sigmoid parameters. This may be possible using a first-order reference-free
# model if it fits the data sufficiently well, or may require a second-order
# model. We also provide functions for performing cross-validation and building
# minimal reference-free models with maximal explanatory power.
#
# Written in R version 4.3.1.
# Requires R packages abind (1.4-5), glmnet (4.1-7), and lbfgs (1.2.1.2).


# Settings (recommended not to change in a run) ---------------------------

# Set working directory.
setwd('/Users/yeonwoopark/Desktop/RFA/Examples/')

# Source the script RFA-fast-functions.R to import function definitions.
source('../Scripts/RFA-fast-functions.R')

# Import genotype and phenotype data.
import_data('Data/ParD3-ParE3-Aarke/Genotypes.txt',
            'Data/ParD3-ParE3-Aarke/Phenotypes.txt')
# Arguments:
#   Path to a csv file listing the genotypes (no header)
#   Path to a text file listing the phenotypes
# The following global variables are created:
#   genotypes : (numeric matrix) Genotype matrix
#   y : (numeric vector) Phenotype vector
#   n : (int) Number of sites
#   q_list : (int vector) Number of states in each site


# Defining the model ------------------------------------------------------

# Specify the lower bound (L) and range (R) of the sigmoid model of nonspecific
# epistasis. Example: the lower bound and range of observed phenotype.
L <- min(y)
R <- diff(range(y))

# To include reference-free effects for all site-combinations up to some order:
model_ord <- 1
site_combn <- generate_site_combn(model_ord)

# Alternatively, to only include specific site-combinations:
# site_combn <- import_site_combn('Data/ParD3-ParE3-Aarke/Site-combinations.txt')
# Argument: Path to a text file listing the site-combinations.


# Inference without regularization ----------------------------------------

model <- infer_model(genotypes, y, site_combn)
# The following arguments are optional:
#   lambda : Lasso penalty (default 0, no regularization)
#   norm_e : Whether to enforce the zero-mean constraint on the returned 
#     effects (default TRUE)
#
# A list with the following elements is returned:
#   sigmoid : The sigmoid parameters used
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
p <- apply_sigmoid(s, L, R)

# Check for any systematic misfit.
plot(p, y, cex = 0.5, col = rgb(0, 0, 0, 0.5), pch = 16L)

# Optimize the sigmoid while fixing the reference-free effects.
refine_sigmoid(genotypes, y, site_combn, model$e, update = TRUE)

# Reinfer the reference-free effects based on the optimized sigmoid.
model <- infer_model(genotypes, y, site_combn)

# This stepwise optimization can be continued until the model fit does not
# increase any more.
model <- sequential_update(genotypes, y, site_combn, model, epsilon = 0.001)
# Update is terminated when R2 does not improve by at least epsilon.

# Plot the distributon of genetic score and the sigmoid.
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
#   n, q_list, site_combn, L, R
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
lambda <- cvR2[which.max(apply(cvR2[, -1], 1, mean)), 1]

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
