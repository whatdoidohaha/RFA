
# This script implements the functions for fast, approximate inference of
# reference-free effects for binary genotype spaces. See the example script
# RFA-binary-fast for more information.
#
# Written in R version 4.3.1.
# Packages used and their version:

library(glmnet) # 4.1-7


# Import genotype and phenotype data.
import_data <- function(dir_genotypes, dir_phenotypes) {
  
  # dir_genotypes : character
  #   Path to a csv file listing the genotypes (no header).
  # dir_phenotypes : character
  #   Path to a text file listing the phenotypes.
  #
  # The following global variables are created:
  #
  # genotypes (int matrix) : Genotype matrix 
  # y (numeric vector) : Phenotype vector
  # n (int) : Number of sites
  # q_list (int vector) : Number of states in each site
  #
  # States are renumbered as 1 and 2.
  
  genotypes <<- unname(as.matrix(read.csv(dir_genotypes, header = FALSE)))
  y <<- as.numeric(readLines(dir_phenotypes))
  
  if(!all(is.finite(genotypes)))
    stop('Genotype matrix contains missing entries.\n')
  
  if(!all(is.finite(y))) stop('Phenotype vector contains missing entries\n')
  
  if(nrow(genotypes) != length(y))
    stop('Genotype matrix and phenotype vector of different lengths.\n')
  
  n <<- ncol(genotypes) # Number of sites
  
  # Number of states in each site.
  q_list <<- apply(genotypes, 2, function(x) length(unique(x)))
  
  if(!all(q_list == 2)) stop('Genotype space is not binary!\n')
  
  # Recoding states if they are not 1 and 2.
  
  # States in each site.
  states_list <- apply(genotypes, 2, function(x) sort(unique(x)), simplify = FALSE)
  
  # Checking for unused states.
  reducible <- vapply(states_list, function(x) length(x) != max(x), logical(1))
  
  if(any(reducible)) {
    
    cat('States should be numbered 1 and 2.\n')
    cat('States renumbered for sites', which(reducible), '\n')
    genotypes <<- vapply(1:ncol(genotypes),
                         function(j) match(genotypes[, j], states_list[[j]]), genotypes[, 1])
  }
}

# Import a list of site-combinations to include in the model.
import_site_combn <- function(dir) {
  
  # dir : character
  #   Path to a text file listing the site-combinations, line by line.
  
  c(list(integer(0L)), lapply(strsplit(readLines(dir), split = ','), as.integer))
}

# Enumerate all site-combinations up to a given order.
generate_site_combn <- function(model_ord) {
  
  # model_ord : int
  #   Maximal site-combination order.
  
  do.call(c, lapply(0L:model_ord, function(k) combn(n, k, simplify = FALSE)))
}

# Construct the phenotype operator for given genotypes and reference-free
# effects. This is a matrix that, when multiplied to a vector of reference-free
# effects, computes the phenotypes.
construct_phenotype_operator <- function(genotypes, site_combn, intercept = TRUE) {
  
  # genotypes : int matrix
  #   Genotype matrix.
  # site_combn : list
  #   Site-combinations included in the model.
  # intercept : logical
  #   If FALSE, the intercept is excluded from the operator. This is needed to
  #   use the operator in glmnet.
  #
  # Unlike in reference-free analysis, the phenotype operator in simplex
  # encoding is a dense matrix where every element is either 1 or -1.
  
  if(!intercept) site_combn <- site_combn[-1]
  
  do.call(cbind, lapply(site_combn, function(sites) {
    
    if(length(sites) == 0L)
      rep(1, nrow(genotypes))
    else if(length(sites) == 1L)
      2 * (genotypes[, sites] - 1) - 1
    else
      -2 * ((rowSums(genotypes[, sites]) - length(sites)) %% 2) + 1
  }))
}

# Transform genetic score (s) by the sigmoid nonspecific epistasis.
apply_sigmoid <- function(s, L, R) L + R / (1 + exp(-s))

# Construct the generalized logit link function to use in glmnet.
glogit <- function(L, R) {
  
  # L, R : numeric
  #   Lower bound and range of the sigmoid.
  
  # In the formalism of generalized linear model, mu is the predicted expected
  # value of the response variable (phenotype), and eta is the linear predictor
  # (genetic score). They are related by
  #
  #   mu = linkinv(eta),
  #   eta = linkfun(mu).
  #
  # Note that what we call the sigmoid link function for incorporating
  # nonspecific epistasis is the inverse link function (linkinv) in the formalism
  # of generalized linear model.
  
  linkfun <- function(mu) {
    
    x <- (mu - L) / ((L + R) - mu)
    x[mu < L] <- 0
    x[mu > (L + R)] <- Inf
    log(x)
  }
  
  linkinv <- function(eta) L + R / (1 + exp(-eta))
  
  # Derivative of mu with respect to eta.
  mu.eta <- function(eta) R / ((1 + exp(-eta)) * (1 + exp(eta)))
  
  # Indicator for whether eta is in the domain of linkinv; this is always TRUE
  # for the sigmoid.
  valideta <- function(eta) TRUE
  
  # Creat an instance of class link-glm.
  structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
                 valideta = valideta, name = 'glogit(L, R)'),
            class = 'link-glm')
}

# Infer reference-free effects.
infer_model <- function(genotypes, y, site_combn, lambda = 0) {
  
  # genotypes : numeric matrix
  #   Genotype matrix.
  # y : numeric vector
  #   Phenotype vector.
  # site_combn : list
  #   Site-combinations to model.
  # lambda : numeric
  #   Lasso penalty.
  #
  # Returns a list with the following elements.
  # sigmoid : numeric vector
  #   The sigmoid parameters used.
  # e : (numeric vector) The inferred reference-free effects.
  # R2 : (numeric) Model fit.
  # residual : (numeric vector) Residual for each genotype.
  
  # Phenotype operator matrix without the intercept.
  X <- construct_phenotype_operator(genotypes, site_combn, intercept = FALSE)
  
  model <- glmnet(X, y, family = gaussian(link = glogit(L, R)), lambda = lambda)
  
  e <- unname(c(model$a0, as.vector(model$beta))) # Inferred effects
  
  s <- as.vector(X %*% model$beta) + model$a0 # Genetic score
  p <- apply_sigmoid(s, L, R) # Predicted phenotype
  
  list(sigmoid = c(L, R),
       e = e,
       R2 = 1 - sum((y - p) ^ 2) / sum((y - mean(y)) ^ 2),
       residual = y - p)
}

# Optimize the sigmoid link given a set of reference-free effects.
refine_sigmoid <- function(genotypes, y, site_combn, e, update = FALSE) {
  
  # genotypes, y, site_combn : as defined in `infer_model`
  # e : numeric vector
  #   Reference-free effects.
  # update : logical
  #   If TRUE, the global variables L and R are updated.
  
  # Phenotype operator and the genetic score of each genotype.
  X <- construct_phenotype_operator(genotypes, site_combn)
  s <- as.numeric(X %*% e)
  
  cat('\tInitial sigmoid: L = ', round(L, 3), ', R = ', round(R, 3), '\n',
      sep = '')
  
  # Initial R2.
  initR2 <- 1 - sum((y - apply_sigmoid(s, L, R)) ^ 2L) / sum((y - mean(y)) ^ 2L)
  
  # Objective function to minimize: the sum of squared errors.
  objective <- function(param) { # param = c(L, R)

    sum((y - param[1] - param[2] / (1 + exp(-s))) ^ 2)
  }
  
  # Gradient of the objective with respect to L and R.
  gradient <- function(param) {
    
    z <- (y - param[1] - param[2] / (1 + exp(-s)))
    c(sum(-2 * z), sum(-2 * z / (1 + exp(-s))))
  }
  
  res <- lbfgs(objective, gradient, c(L, R), invisible = 1)$par
  
  # Optimized R2.
  optR2 <- 1 - sum((y - apply_sigmoid(s, res[1], res[2])) ^ 2L) /
    sum((y - mean(y)) ^ 2L)
  
  cat('\tOptimal sigmoid: L = ', round(res[1], 3), ', R = ', round(res[2], 3),
      '\n', sep = '')
  cat('\tR2 improved from', round(initR2, 3), 'to', round(optR2, 3), '\n')
  
  if(update) {
    cat('\tGlobal variables L and R updated.\n')
    L <<- res[1]
    R <<- res[2]
  }
}

# Sequentially update the sigmoid and reference-free effects.
sequential_update <- function(genotypes, y, site_combn, model, epsilon) {
  
  # genotypes, y, site_combn : As defined in `infer_model`
  # model : list
  #   Initial model to begin the sequential update with.
  # epsilon : numeric
  #   Update is terminated when the model R2 does not increase by at least
  #   epsilon.
  
  iter <- 1L
  
  repeat {
    
    cat('Iteration', iter, '\n')
    
    init_model <- model
    init_L <- L
    init_R <- R
    
    refine_sigmoid(genotypes, y, site_combn, model$e, update = TRUE)
    model <- infer_model(genotypes, y, site_combn)
    
    if(model$R2 < init_model$R2) {
      
      cat('Model R2 decreased; update reverted and terminated.\n')
      L <<- init_L
      R <<- init_R
      break
      
    } else if(model$R2 - init_model$R2 < epsilon) {
      
      cat('Model R2 did not increase sufficiently; terminating.\n')
      break
      
    } else {
      iter <- iter + 1
    }
  }
  
  infer_model(genotypes, y, site_combn)
}

# Show nonspecific epistasis and the distribution of genetic score.
show_latent_space <- function(model, genotypes, y, site_combn) {
  
  # model : list
  #   Model to evaluate.
  # genotypes, y, site_combn : As defined in `infer_model`
  
  X <- construct_phenotype_operator(genotypes, site_combn) # Phenotype operator
  s <- as.vector(X %*% model$e) # Genetic score
  
  h <- hist(s, breaks = 50L, plot = FALSE)
  plot(h$mids, h$density, type = 'l', xlim = range(s), xlab = 'Genetic score',
       ylab = 'Fraction of genotypes')
  
  par(new = TRUE)
  t <- seq(min(s), max(s), length.out = 250L)
  plot(t, apply_sigmoid(t, model$sigmoid[1], model$sigmoid[2]), xlim = range(s),
       type = 'l', col = 'red', axes = FALSE, ylab = NA)
  axis(4, signif(apply_sigmoid(c(min(s), max(s)), model$sigmoid[1], model$sigmoid[2]), 2L),
       col = 'red', col.axis = 'red')
  
  legend('topleft', legend = c('Genetic score', 'Nonspecific epistasis'),
         text.col = c('black', 'red'))
}

# Export model.
export_model <- function(path, n, q_list, site_combn, model) {
  
  # path : character
  #   Path to the export file. 
  
  f <- file(path)
  open(f, 'w')
  
  writeLines('n', f)
  write(n, f)
  write('q_list', f)
  write(paste(q_list, collapse = ' '), f)
  write('site_combn', f)
  for(x in site_combn[-1L]) write(paste(x, collapse = ' '), f)
  write('nonspec_type', f)
  write(2, f)
  write('nonspec_param', f)
  write(paste(signif(model$sigmoid, 5L), collapse = ' '), f)
  write('e', f)
  write(as.character(signif(model$e, 5L)), f)
  
  close(f)
}

# Import model.
import_model <- function(path) {
  
  # path : character
  #   Directory to the model file.
  #
  # The following global variables are overwritten and therefore must be
  # compatible with the genotype and phenotype data:
  # n, q_list, site_combn
  
  file <- readLines(path)
  
  n <<- as.integer(file[which(file == 'n') + 1L]) # Number of sites
  
  # Number of states per site.
  q_list <<- as.integer(strsplit(file[which(file == 'q_list') + 1L], split = ' ')[[1L]])
  
  # Site-combinations.
  site_combn <- strsplit(file[(which(file == 'site_combn') + 1L):(which(file == 'nonspec_type') - 1L)], split = ' ')
  site_combn <<- c(list(integer(0L)), lapply(site_combn, function(x) as.integer(x)))
  
  # Type of nonspecific epistasis.
  nonspec_type <- as.integer(file[which(file == 'nonspec_type') + 1L])
  
  # Sigmoid lower bound (L) and range (R).
  nonspec_param <- as.numeric(strsplit(file[which(file == 'nonspec_param') + 1L], split = ' ')[[1L]])
  L <<- nonspec_param[1]
  R <<- nonspec_param[2]
  
  # Reference-free effects.
  e <- as.numeric(file[(which(file == 'e') + 1L):length(file)])
  
  # Returning model.
  list(sigmoid = nonspec_param, e = e)
}

# Perform k-fold cross-validation.
cross_validation <- function(genotypes, y, site_combn, k, niter) {
  
  # genotypes, y, site_combn : As defined in `infer_model`
  # k : int
  #   Number of genotype partitions for cross-validation.
  # niter : int
  #   Number of times the k-fold cross-validation should be performed.
  #
  # Returns the out-of-sample R2 values for each lambda and iteration.
  
  # Phenotype operator matrix without the intercept.
  X <- construct_phenotype_operator(genotypes, site_combn, intercept = FALSE)
  
  cvR2 <- list()
  
  for(iter in 1:niter) {
    
    cat('Iteration', iter, '\n')
    
    cv <- cv.glmnet(X, y, family = gaussian(link = glogit(L, R)), trace.it = 1,
                    type.measure = 'mse', nfolds = k)
    
    cvR2[[iter]] <- 1 - cv$cvm / mean((y - mean(y)) ^ 2) # Out-of-sample R2
  }

  cvR2 <- cbind(cv$lambda, do.call(cbind, cvR2))
  colnames(cvR2) <- c('lambda', 1:niter)
  
  cvR2
}

# Perform cross-validation to evaluate a series of nested models constructed by
# including reference-free effects in the order of their phenotypic contribution.
minimal_models <- function(genotypes, y, site_combn, model, max_size, k) {
  
  # genotypes, y, site_combn : As defined in `infer_model`
  # model : list
  #   Model inferred from genotypes, y, and site_combn; used to rank the effects
  #   by their phenotypic contribution.
  # max_size : int
  #   Maximum number of effects to evaluate.
  # k : int
  #   Number of genotype partitions to use for cross-validation.
  #
  # Returns the out-of-sample R2 of the minimal models.
  
  # Ordering reference-free effects by phenotypic contribution.
  
  # Fraction of phenotypic variance due to each effect.
  v <- model$e[-1] ^ 2 / sum(model$e[-1] ^ 2)
  
  # Order of effects to include.
  ord <- order(v, decreasing = TRUE)
  
  # Full phenotype operator.
  X <- construct_phenotype_operator(genotypes, site_combn, intercept = FALSE)
  
  # Lasso is activated only when the number of parameters is sufficiently large.
  perform_lasso <- FALSE
  
  cat('Intercept: R2 = 0\n') # Intercept-only model
  
  # glmnet does not support a single-parameter model.
  cat('Effect 1: R2 = NA (glmnet cannot fit just two parameters)\n')
  minR2 <- c(0, NA_real_)
  
  for(i in 2:min(length(v), max_size)) {
    
    cat('Effect number ', i, ':', sep = '')
    
    if(perform_lasso) {
      
      cv <- cv.glmnet(X[, ord[1:i]], y, family = gaussian(link = glogit(L, R)),
                      type.measure = 'mse', nfolds = k)
      
      minR2 <- c(minR2, max(1 - cv$cvm / mean((y - mean(y)) ^ 2)))
      
    } else {
      
      if(i %% 20L == 0L) { # Perform lasso with each 20 additional parameters
        
        cv <- cv.glmnet(X[, ord[1:i]], y, family = gaussian(link = glogit(L, R)),
                        type.measure = 'mse', nfolds = k)
        
        cvR2 <- 1 - cv$cvm / mean((y - mean(y)) ^ 2)
        minR2 <- c(minR2, max(cvR2))
        
        # If better model fit is found for a non-minimal lambda, begin lasso.
        if(max(cvR2) - cvR2[length(cvR2)] > 0.002) {
          
          cat('  Performing lasso regression from now on\n')
          perform_lasso <- TRUE
        }
        
      } else {
        
        cv <- cv.glmnet(X[, ord[1:i]], y, family = gaussian(link = glogit(L, R)),
                        lambda = c(1e-5, 0), type.measure = 'mse', nfolds = k)
        
        minR2 <- c(minR2, 1 - cv$cvm[2] / mean((y - mean(y)) ^ 2))
      }
    }
    
    cat(' R2 = ', round(minR2[length(minR2)], 3), '\n')
  }
  
  unname(minR2)
}

