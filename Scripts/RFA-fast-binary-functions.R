
# This script implements the functions for fast, approximate inference of
# reference-free effects for binary genotype spaces. See the example script
# RFA-binary-fast for more information.
#
# Written in R version 4.3.1.
# Packages used and their version:

library(glmnet) # 4.1-7

# To add more link functions for modeling nonspecific epistasis, edit the
# following functions as instructed in each function implementation:
#
# custom_link, apply_nonspecific_epistasis, refine_nonspecific_epistasis.

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
  reducible <- vapply(states_list, function(x) length(x) != max(x), logical(1L))
  
  if(any(reducible)) {
    
    cat('States should be numbered 1 and 2.\n')
    cat('States renumbered for sites', which(reducible), '\n')
    genotypes <<- vapply(1:ncol(genotypes),
                         function(j) match(genotypes[, j], states_list[[j]]), genotypes[, 1L])
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

# Custom link functions for use in glmnet.
custom_link <- function() {
  
  # To add more link functions, create a function analogous to 'glogit'
  # following the instructions therein. Then return a function call using the
  # if-else loop at the end.
  
  # The sigmoid link function.
  glogit <- function(L, R) {
    
    # L, R : numeric
    #   Lower bound and range of the sigmoid.
    #
    # In the formalism of generalized linear model, mu is the predicted expected
    # value of the response variable (phenotype), and eta is the
    # linear predictor (genetic score). They are related by
    #
    #   mu = linkinv(eta),
    #   eta = linkfun(mu).
    #
    # Note that what we call the sigmoid link function for incorporating
    # nonspecific epistasis is the inverse link function (linkinv) in the
    # formalism of generalized linear model.
    
    # Link function (inverse of the sigmoid).
    linkfun <- function(mu) {
      
      x <- (mu - L) / ((L + R) - mu)
      x[mu < L] <- 0
      x[mu > (L + R)] <- Inf
      log(x)
    }
    
    # Inverse link function (sigmoid).
    linkinv <- function(eta) L + R / (1 + exp(-eta))
    
    # Derivative of mu with respect to eta.
    mu.eta <- function(eta) R / ((1 + exp(-eta)) * (1 + exp(eta)))
    
    # Indicator for whether eta is in the domain of linkinv; this is always TRUE
    # for the sigmoid.
    valideta <- function(eta) TRUE
    
    # Creat an instance of class link-glm.
    structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
                   valideta = valideta, name = 'glogit'),
              class = 'link-glm')
  }
  
  if(nonspec_type == 'sigmoid')
    glogit(nonspec_param[1L], nonspec_param[2L])
  else
    stop(paste("Link function '", nonspec_type, "' unavailable.", sep = ''))
}

# Transform genetic score by nonspecific epistasis.
apply_nonspecific_epistasis <- function(s, nonspec_param) {
  
  # s : numeric vector
  #   Genetic scores.
  # nonspec_param : numeric vector
  #   Nonspecific epistasis parameters to use.
  
  # Add more link functions in this if-else loop.
  if(nonspec_type == 'none')
    s
  else if(nonspec_type == 'sigmoid')
    nonspec_param[1L] + nonspec_param[2L] / (1 + exp(-s))
  else
    stop(paste("Link function '", nonspec_type, "' unavailable.", sep = '')) 
}

# Infer reference-free effects.
infer_model <- function(genotypes, y, site_combn, lambda = 0, norm_e = TRUE) {
  
  # genotypes : numeric matrix
  #   Genotype matrix.
  # y : numeric vector
  #   Phenotype vector.
  # site_combn : list
  #   Site-combinations to model.
  # lambda : numeric
  #   Lasso penalty.
  # norm_e : logical
  #   If TRUE, the zero-mean constraint is enforced on the returned model.
  #
  # Returns a list with the following elements.
  # nonspec_param : (numeric vector) The nonspecific epistasis parameters used
  # e : (numeric vector) The inferred reference-free effects
  # R2 : (numeric) Model fit
  # residual : (numeric vector) Residual for each genotype
  
  # Defining the glmnet family.
  if(nonspec_type == 'none')
    family <- 'gaussian'
  else
    family <- gaussian(link = custom_link())
  
  # Phenotype operator matrix without the intercept.
  X <- construct_phenotype_operator(genotypes, site_combn, intercept = FALSE)
  
  model <- glmnet(X, y, family = family, lambda = lambda)
  
  e <- unname(c(model$a0, as.vector(model$beta))) # Inferred effects
  
  if(nonspec_type == 'none') nonspec_param <<- NULL
  
  s <- as.vector(X %*% model$beta) + model$a0 # Genetic score
  p <- apply_nonspecific_epistasis(s, nonspec_param) # Predicted phenotype
  
  list(nonspec_param = nonspec_param,
       e = e,
       R2 = 1 - sum((y - p) ^ 2) / sum((y - mean(y)) ^ 2),
       residual = y - p)
}

# Optimize the nonspecific epistasis parameters while fixing the reference-free
# effects.
refine_nonspecific_epistasis <- function(genotypes, y, site_combn, e, update) {
  
  # genotypes, y, site_combn : as defined in `infer_model`
  # e : numeric vector
  #   Reference-free effects.
  # update : logical
  #   If TRUE, the global variable 'nonspec_param' is updated.
  
  if(nonspec_type == 'none') {
    
    cat('No nonspecific epistasis to optimize.\n')
    return(NULL)
  }
  
  # Phenotype operator and the genetic score of each genotype.
  X <- construct_phenotype_operator(genotypes, site_combn)
  s <- as.numeric(X %*% e)
  
  cat('\tInitial parameters: '); cat(round(nonspec_param, 3), sep = ', ') 
  
  # Initial R2.
  initR2 <- 1 - sum((y - apply_nonspecific_epistasis(s, nonspec_param)) ^ 2) / 
    sum((y - mean(y)) ^ 2)
  
  # Objective function to minimize: the sum of squared errors.
  objective <- function(param) {
    
    if(nonspec_type == 'sigmoid') {
      
      sum((y - param[1L] - param[2L] / (1 + exp(-s))) ^ 2)
      
    } else { # Add more link functions in this if-else loop.
      
      stop(paste("Link function '", nonspec_type, "' unavailable.", sep = '')) 
    }
  }
  
  # Gradient of the objective with respect to the parameters.
  gradient <- function(param) {
    
    if(nonspec_type == 'sigmoid') {
      
      z <- (y - param[1L] - param[2L] / (1 + exp(-s)))
      c(sum(-2 * z), sum(-2 * z / (1 + exp(-s))))
      
    } else { # Add more link functions in this if-else loop.
      
      stop(paste("Link function '", nonspec_type, "' unavailable.", sep = ''))
    }
  }
  
  res <- lbfgs(objective, gradient, nonspec_param, invisible = 1)$par
  
  # Optimized R2.
  optR2 <- 1 - sum((y - apply_nonspecific_epistasis(s, res)) ^ 2) /
    sum((y - mean(y)) ^ 2)
  
  cat('\n\tOptimized parameters: '); cat(round(res, 3), sep = ', ')
  cat('\n\tR2 improved from', round(initR2, 3), 'to', round(optR2, 3), '\n')
  
  if(update) {
    cat('\tGlobal variable "nonspec_param" updated.\n')
    nonspec_param <<- res
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
    init_nonspec_param <- nonspec_param
    
    refine_nonspecific_epistasis(genotypes, y, site_combn, model$e, TRUE)
    model <- infer_model(genotypes, y, site_combn)
    
    if(model$R2 < init_model$R2) {
      
      cat('Model R2 decreased; update reverted and terminated.\n')
      nonspec_param <<- init_nonspec_param
      return(init_model)
      
    } else if(model$R2 - init_model$R2 < epsilon) {
      
      cat('Model R2 did not increase sufficiently; terminating.\n')
      return(model)
      
    } else {
      iter <- iter + 1
    }
  }
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
  plot(t, apply_nonspecific_epistasis(t, model$nonspec_param), xlim = range(s),
       type = 'l', col = 'red', axes = FALSE, ylab = NA)
  axis(4, signif(apply_nonspecific_epistasis(c(min(s), max(s)), model$nonspec_param), 2L),
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
  write(nonspec_type, f)
  
  write('nonspec_param', f)
  if(nonspec_type == 'none')
    write('NA', f)
  else
    write(paste(signif(model$nonspec_param, 5L), collapse = ' '), f)
  
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
  # n, q_list, site_combn, nonspec_type.
  
  file <- readLines(path)
  
  n <<- as.integer(file[which(file == 'n') + 1L]) # Number of sites
  
  # Number of states per site.
  q_list <<- as.integer(strsplit(file[which(file == 'q_list') + 1L], split = ' ')[[1L]])
  
  # Site-combinations.
  site_combn <- strsplit(file[(which(file == 'site_combn') + 1L):(which(file == 'nonspec_type') - 1L)], split = ' ')
  site_combn <<- c(list(integer(0L)), lapply(site_combn, function(x) as.integer(x)))
  
  # Type of nonspecific epistasis.
  nonspec_type <<- file[which(file == 'nonspec_type') + 1L]
  
  # Nonspecific epistasis parameters.
  if(nonspec_type == 'none')
    nonspec_param <<- NULL
  else
    nonspec_param <<- as.numeric(strsplit(file[which(file == 'nonspec_param') + 1L], split = ' ')[[1L]])
  
  # Reference-free effects.
  e <- as.numeric(file[(which(file == 'e') + 1L):length(file)])
  
  # Returning model.
  list(nonspec_param = nonspec_param, e = e)
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
  
  # Defining the glmnet family.
  if(nonspec_type == 'none')
    family <- 'gaussian'
  else
    family <- gaussian(link = custom_link())
  
  # Phenotype operator matrix without the intercept.
  X <- construct_phenotype_operator(genotypes, site_combn, intercept = FALSE)
  
  cvR2 <- list()
  
  for(iter in 1:niter) {
    
    cat('Iteration', iter, '\n')
    
    
    cv <- cv.glmnet(X, y, family = family, trace.it = 1, type.measure = 'mse',
                    nfolds = k)
    
    cvR2[[iter]] <- 1 - cv$cvm / mean((y - mean(y)) ^ 2) # Out-of-sample R2
  }
  
  cvR2 <- cbind(cv$lambda, do.call(cbind, cvR2))
  colnames(cvR2) <- c('lambda', 1L:niter)
  
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
  
  # Defining the glmnet family.
  if(nonspec_type == 'none')
    family <- 'gaussian'
  else
    family <- gaussian(link = custom_link())
  
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
      
      cv <- cv.glmnet(X[, ord[1:i]], y, family = family, type.measure = 'mse',
                      nfolds = k)
      
      minR2 <- c(minR2, max(1 - cv$cvm / mean((y - mean(y)) ^ 2)))
      
    } else {
      
      if(i %% 20L == 0L) { # Perform lasso with each 20 additional parameters
        
        cv <- cv.glmnet(X[, ord[1:i]], y, family = family, type.measure = 'mse',
                        nfolds = k)
        
        cvR2 <- 1 - cv$cvm / mean((y - mean(y)) ^ 2)
        minR2 <- c(minR2, max(cvR2))
        
        # If better model fit is found for a non-minimal lambda, begin lasso.
        if(max(cvR2) - cvR2[length(cvR2)] > 0.002) {
          
          cat('  Performing lasso regression from now on\n')
          perform_lasso <- TRUE
        }
        
      } else {
        
        cv <- cv.glmnet(X[, ord[1:i]], y, family = family, lambda = c(1e-5, 0),
                        type.measure = 'mse', nfolds = k)
        
        minR2 <- c(minR2, 1 - cv$cvm[2] / mean((y - mean(y)) ^ 2))
      }
    }
    
    cat(' R2 = ', round(minR2[length(minR2)], 3), '\n')
  }
  
  unname(minR2)
}

