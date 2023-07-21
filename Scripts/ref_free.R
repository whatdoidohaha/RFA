
# This script implements functions for reference-free analysis of any form of sequence space.

# Libraries required.

library(lbfgs)
library(Matrix)
library(abind)
library(gplots)


# Import model settings from a control file.
import_control_file <- function(dir) {
  
  # Control file must contain four lines:
  #   1. Path to a csv file listing the genotypes.
  #   2. Path to a text file listing the phenotypes (line by line for each genotype).
  #   3. Integer specifycing the model order or path to a text file listing the site-combinations to model.
  #   4. Integer specifying the type of nonspecific epistasis (1 = None; 2 = Logistic).
  
  input <- readLines(dir)
  control <- lapply(strsplit(input, split = '#'), function(x) trimws(x[1L], 'right'))
  control[-c(1L:3L)] <- lapply(control[-c(1L:3L)], as.integer)
  control
}

# Set model as specified in a control file.
set_model <- function(control) {
  
  # The following global variables are created:
  #
  # genotypes : Genotype matrix.
  # y : Phenotype vector.
  # n : Number of sites.
  # q_list : Number of states in each site.
  # site_combn : List of site-combinations for which effects are inferred.
  # nonspec_type : Type of nonspecific epistasis (1 = none, 2 = logistic).
  #
  # States are renumbered to be 1, 2, 3, ...
  
  genotypes <<- unname(as.matrix(read.csv(control[[1L]], header = FALSE)))
  if(!all(is.finite(genotypes))) stop('Genotype matrix contains missing entries.\n')
  
  y <<- as.numeric(readLines(control[[2L]]))
  if(!all(is.finite(y))) stop('Phenotype vector contains missing entries\n')
  if(nrow(genotypes) != length(y)) stop('Genotype matrix and phenotype vector of different lengths.\n')
  
  n <<- ncol(genotypes) # Number of sites
  q_list <<- apply(genotypes, 2L, function(x) length(unique(x))) # Number of states in each site
  
  # Site-combinations to model.
  
  # All site-combinations up to the specified model order:
  if(all(strsplit(control[[3L]], '')[[1L]] %in% as.character(0L:9L))) {
    
    site_combn <<- do.call(c, lapply(0L:as.integer(control[[3L]]), function(k) combn(n, k, simplify = FALSE)))
    
  } else { # Specific site-combinations.
    
    site_combn <<- c(list(integer(0L)), lapply(strsplit(readLines(control[[3L]]), split = ','), as.integer))
  }
  
  # Renumbering states if necessary.
  
  ori_q_list <- lapply(1L:ncol(genotypes), function(j) sort(unique(genotypes[, j]))) # States in each site
  
  reducible <- vapply(ori_q_list, function(x) length(x) != max(x), TRUE) # Sites to renumber states.
  
  if(any(reducible)) {
    
    cat('States should be numbered 1, 2, 3, ...\n')
    cat('States renumbered for sites', which(reducible), '\n')
    genotypes <<- vapply(1L:ncol(genotypes), function(j) match(genotypes[, j], ori_q_list[[j]]), genotypes[, 1L])
  }
  
  # Type of global nonlinearity.
  nonspec_type <<- control[[4L]]
}

# Enumerate all site-combinations up to a given order, including the intercept.
generate_site_combn <- function(model_ord) {
  
  do.call(c, lapply(0L:model_ord, function(k) combn(n, k, simplify = FALSE)))
}

# Construct the phenotype operator matrix for the given genotypes and reference-free effects.
# This is the matrix that, when multiplied to the vector of reference-free effects, yields the phenotypes.
construct_phenotype_operator <- function(genotypes, site_combn) {
  
  # genotypes : numeric matrix
  #   Genotype matrix.
  # site_combn : list
  #   Site-combinations to model.
  #
  # The phenotype operator matrix is a sparse matrix in which elements are either 0 or 1.
  # The sparse matrix implemention of the library Matrix is used.
  
  ord <- vapply(site_combn, length, 0L) # Order of each site-combination
  n_terms <- as.integer(vapply(site_combn, function(x) prod(q_list[x]), 0)) # Number of terms in each site-combination
  cum_terms <- cumsum(c(0L, n_terms[-length(n_terms)])) # Cumulative number of terms
  
  # Row indices of the elements to be set as 1.
  a <- rep(1L:nrow(genotypes), length(site_combn))
  
  # Column indices of the elements to be set as 1.
  b <- vapply(1L:length(site_combn), function(k) {
    
    if(ord[k] == 0L) {
      
      rep(1L, nrow(genotypes))
      
    } else if(ord[k] == 1L) {
      
      cum_terms[k] + genotypes[, site_combn[[k]]]
      
    } else {
      
      # Within a site-combination, state-combinations are ordered by last site moving fastest.
      q_k <- q_list[site_combn[[k]]]
      cum_q <- as.integer(c(rev(cumprod(rev(q_k[-1L]))), 1L))
      cum_terms[k] + as.integer(genotypes[, site_combn[[k]]] %*% cum_q) + (1L - sum(cum_q))
    }
    
  }, integer(nrow(genotypes)))
  
  sparseMatrix(a, b, dims = c(nrow(genotypes), sum(n_terms)))
}

# Transform genetic score by nonspecific epistasis.
apply_nonspecific_epistasis <- function(s, nonspec_param) {
  
  # s : numeric vector
  #   Genetic score.
  # nonspec_param : vector
  #   Parameters of the nonspecific epistasis function.
  
  if(nonspec_type == 1L)
    s
  else if(nonspec_type == 2L)
    nonspec_param[1L] + nonspec_param[2L] / (1 + exp(-s))
}

# Sort reference-free effects by their site-combination.
parse_effects_by_site_combn <- function(e, site_combn, as_array = FALSE) {
  
  # e : numeric vector
  #   Vector of reference-free effects to be sorted.
  # site_combn : list
  #   List of site-combinations corresponding to e.
  # as_array : logical
  #   If TRUE, state-combinations within each site-combination are arranged in an array.
  #   Otherwise, they are arranged in a vector.
  
  # Identifying reference-free effects for each site-combination.
  
  ord <- vapply(site_combn, length, 0L) # Order of each site-combination
  n_terms <- vapply(site_combn, function(x) prod(q_list[x]), 0) # Number of terms in each site-combination
  pos <- cumsum(c(0L, n_terms)) # Index of the last term of each site-combination
  
  e_parsed <- lapply(1L:length(site_combn), function(i) e[(pos[i] + 1L):pos[i + 1L]])
  
  if(as_array) { # Arranging state-combinations in an array, indexed in the order of sites.
    
    for(i in which(ord > 1L)) {
      
      sites <- site_combn[[i]]
      qs <- q_list[sites]  
      E <- array(e_parsed[[i]], rev(qs))
      e_parsed[[i]] <- aperm(E, length(sites):1L)
    }
  }
  
  e_parsed
}

# Enforce the zero-mean constraint.
normalize_e <- function(e, site_combn, genotypes) {
  
  # e : numeric vector
  #   Reference-free effects to be normalized.
  # site_combn : list
  #   List of site-combinations corresponding to e.
  # genotypes : integer matrix
  #   Genotypes from which e was inferred.
  #   This information is needed for the numerical algorithm implemented here.
  #
  # The algorithm is explained in the supplementary text appendix 5.
  
  # Enforce zero-mean constraint on all k-th order terms, leaving all others constant.
  normalize_e_k <- function(e, site_combn, k) {
    
    ord <- vapply(site_combn, length, 0L) # Order of each site-combination
    if(sum(ord == k) == 0L) return (e)
    
    e_parsed <- parse_effects_by_site_combn(e, site_combn)
    
    if(k == 1L) {
      
      e_parsed[[1L]] <- e_parsed[[1L]] + sum(vapply(which(ord == 1L), function(i) mean(e_parsed[[i]]), 0))
      
      for(i in which(ord == 1L))
        e_parsed[[i]] <- e_parsed[[i]] - mean(e_parsed[[i]])
      
    } else {
      
      for(i in which(ord == k)) {
        
        sites <- site_combn[[i]]
        qi <- q_list[sites]
        rev_qi <- rev(qi)
        
        E <- array(e_parsed[[i]], rev_qi)
        D <- E
        sgn <- 1
        
        for(low_ord in (k - 1L):1L) {
          
          correction_matrices <- lapply(combn(k, low_ord, simplify = FALSE), function(margin) {
            
            reduced <- apply(E, margin, mean)
            expanded <- array(reduced, c(rev_qi[margin], rev_qi[-margin]))
            aperm(expanded, order(c(margin, setdiff(1L:k, margin))))
          })
          
          D <- D - sgn * Reduce(`+`, correction_matrices)
          sgn <- sgn * -1
        }
        
        e_parsed[[i]] <- as.vector(D) - sgn * mean(E)
      }
    }
    
    unlist(e_parsed)
  }
  
  ord <- vapply(site_combn, length, 0L) # Order of each site combination
  
  if(max(ord) == 1L) return(normalize_e_k(e, site_combn, 1L))
  
  G <- construct_phenotype_operator(genotypes, site_combn)
  s <- as.vector(G %*% e) # Genetic score
  n_terms <- vapply(site_combn, function(x) prod(q_list[x]), 0) # Number of terms in each site-combination.
  e_ord <- unlist(lapply(1L:length(site_combn), function(i) rep(ord[i], n_terms[i]))) # Order of each effect
  
  for(k in max(ord):2L) {
    
    e <- normalize_e_k(e, site_combn, k)

    G_l_k <- G[, e_ord < k]
    G_ge_k <- G[, e_ord >= k]
    e_l_k <- e[e_ord < k]
    e_ge_k <- e[e_ord >= k]
    
    # Objective function for correcting the lower-order effects.
    objective <- function(e_l_k) {
      
      z0 <<- s - as.vector(G_l_k %*% e_l_k + G_ge_k %*% e_ge_k)
      sum(z0 ^ 2L)
    }
    
    # Gradient function for correcting the lower-order effects.
    gradient <- function(e_l_k) as.vector(-2 * z0 %*% G_l_k)
    
    e[e_ord < k] <- lbfgs(objective, gradient, e_l_k, invisible = 1)$par
  }

  normalize_e_k(e, site_combn, 1L)
}

# Objective function to minimize for inferring reference-free effects and nonspecific epistasis parameters.
objective <- function(param) {
  
  # param : numeric vector
  #   Parameters to evaluate; reference-free effects followed by the nonspecific epistasis parameters.
  #
  # Computes the sum of squared residuals and the lasso penalty.
  # The phenotype vector and phenotype operator matrix are globally stored (`y_active` and `G_active`).
  
  if(nonspec_type == 1L) {
    
    e <- param
    
    # Residuals; stored as a global variable to be used in the gradient calculation by the function `gradient`.
    z0 <<- y_active - as.vector(G_active %*% e)
    
    sum(z0 ^ 2L)
    
  } else if(nonspec_type == 2L) {
    
    nonspec_param <- param[(length(param) - 1L):length(param)]
    e <- param[1L:(length(param) - 2L)]
    
    # Partial residuals; stored as a global variable to be used in the gradient calculation by the function `gradient`.
    z0 <<- nonspec_param[2L] / (1 + exp(-as.vector(G_active %*% e)))
    
    sum((y_active - nonspec_param[1L] - z0) ^ 2L)
  }
}

# Gradient function to use for OWL-QN (numerical optimizer for lasso regression).
gradient <- function(param) {
  
  # The phenotype vector, phenotype operator matrix, and residuals are globally stored (`y_active`, `G_active`, and `z0`).
  
  if(nonspec_type == 1L) {
    
    -2 * as.vector(z0 %*% G_active)
    
  } else if(nonspec_type == 2) {
    
    nonspec_param <- param[(length(param) - 1L):length(param)]
    e <- param[1L:(length(param) - 2L)]
    
    z1 <- -2 * (y_active - nonspec_param[1L] - z0)
    z2 <- z1 * z0 / nonspec_param[2L]
    z3 <- z2 * (nonspec_param[2L] - z0)
    c(as.vector(z3 %*% G_active), sum(z1), sum(z2))
  }
}

# Jointly infer nonspecific epistasis and reference-free effects.
infer_model <- function(genotypes, y, site_combn, lambda = 0, init_model = NULL,
                        norm_e = TRUE, silent = 1L, epsilon = 1e-5) {
  
  # genotypes : numeric matrix
  #   Genotype matrix.
  # y : numeric vector
  #   Phenotype vector.
  # site_combn : list
  #   Site-combinations to model.
  # lambda : numeric
  #   Lasso penalty.
  # init_model : list (optional)
  #   Model to use as initial value for L-BFGS.
  # norm_e : logical
  #   If TRUE, zero-mean constraint is enforced on the returned reference-free effects.
  # silent : 0 or 1 (default)
  #   If 1, the console output of lbfgs::lbfgs is suppressed.
  # epsilon : numeric
  #   Optimization precision; see the documentation for lbfgs::lbfgs for definition.
  
  # Phenotye vector and phenotype operator matrix to use for calculating the objective and gradient.
  G_active <<- construct_phenotype_operator(genotypes, site_combn)
  y_active <<- y
  
  if(var(y) < 1e-10) { # No phenotypic variation to model
    
    if(nonspec_type == 1L) {
      
      return(list(nonspec_param = NULL, e = c(mean(y), rep(0, ncol(G_active) - 1L)),
                  R2 = NA_real_, residual = rep(0, length(y))))
      
    } else if(nonspec_type == 2L) {
      
      return(list(nonspec_param = c(mean(y), 0), e = rep(0, ncol(G_active)),
                  R2 = NA_real_, residual = rep(0, length(y))))
    }
  }
  
  if(!is.null(init_model)) {
    
    init_param <- c(init_model$e, init_model$nonspec_param)
    
  } else {
    
    init_param <- rep(0, ncol(G_active)) # Initial value for every reference-free effect is zero.
    
    if(nonspec_type == 2L) # Initial nonspecific epistasis minimally covers the data range.
      init_param <- c(init_param, min(y), max(y) - min(y))
  }
  
  res <- lbfgs(objective, gradient, init_param, orthantwise_c = lambda, orthantwise_start = 1L,
               invisible = silent, epsilon = epsilon)$par
  # Nonspecific epistasis parameters are subject to regularization to avoid identifiability problem in some cases.
  
  if(norm_e)
    e <- normalize_e(res[1L:ncol(G_active)], site_combn, genotypes)
  else
    e <- res[1L:ncol(G_active)]
  
  if(nonspec_type == 1L)
    nonspec_param <- NULL
  else
    nonspec_param <- res[(ncol(G_active) + 1L):length(res)]
  
  p <- apply_nonspecific_epistasis(as.vector(G_active %*% e), nonspec_param)
  
  list(nonspec_param = nonspec_param, e = e, R2 = 1 - sum((y - p) ^ 2L) / sum((y - mean(y)) ^ 2L), residual = y - p)
}

# Infer reference-free effects under a fixed model of nonspecific epistasis.
infer_specific_effects <- function(genotypes, y, site_combn, nonspec_param, lambda = 0, init_e = NULL,
                                   norm_e = TRUE, silent = 1L, epsilon = 1e-5) {
  
  # genotypes : numeric matrix
  #   Genotype matrix.
  # y : numeric vector
  #   Phenotype vector.
  # site_combn : list
  #   Site-combinations to model.
  # nonspec_param : vector
  #   Parameters of the nonspecific epistasis function.
  # lambda : numeric
  #   Lasso penalty.
  # init_e : numeric vector
  #   Optional initial values to use in L-BFGS inference.
  # norm_e : logical
  #   If TRUE, zero-mean constraint is enforced on the returned reference-free effects.
  # silent : 0 or 1
  #   If 1, the console output of lbfgs is suppressed.
  # epsilon : numeric
  #   Optimization precision; see the documentation for lbfgs::lbfgs for definition.
  
  G_active <<- construct_phenotype_operator(genotypes, site_combn)
  y_active <<- y
  
  # Objective and gradient functions to use for this purpose.
  
  objective_specific <- function(e) objective(c(e, nonspec_param))
  
  gradient_specific <- function(e) {
    
    grad <- gradient(c(e, nonspec_param))
    
    if(nonspec_type == 1L)
      grad
    else
      grad[1L:ncol(G_active)]
  }
  
  if(is.null(init_e)) init_e <- rep(0, ncol(G_active))
  
  res <- lbfgs(objective_specific, gradient_specific, init_e,
               orthantwise_c = lambda, orthantwise_start = 1L, invisible = silent, epsilon = epsilon)$par
  
  if(norm_e) e <- normalize_e(res, site_combn, genotypes)
  
  p <- apply_nonspecific_epistasis(as.vector(G_active %*% e), nonspec_param)
  
  list(nonspec_param = nonspec_param, e = e, R2 = 1 - sum((y - p) ^ 2L) / sum((y - mean(y)) ^ 2L), residual = y - p)
}

# Evaluate a model by plotting predicted versus observed phenotype and calculating the model fit.
evaluate_model <- function(genotypes, y, site_combn, model, show_plot = TRUE, ...) {
  
  # genotypes : numeric matrix
  #   Genotype matrix.
  # y : numeric vector
  #   Phenotype vector.
  # site_combn : list
  #   Site-combinations to model.
  # model : list
  #   Model to evaluate.
  # show_plot : logical
  #   Whether a plot of predicted versus observed phenotype should be generated.
  # ... : optional parameters for `plot`.
  
  G <- construct_phenotype_operator(genotypes, site_combn)
  p <- apply_nonspecific_epistasis(as.vector(G %*% model$e), model$nonspec_param)
  
  if(show_plot) plot(p, y, cex = 0.2, col = rgb(0, 0, 0, 0.2), xlab = 'Predicted', ylab = 'Observed', ...)
  
  1 - sum((y - p) ^ 2L) / sum((y - mean(y)) ^ 2L)
}

# Show the distribution of genetic score and inferred nonspecific epistasis.
show_latent_space <- function(genotypes, y, site_combn, model, ...) {
  
  # genotypes : numeric matrix
  #   Genotype matrix.
  # y : numeric vector
  #   Phenotype vector.
  # site_combn : list
  #   Site-combinations to model.
  # model : list
  #   Model to evaluate.
  # ... : optional parameters to `plot`.
  
  G <- construct_phenotype_operator(genotypes, site_combn)
  s <- as.vector(G %*% model$e) # Genetic score
  
  t <- seq(min(s), max(s), length.out = 250L)
  plot(t, apply_nonspecific_epistasis(t, model$nonspec_param), type = 'l', xlab = 'Latent', xlim = range(s))
  
  h <- hist(s, breaks = 50L, plot = FALSE)
  par(new = TRUE)
  plot(h$mids, h$density, type = 'l', col = 'red', axes = FALSE, xlim = range(s))
  axis(4L, c(0, signif(max(h$density), 2L)), col = 'red', col.axis = 'red')
  abline(v = quantile(s, c(0.25, 0.5, 0.75)), lty = 2L)
  
  legend('topleft', legend = c('Nonspecific epistasis', 'Distr. of latent phenotypes'), text.col = c('black', 'red'))
}

# Calculate the fraction of variance due to each reference-free effect.
partition_variance <- function(model, site_combn) {
  
  # model : list
  #   Model subject to variance partition.
  # site_combn : list
  #   List of site-combinations modeled.
  
  e_parsed <- parse_effects_by_site_combn(model$e, site_combn, as_array = TRUE)[-1L]
  v <- lapply(e_parsed, function(x) x ^ 2L / length(x))
  lapply(v, function(x) x / sum(unlist(v)))
}

# Perform k-fold cross-validation by partitioning genotypes into train and test set.
cross_validation <- function(genotypes, y, site_combn, k, lambda_list,
                             niter = 1L, init_model = NULL, silent = 0L, epsilon = 1e-5) {
  
  # genotypes : numeric matrix
  #   Genotype matrix.
  # y : numeric vector
  #   Phenotype vector.
  # site_combn : list
  #   Site-combinations to model.
  # k : integer
  #   Number of partitions to use for cross-validation.
  # lambda_list : numeric vector
  #   Lasso penalty values to evaluate.
  # niter : integer
  #   Number of times cross-validation should be performed for each lasso penalty value.
  # init_model : list (optional)
  #   Model to use as initial value for L-BFGS.
  # silent : 0 (default) or 1
  #   If 1, cross-validation progression is not reported.
  # epsilon : numeric
  #   Optimization precision; see the documentation for lbfgs::lbfgs for definition.
  
  G <- construct_phenotype_operator(genotypes, site_combn)
  
  cv_cor <- matrix(NA_real_, k * niter, length(lambda_list))
  
  for(iter in 1L:niter) {
    
    membership <- sample(1L:nrow(G) %% k + 1L) # Random partition of genotypes into k categories
    
    for(j in seq_along(lambda_list)) { # j indexes lasso penalty values
      
      for(i in 1L:k) { # i indexes the k partitions
        
        if(var(y[membership == i]) < 1e-10) next # If the test set has no phenotypic variation to explain
        
        model <- infer_model(genotypes[membership != i, ], y[membership != i], site_combn, lambda_list[j],
                             init_model = init_model, norm_e = FALSE, epsilon = epsilon)
        
        # Predicted phenotype for the test set.
        p <- apply_nonspecific_epistasis(as.vector(G[membership == i, ] %*% model$e), model$nonspec_param)
        
        # Out-of-sample R2.
        cv_cor[k * (iter - 1L) + i, j] <-
          1 - sum((y[membership == i] - p) ^ 2L) / sum((y[membership == i] - mean(y[membership == i])) ^ 2L)
        
        if(!silent) cat('CV iteration ', iter, ', lambda = ', lambda_list[j], ', partition ', i,
                        ', complete; ', 'R2 = ', round(cv_cor[k * (iter - 1L) + i, j], 3L), '\n', sep = '')
      }
    }
  }
  
  cv_cor
}

# Perform cross-validation across measurement replicates.
# This is useful when there are too few genotypes to partition them into train and test set.
cross_validation_across_replicates <- function(genotypes, replicates, site_combn, lambda_list,
                                               init_model = NULL, silent = 0L, epsilon = 1e-5) {
  
  # genotypes : numeric matrix
  #   Genotype matrix.
  # replicates : numeric matrix
  #   Each column of this matrix is a phenotype vector from each replicate.
  #   This matrix must be separated imported.
  # site_combn : list
  #   Site-combinations to model.
  # lambda_list : numeric vector
  #   Lasso penalty values to evaluate.
  # init_model : list (optional)
  #   Model to use as initial value for L-BFGS.
  # silent : 0 (default) or 1
  #   If 1, cross-validation progression is not reported.
  # epsilon : numeric
  #   Optimization precision; see the documentation for lbfgs::lbfgs for definition.
  
  # Only genotypes for which all replicate measurements are available are analyzed.
  exclude <- rowSums(is.finite(replicates)) < ncol(replicates)
  genotypes <- genotypes[!exclude, ]
  replicates <- replicates[!exclude, ]
  
  G <- construct_phenotype_operator(genotypes, site_combn)
  cv_cor <- matrix(NA_real_, ncol(replicates), length(lambda_list))
  
  for(i in 1L:ncol(replicates)) {
    
    # Mean phenotype excluding the test replicate set.
    y_minus_i <- apply(replicates[, -i, drop = FALSE], 1L, mean, na.rm = TRUE)
    
    for(j in seq_along(lambda_list)) {
      
      model <- infer_model(genotypes, y_minus_i, site_combn, lambda_list[j],
                           init_model = init_model, norm_e = FALSE, epsilon = epsilon)
      
      # Predicted phenotype.
      p <- apply_nonspecific_epistasis(as.vector(G %*% model$e), model$nonspec_param)
      
      # Out-of-sample R2.
      cv_cor[i, j] <-
        1 - sum((replicates[, i] - p) ^ 2L) / sum((replicates[, i] - mean(replicates[, i])) ^ 2L)
      
      if(!silent) cat('CV replicate ', i, ', lambda = ', lambda_list[j],
                      ' complete; R2 = ', round(cv_cor[i, j], 3L), '\n', sep = '')
    }
  }
  
  cv_cor 
}

# Identify the regularization strength corresponding to the sparsest among equally fit models.
identify_lambda <- function(lambda = NULL, cv = NULL, dir = NULL,
                            p_cutoff = 0.1, also_return_R2 = FALSE, show_plot = FALSE, ...) {
  
  # lambda : numeric vector
  #   Regularization strengths.
  # cv : numeric matrix
  #   Out-of-sample R2 values; each column corresponds to a regularization strength.
  # dir : string
  #   Directory to out-of-sample R2 matrix file.
  # p_cutoff : numeric
  #   Out-of-sample R2 must differ by this threshold to be considered significantly different.
  # also_return_R2 : logical
  #   If FALSE, only returns the chosen regularization strength.
  #   Otherwise, also returns the corresponding out-of-sample R2 value.
  # show_plot : logical
  #   Whether a plot of out-of-sample R2 as a function of regularization strength should be displayed.
  # ... : optional parameters for `plot`.
  
  # Importing file if given; column name should be the regularization strength.
  if(!is.null(dir)) {
    cv <- as.matrix(read.csv(dir, header = FALSE))
    lambda <- unname(cv[1L, ])
    cv <- unname(cv[-1L, ])
  }
  
  # Mean out-of-sample R2.
  cv_mean <- apply(cv, 2L, function(x) if (sum(is.finite(x)) > 2L) mean(x, na.rm = TRUE) else NA_real_)
  
  if(all(!is.finite(cv_mean))) if(also_return_R2) return(c(NA_real_, NA_real_)) else return(NA_real_)
  
  # Finding the maximal regularization strength among equally best-fit models.
  
  max_ind <- which.max(cv_mean)
  
  # Test of significance for difference in mean out-of-sample R2 to the best-fit model.
  p <- apply(cv, 2L, function(x) {
    
    if(sum(is.finite(x)) <= 2L) return(NA_real_)
    if(identical(x, cv[, max_ind])) return(1)
    if(var(x, na.rm = TRUE) < 1e-10 | var(cv[, max_ind], na.rm = TRUE) < 1e-10) 0 else t.test(x, cv[, max_ind])$p.value
  })
  
  best_lambda <- lambda[max(which(p > p_cutoff))]
  
  if(show_plot) {
    dotplot_with_sd(lambda, t(cv), plot_SEM = TRUE, ...)
    points(lambda, apply(t(cv), 1L, mean, na.rm = TRUE), type = 'l', ...)
    abline(v = best_lambda)
  }
  
  if(also_return_R2)
    c(best_lambda, cv_mean[max(which(p > p_cutoff))])
  else
    best_lambda
}

# Analyze the sparsity of genetic architecture.
analyze_sparsity <- function(genotypes, y, model, site_combn, k, lambda_list) {
  
  # genotypes : numeric matrix
  #   Genotype matrix.
  # y : numeric vector
  #   Phenotype vector.
  # model : list
  #   Model to detemrine the order of parameter inclusion.
  # site_combn : list
  #   List of site-combinations to model; must be consistent with the above model.
  # k : integer
  #   Number of partitions for cross-validation.
  # lambda_list : numeric vector
  #   Lasso penalty values to try for the cross-validation of the model to the subsample.
  #
  # Output is written real-time in a text file sparsity.txt.
  
  # Cross-validation function customized for analysis of sparsity.
  # This function takes as input a pre-calculated phenotype operator matrix, not a list of site-combinations.
  cross_validation_modified <- function(G, y, k, lambda_list, epsilon = 1e-5) {
    
    cv_cor <- matrix(NA_real_, k, length(lambda_list))
    
    membership <- sample(1L:nrow(G) %% k + 1L)
    
    for(i in 1L:k) { # i indexes the data partitions
      
      if(var(y[membership == i]) < 1e-10) next # No phenotypic variation to explain in the test set
      
      G_active <<- G[membership != i, ]
      y_active <<- y[membership != i]
      
      for(j in seq_along(lambda_list)) { # j indexes the regularization strengths
        
        # Initial parameter values.
        
        init_param <- rep(0, ncol(G_active))
        
        if(nonspec_type == 2L)
          init_param <- c(init_param, min(y_active), max(y_active) - min(y_active))
        
        res <- lbfgs(objective, gradient, init_param, orthantwise_c = lambda_list[j], orthantwise_start = 1L,
                     invisible = 1L, epsilon = epsilon)$par
        
        e <- res[1L:ncol(G_active)]
        
        if(nonspec_type == 1L)
          nonspec_param <- NULL
        else
          nonspec_param <- res[(ncol(G_active) + 1L):length(res)]
        
        p <- apply_nonspecific_epistasis(as.vector(G[membership == i, ] %*% e), nonspec_param)
        
        # Out-of-sample R2.
        cv_cor[i, j] <-
          1 - sum((y[membership == i] - p) ^ 2L) / sum((y[membership == i] - mean(y[membership == i])) ^ 2L)
      }
    }
    
    cv_cor
  }
  
  # Determining the order of parameter expansion.
  
  # Variance due to each effect.
  v <- unlist(lapply(parse_effects_by_site_combn(model$e, site_combn)[-1L], function(x) x ^ 2L / length(x)))
  
  # Order of effects to include; intercept is always included.
  ord <- c(1L, order(v, decreasing = TRUE) + 1L)
  
  # Full phenotype operator matrix.
  G <- construct_phenotype_operator(genotypes, site_combn)
  
  perform_lasso <- FALSE # Lasso penalty is activated when the number of parameters is sufficiently large.
  
  # Intercept-only model explains no variance.
  write(paste(rep(0, k), collapse = ','), 'sparsity.txt')
  df <- 1L # Degree of freedom
  
  for(i in 2L:min(250L, length(model$e))) {
    
    cat('Parameter ', i, ';', sep = '')
    
    # Calculating the degree of freedom of the model.
    rank_i <- rankMatrix(G[, ord[1L:i]], method = 'qr')
    
    if(!(rank_i > df)) { # Last parameter added did not increase the degree of freedom
      
      cat(' skipped because it is redundant\n')
      next
      
    } else {
      
      df <- df + 1L
      cat(' d.f. = ', df, ';', sep = '')
    }
    
    if(perform_lasso) {
      
      cv_R2 <- cross_validation_modified(G[, ord[1L:i]], y, k, lambda_list, epsilon = 1e-3)
      cv_R2_mean <- apply(cv_R2, 2L, mean)
      cat(' R2 =', max(cv_R2_mean), '\n')
      write(paste(signif(cv_R2[, which.max(cv_R2_mean)], 3L), collapse = ','), 'sparsity.txt', append = TRUE)
      
    } else {
      
      if(i %% 20L == 0L) { # Perform lasso regression with each 20 additional parameters.
        
        cv_R2 <- cross_validation_modified(G[, ord[1L:i]], y, k, c(0, lambda_list), epsilon = 1e-3)
        cv_R2_mean <- apply(cv_R2, 2L, mean)
        cat(' R2 =', max(cv_R2_mean), '\n')
        write(paste(signif(cv_R2[, which.max(cv_R2_mean)], 3L), collapse = ','), 'sparsity.txt', append = TRUE)
        
        # If better model fit is observed for a nonzero lambda, begin performing lasso regression.
        if(max(cv_R2_mean) - cv_R2_mean[1L] > 0.002) {
          
          cat('  Performing lasso regression from now on\n')
          perform_lasso <- TRUE
        }
        
      } else {
        
        cv_R2 <- cross_validation_modified(G[, ord[1L:i]], y, k, 0, epsilon = 1e-3)
        cat(' R2 =', mean(cv_R2), '\n')
        write(paste(signif(cv_R2, 3L), collapse = ','), 'sparsity.txt', append = TRUE)
      }
    }
  }
}

# Analyze the sparsity of genetic architecture.
# Performs cross-validation across measurement replicates to calculate the out-of-sample R2.
# This is useful when there are too few genotypes to partition them into train and test set.
analyze_sparsity_replicates <- function(genotypes, replicates, model, site_combn, lambda_list) {
  
  # genotypes : numeric matrix
  #   Genotype matrix.
  # replicates : numeric matrix
  #   Phenotype replicate measurements.
  # model : list
  #   Model to detemrine the order of parameter inclusion.
  # site_combn : list
  #   List of site-combinations to model; must be consistent with the above model.
  # lambda_list : numeric vector
  #   Lasso penalty values to try for the cross-validation of the model to the subsample.
  #
  # Output is written real-time in a text file sparsity.txt.
  
  # Cross-validation function customized for analysis of sparsity.
  # This function takes as input a pre-calculated phenotype operator matrix, not a list of site-combinations.
  # Cross-validation is performed across replicate measurements.
  cross_validation_modified <- function(G, replicates, lambda_list, epsilon = 1e-5) {
    
    # Only genotypes for which all replicate measurements are available are analyzed.
    exclude <- rowSums(is.finite(replicates)) < ncol(replicates)
    genotypes <- genotypes[!exclude, ]
    replicates <- replicates[!exclude, ]
    G <- G[!exclude, ]
    
    cv_cor <- matrix(NA_real_, ncol(replicates), length(lambda_list))
    
    G_active <<- G
    
    for(i in 1L:ncol(replicates)) {
      
      # Mean phenotype excluding the test replicate set.
      y_active <<- apply(replicates[, -i, drop = FALSE], 1L, mean, na.rm = TRUE)
      
      for(j in seq_along(lambda_list)) {
        
        # Initial parameter values.
        
        init_param <- rep(0, ncol(G_active))
        
        if(nonspec_type == 2L)
          init_param <- c(init_param, min(y_active), max(y_active) - min(y_active))
        
        res <- lbfgs(objective, gradient, init_param, orthantwise_c = lambda_list[j], orthantwise_start = 1L,
                     invisible = 1L, epsilon = epsilon)$par
        
        e <- res[1L:ncol(G_active)]
        
        if(nonspec_type == 1L)
          nonspec_param <- NULL
        else
          nonspec_param <- res[(ncol(G_active) + 1L):length(res)]
        
        # Predicted phenotype.
        p <- apply_nonspecific_epistasis(as.vector(G %*% e), nonspec_param)
        
        # Out-of-sample R2.
        cv_cor[i, j] <-
          1 - sum((replicates[, i] - p) ^ 2L) / sum((replicates[, i] - mean(replicates[, i])) ^ 2L)
      }
    }
    
    cv_cor
  }
  
  # Determining the order of parameter expansion.
  
  # Variance due to each effect.
  v <- unlist(lapply(parse_effects_by_site_combn(model$e, site_combn)[-1L], function(x) x ^ 2L / length(x)))
  
  # Order of effects to include; intercept is always included.
  ord <- c(1L, order(v, decreasing = TRUE) + 1L)
  
  # Full phenotype operator matrix.
  G <- construct_phenotype_operator(genotypes, site_combn)
  
  perform_lasso <- FALSE # Lasso penalty is activated when the number of parameters is sufficiently large.
  
  # Intercept-only model explains no variance.
  write(paste(rep(0, ncol(replicates)), collapse = ','), 'sparsity.txt')
  df <- 1L # Degree of freedom
  
  for(i in 2L:min(250L, length(model$e))) {
    
    cat('Parameter ', i, ';', sep = '')
    
    # Calculating the degree of freedom of the model.
    rank_i <- rankMatrix(G[, ord[1L:i]], method = 'qr')
    
    if(!(rank_i > df)) { # Last parameter added did not increase the degree of freedom
      
      cat(' skipped because it is redundant\n')
      next
      
    } else {
      
      df <- df + 1L
      cat(' d.f. = ', df, ';', sep = '')
    }
    
    if(perform_lasso) {
      
      cv_R2 <- cross_validation_modified(G[, ord[1L:i]], replicates, lambda_list, epsilon = 1e-3)
      cv_R2_mean <- apply(cv_R2, 2L, mean)
      cat(' R2 =', max(cv_R2_mean), '\n')
      write(paste(signif(cv_R2[, which.max(cv_R2_mean)], 3L), collapse = ','), 'sparsity.txt', append = TRUE)
      
    } else {
      
      if(i %% 20L == 0L) { # Perform lasso regression with each 20 additional parameters.
        
        cv_R2 <- cross_validation_modified(G[, ord[1L:i]], replicates, c(0, lambda_list), epsilon = 1e-3)
        cv_R2_mean <- apply(cv_R2, 2L, mean)
        cat(' R2 =', max(cv_R2_mean), '\n')
        write(paste(signif(cv_R2[, which.max(cv_R2_mean)], 3L), collapse = ','), 'sparsity.txt', append = TRUE)
        
        # If better model fit is observed for a nonzero lambda, begin performing lasso regression.
        if(max(cv_R2_mean) - cv_R2_mean[1L] > 0.002) {
          
          cat('  Performing lasso regression from now on\n')
          perform_lasso <- TRUE
        }
        
      } else {
        
        cv_R2 <- cross_validation_modified(G[, ord[1L:i]], replicates, 0, epsilon = 1e-3)
        cat(' R2 =', mean(cv_R2), '\n')
        write(paste(signif(cv_R2, 3L), collapse = ','), 'sparsity.txt', append = TRUE)
      }
    }
  }
}

# Infer model from a random subsample of specified size and predict the phenotype of all other genotypes.
predict_from_subsample <- function(genotypes, y, order_list, frac_list, niter, lambda_list, CV_k) {
  
  # genotypes : numeric matrix
  #   Genotype matrix.
  # y : numeric vector
  #   Phenotype vector.
  # order_list : integer vector
  #   Models to fit on the subsample; model that best describes the subsample is used for prediction.
  # frac_list : numeric vector
  #   Sample size values to evaluate.
  #   This is interpreted as the fraction of all possible genotypes, not that of experimentally characterized genotypes.
  # niter : integer
  #   Number of samples for each sample size.
  # lambda_list : numeric vector
  #   Lasso penalty values to try for the cross-validation of the model to the subsample.
  # CV_k : integer
  #   Number of partitions for cross-validation.

  prediction <- matrix(NA_real_, length(frac_list), niter, dimnames = list(frac_list))
  
  site_combn_list <- lapply(order_list, generate_site_combn)
  G_list <- lapply(order_list, function(ord) construct_phenotype_operator(genotypes, generate_site_combn(ord)))
  
  for(i in seq_along(frac_list)) { # i indexes sample sizes
    
    cat('Prediction with sample size', frac_list[i], '\n')
    
    for(j in 1L:niter) { # j indexes iterations for each sample size
      
      cat(' Sample ', j, '; ', sep = '')
      
      subsample <- sample(nrow(genotypes), frac_list[i] * prod(q_list))
      
      if(var(y[-subsample]) < 1e-10) next # Skips when there is no phenotypic variation to predict
      
      # Determining the best-fit model by cross-validation.
      # Cross-validation is repeated until the standard error of best out-of-sample R2 is less than 0.05.
      
      cv <- matrix(NA_real_, length(order_list), 2L)
      
      for(k in seq_along(order_list)) {
        
        cv_ord <- matrix(nrow = 0L, ncol = length(lambda_list))
        
        cv_ord <- rbind(cv_ord, cross_validation(genotypes[subsample, ], y[subsample], site_combn_list[[k]],
                                                 CV_k, lambda_list, silent = 1L, epsilon = 1e-2))
        
        # Out-of-sample R2 values for the best lasso penalty.
        best_r2 <- cv_ord[, lambda_list == identify_lambda(lambda_list, cv_ord, p_cutoff = 0.5)]
        
        # Repeat cross-validation if standard error of the best out-of-sample R2 is greater than 0.05.
        while(all(!is.finite(best_r2)) | sd(best_r2, na.rm = TRUE) / sqrt(sum(is.finite(best_r2))) > 0.05) {
          
          if(nrow(cv_ord) > 100L) break # Capped at 100 iterations
          
          cv_ord <- rbind(cv_ord, cross_validation(genotypes[subsample, ], y[subsample], site_combn_list[[k]],
                                                   CV_k, lambda_list, silent = 1L, epsilon = 1e-2))
          best_r2 <- cv_ord[, lambda_list == identify_lambda(lambda_list, cv_ord, p_cutoff = 0.5)]
        }
        
        cv[k, ] <- identify_lambda(lambda_list, cv_ord, p_cutoff = 0.5, also_return_R2 = TRUE)
      }
      
      # Prediction.
      
      if(all(!is.finite(cv))) { # If every iteration of cross-validation failed
        
        # This happens when there is no phenotypic variation within the subsample.
        p <- rep(mean(y[subsample]), length(y))
      
      } else {
        
        k <- which.max(cv[, 2L]) # Model to use.
        
        model <- infer_model(genotypes[subsample, ], y[subsample], site_combn_list[[k]],
                             cv[k, 1L], norm_e = FALSE, epsilon = 1e-3)
        
        p <- apply_nonspecific_epistasis(as.vector(G_list[[k]] %*% model$e), model$nonspec_param)
      }
      
      # Out-of-sample R2.
      prediction[i, j] <-
        1 - sum((y[-subsample] - p[-subsample]) ^ 2L) / sum((y[-subsample] - mean(y[-subsample])) ^ 2L)
      
      cat('best-fit model order: ', order_list[[k]], ', prediction R2: ', prediction[i, j], '\n', sep = '')
    }
  }
  
  prediction
}

# Infer model from a random subsample of specified size and predict the phenotype of all other genotypes.
# Performs cross-validation across measurement replicates to calculate the out-of-sample R2.
# This is useful when there are too few genotypes to partition them into train and test set.
predict_from_subsample_replicates <- function(genotypes, y, replicates, order_list, frac_list,
                                              niter, lambda_list) {
  
  # genotypes : numeric matrix
  #   Genotype matrix.
  # y : numeric vector
  #   Phenotype vector.
  # replicates : numeric matrix
  #   Measurement replicates.
  # order_list : integer vector
  #   Models to fit on the subsample; model that best describes the subsample is used for prediction.
  # frac_list : numeric vector
  #   Sample size values to evaluate.
  #   This is interpreted as the fraction of all possible genotypes, not that of experimentally characterized genotypes.
  # niter : integer
  #   Number of samples for each sample size.
  # lambda_list : numeric vector
  #   Lasso penalty values to try for the cross-validation of the model to the subsample.
  # CV_k : integer
  #   Number of partitions for cross-validation.
  
  prediction <- matrix(NA_real_, length(frac_list), niter, dimnames = list(frac_list))
  
  site_combn_list <- lapply(order_list, generate_site_combn)
  G_list <- lapply(order_list, function(ord) construct_phenotype_operator(genotypes, generate_site_combn(ord)))
  
  for(i in seq_along(frac_list)) { # i indexes sample sizes
    
    cat('Prediction with sample size', frac_list[i], '\n')
    
    for(j in 1L:niter) { # j indexes iterations for each sample size
      
      cat(' Sample ', j, '; ', sep = '')
      
      subsample <- sample(nrow(genotypes), frac_list[i] * prod(q_list))
      
      if(var(y[-subsample]) < 1e-10) next # Skips when there is no phenotypic variation to predict
      
      # Determining the best-fit model by cross-validation.
      
      cv <- matrix(NA_real_, length(order_list), 2L)
      
      for(k in seq_along(order_list)) {
        
        cv_ord <- cross_validation_across_replicates(genotypes[subsample, ], replicates[subsample, ],
                                                     site_combn_list[[k]], lambda_list, silent = 1L, epsilon = 1e-2)
        
        cv[k, ] <- identify_lambda(lambda = lambda_list, cv = cv_ord, p_cutoff = 0.5, also_return_R2 = TRUE)
      }
      
      # Prediction.
      
      if(all(!is.finite(cv))) { # If every iteration of cross-validation failed
        
        # This happens when there is no phenotypic variation within the subsample.
        p <- rep(mean(y[subsample]), length(y))
        
      } else {
        
        k <- which.max(cv[, 2L]) # Model to use.
        
        model <- infer_model(genotypes[subsample, ], y[subsample], site_combn_list[[k]],
                             cv[k, 1L], norm_e = FALSE, epsilon = 1e-3)
        
        p <- apply_nonspecific_epistasis(as.vector(G_list[[k]] %*% model$e), model$nonspec_param)
      }
      
      # Out-of-sample R2.
      prediction[i, j] <-
        1 - sum((y[-subsample] - p[-subsample]) ^ 2L) / sum((y[-subsample] - mean(y[-subsample])) ^ 2L)
      
      cat('best-fit model order: ', order_list[[k]], ', prediction R2: ', prediction[i, j], '\n', sep = '')
    }
  }
  
  prediction
}

# Export model.
export_model <- function(dir, n, q_list, site_combn, nonspec_type, nonspec_param, e) {
  
  f <- file(dir)
  open(f, 'w')
  
  writeLines('n', f)
  write(n, f)
  write('q_list', f)
  write(paste(q_list, collapse = ' '), f)
  write('site_combn', f)
  for(x in site_combn[-1L]) write(paste(x, collapse = ' '), f) # Intercept not printed
  write('nonspec_type', f)
  write(nonspec_type, f)
  write('nonspec_param', f)
  write(paste(signif(nonspec_param, 5L), collapse = ' '), f)
  write('e', f)
  write(as.character(signif(e, 5L)), f)
  
  close(f)
}

# Import model.
import_model <- function(dir) {
  
  # dir : character
  #   Directory to a text file specifying a model.
  #
  # The following global parameters are overwritten and must be compatible with other data:
  #   `n`, `q_list`, `site_combn`, `nonspec_type`
  
  file <- readLines(dir)
  
  # Number of sites.
  n <<- as.integer(file[which(file == 'n') + 1L])
  
  # Number of states.
  q_list <<- as.integer(strsplit(file[which(file == 'q_list') + 1L], split = ' ')[[1L]])
  
  # Site-combinations.
  site_combn <- strsplit(file[(which(file == 'site_combn') + 1L):(which(file == 'nonspec_type') - 1L)], split = ' ')
  site_combn <<- c(list(integer(0)), lapply(site_combn, function(x) as.integer(x)))
  
  # Type of global nonlinearity.
  nonspec_type <<- as.integer(file[which(file == 'nonspec_type') + 1L])
  
  # Global nonlinearity parameters
  nonspec_param <- as.numeric(strsplit(file[which(file == 'nonspec_param') + 1L], split = ' ')[[1L]])
  
  # Specific effects.
  e <- as.numeric(file[(which(file == 'e') + 1L):length(file)])
  
  # Returning model.
  list(nonspec_param = nonspec_param, e = e)
}

# Generate a dotplot with standard deviation or standard error.
dotplot_with_sd <- function(X, Y, plot_SEM = FALSE, ...) {
  
  # X, Y : numeric vector, matrix, or list
  #   If a matrix, rows correspond to data points and columns to replicate values.
  #   If a vector, Sd is not plotted.
  # plot_SEM : boolean
  #   If `TRUE`, plots standard error of the mean; otherwise plots standard deviation.
  # ... : optional additional parameters
  #   Passed to `plot` for managing graphic parameters.
  
  if(is.vector(X) & !is.list(X)) X <- vapply(X, as.list, list(1))
  if(is.vector(Y) & !is.list(Y)) Y <- vapply(Y, as.list, list(1))
  if(is.matrix(X)) X <- split(X, rep(1L:nrow(X), ncol(X)))
  if(is.matrix(Y)) Y <- split(Y, rep(1L:nrow(Y), ncol(Y)))
  
  X_mean <- vapply(X, mean, numeric(1L), na.rm = TRUE)
  Y_mean <- vapply(Y, mean, numeric(1L), na.rm = TRUE)
  
  if(plot_SEM) {
    
    X_range <- cbind(X_mean - vapply(X, function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))), 0),
                     X_mean + vapply(X, function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))), 0))
    
    Y_range <- cbind(Y_mean - vapply(Y, function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))), 0),
                     Y_mean + vapply(Y, function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))), 0))
    
  } else {
  
    X_range <- cbind(X_mean - vapply(X, sd, numeric(1L), na.rm = TRUE), X_mean + vapply(X, sd, numeric(1L), na.rm = TRUE))
    Y_range <- cbind(Y_mean - vapply(Y, sd, numeric(1L), na.rm = TRUE), Y_mean + vapply(Y, sd, numeric(1L), na.rm = TRUE)) 
  }
  
  plot(X_mean, Y_mean, ...)
  
  for(k in 1L:length(X)) {
    
    if(all(is.finite(Y_range[k, ]))) arrows(x0 = X_mean[k], y0 = Y_range[k, 1L], x1 = X_mean[k], y1 = Y_range[k, 2L], angle = 90, code = 3, length = 0.1)
    if(all(is.finite(X_range[k, ]))) arrows(x0 = X_range[k, 1L], y0 = Y_mean[k], x1 = X_range[k, 2L], y1 = Y_mean[k], angle = 90, code = 3, length = 0.1)
  }
}

# Generate a dotplot with confidence interval.
dotplot_with_CI <- function(X, Y, CI = 0.95, plot_line = FALSE, ...) {
  
  # X, Y : numeric vector, matrix, or list
  #   If a matrix, rows correspond to data points and columns to replicate values.
  #   If a vector, Sd is not plotted.
  # CI : numeric
  #   Confidence interval width.
  # ... : optional additional parameters
  #   Passed to `plot` for managing graphic parameters.
  
  if(is.vector(X) & !is.list(X)) X <- vapply(X, as.list, list(1))
  if(is.vector(Y) & !is.list(Y)) Y <- vapply(Y, as.list, list(1))
  if(is.matrix(X)) X <- split(X, rep(1L:nrow(X), ncol(X)))
  if(is.matrix(Y)) Y <- split(Y, rep(1L:nrow(Y), ncol(Y)))
  
  X_median <- vapply(X, median, numeric(1L), na.rm = TRUE)
  Y_median <- vapply(Y, median, numeric(1L), na.rm = TRUE)
  
  X_range <- cbind(vapply(X, function(x) quantile(x, 0.5 - CI / 2), 0),
                   vapply(X, function(x) quantile(x, 0.5 + CI / 2), 0))
  
  Y_range <- cbind(vapply(Y, function(x) quantile(x, 0.5 - CI / 2), 0),
                   vapply(Y, function(x) quantile(x, 0.5 + CI / 2), 0))
  
  plot(X_median, Y_median, ...)
  
  if(plot_line) points(X_median, Y_median, type = 'l')
  
  for(k in 1L:length(X)) {
    
    if(all(is.finite(Y_range[k, ]))) arrows(x0 = X_median[k], y0 = Y_range[k, 1L], x1 = X_median[k], y1 = Y_range[k, 2L], angle = 90, code = 3, length = 0.1)
    if(all(is.finite(X_range[k, ]))) arrows(x0 = X_range[k, 1L], y0 = Y_median[k], x1 = X_range[k, 2L], y1 = Y_median[k], angle = 90, code = 3, length = 0.1)
  }
}

# Function for linearly interpolating a value for x given a target value of y.
interpolate <- function(x, y, y_value) {
  
  if(y_value > max(y)) return(NA_real_)
  if(y_value < min(y)) return(NA_real_)
  
  x_upper <- min(x[y >= y_value])
  x_lower <- max(x[y < y_value])
  y_upper <- min(y[y >= y_value])
  y_lower <- max(y[y < y_value])
  
  x_lower + ((x_upper - x_lower) / (y_upper - y_lower)) * (y_value - y_lower)
}

## For simulating genetic architectures:

# Enumerate all genotypes
enumerate_genotypes <- function(n, s) {
  
  # Generate all permutations of elements in a list of vectors, the right-most vector fastest-changing.
  permutations <- function(vecs) {
    
    # If all vectors have only one element.
    if(all(vapply(vecs, length, 0L) == 1L)) return(matrix(unlist(vecs), nrow = 1L))
    
    if(length(vecs) == 1L) {
      matrix(vecs[[1L]], ncol = 1L)
    } else {
      sp <- permutations(vecs[-1L])
      cbind(rep(vecs[[1L]], each = nrow(sp)), do.call(rbind, rep(list(sp), length(vecs[[1L]]))))
    }
  }
  
  permutations(rep(list(1L:s), n))
}

# Simulate a genetic architecture.
simulate_GP_map <- function(n, q_list, var_spec, nose_var, distr = 'normal') {
  
  # n, s : int
  #   Number of sites and states.
  # var_spec : numeric
  #   Variance spectrum, or the absolute genetic variance due to each epistatic order.
  # noise_var : numeric
  #   Phenotypic variance due to noise.
  # distr : character ('normal', 'uniform')
  #   Distribution from which non-zero effects are drawn.
  
  genotypes <- enumerate_genotypes(n, s)
  site_combn <- generate_site_combn(n)
  P <- construct_phenotype_operator(genotypes, site_combn)
  
  # Simulating terms without regards to variance spectrum.
  e <- unlist(lapply(0L:n, function(model_ord) {
    
    if(model_ord == 0L) return(0)
    
    n_terms <- choose(n, model_ord) * s ^ model_ord # Number of terms of the order
    
    if(var_spec[model_ord] == 0) return(rep(0, n_terms))
    
    if(distr == 'normal')
      e_model_ord <- rnorm(n_terms, 0, 1)
    else if(distr == 'uniform')
      e_model_ord <- runif(n_terms, -1, 1)
    
    e_model_ord
  }))
  
  e <- normalize_e(e, site_combn, genotypes) # Enforcing the zero-mean constraint
  e_parsed <- parse_effects_by_site_combn(e, site_combn) # Parsing by site-combination
  
  # Adjusting variance spectrum.
  
  # Variance due to each site-combination.
  v <- c(0, vapply(e_parsed[-1L], function(x) mean(x ^ 2L), numeric(1L)))
  
  # Order of each site-combination.
  ord <- unlist(lapply(0L:n, function(i) rep(i, choose(n, i)))) # Order of each effect set
  
  # Variance due to each order of terms.
  v_ord <- vapply(1L:n, function(i) sum(v[ord == i]), numeric(1L))
  
  for(i in seq_along(e_parsed)) {
    
    if(ord[i] == 0L) { # Intercept is set at 0
      
      e_parsed[[i]] <- 0
      
    } else if(var_spec[ord[i]] == 0) {
      
      e_parsed[[i]] <- e_parsed[[i]] * 0
      
    } else {
      
      e_parsed[[i]] <- e_parsed[[i]] * sqrt(var_spec[ord[i]] / v_ord[ord[i]])
    }
  }
  
  e <- unlist(e_parsed)
  
  # Adding measurement noise.
  y <- as.vector(P %*% e) + rnorm(s ^ n, 0, sqrt(noise_var))
  
  list(genotypes = genotypes, y = y, e = e) 
}


