
# This script implements the functions required for joint inference of
# reference-free effects and nonspecific epistasis. See the example script RFA
# for more information.
#
# Written in R version 4.3.1.
# Packages required and the version used:

library(abind) # 1.4-5
library(lbfgs) # 1.2.1.2
library(Matrix) # 1.6-0


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
  # States are renumbered as 1, 2, 3, ...
  
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
  
  # Recoding states if they are not 1, 2, 3, ... consecutively.
  
  # States in each site.
  states_list <- apply(genotypes, 2, function(x) sort(unique(x)), simplify = FALSE)
  
  # Checking for unused states.
  reducible <- vapply(states_list, function(x) length(x) != max(x), logical(1))
  
  if(any(reducible)) {
    
    cat('States should be numbered 1, 2, 3, ... consecutively\n')
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
construct_phenotype_operator <- function(genotypes, site_combn) {
  
  # genotypes : int matrix
  #   Genotype matrix.
  # site_combn : list
  #   Site-combinations included in the model.
  #
  # The phenotype operator matrix is a sparse pattern matrix. The sparse matrix
  # implemention of the package Matrix is used.
  
  ord <- sapply(site_combn, length) # Order of each site-combination
  
  # Number of terms in each site-combination.
  n_terms <- as.integer(sapply(site_combn, function(x) prod(q_list[x])))
  
  # Cumulative number of terms.
  cum_terms <- cumsum(c(0L, n_terms[-length(n_terms)]))
  
  # Row indices of the elements to be set as 1 in the phenotype matrix.
  a <- rep(1L:nrow(genotypes), length(site_combn))
  
  # Column indices of the elements to be set as 1 in the phenotype matrix.
  b <- vapply(1L:length(site_combn), function(k) {
    
    if(ord[k] == 0L) { # Intercept
      
      rep(1L, nrow(genotypes))
      
    } else if(ord[k] == 1L) { # First-order effects
      
      cum_terms[k] + genotypes[, site_combn[[k]]]
      
    } else {
      
      # Within a site-combination, state-combinations are ordered by the last
      # site moving fastest.
      q_k <- q_list[site_combn[[k]]]
      cum_q <- as.integer(c(rev(cumprod(rev(q_k[-1L]))), 1L))
      
      cum_terms[k] + as.integer(genotypes[, site_combn[[k]]] %*% cum_q) + 
        (1L - sum(cum_q))
    }
    
  }, integer(nrow(genotypes)))
  
  sparseMatrix(a, b, dims = c(nrow(genotypes), sum(n_terms)))
}

# Sort reference-free effects by site-combination.
parse_effects_by_site_combn <- function(site_combn, e, as_array = FALSE) {
  
  # site_combn : list
  #   List of site-combinations being analyzed.
  # e : numeric vector
  #   Reference-free effects for the site-combinations.
  # as_array : logical
  #   If TRUE, state-combinations within each site-combination are arranged in
  #   an array. Otherwise, they are arranged in a vector.
  
  # Identifying reference-free effects for each site-combination.
  
  ord <- sapply(site_combn, length) # Order of each site-combination
  
  # Number of terms in each site-combination.
  n_terms <- sapply(site_combn, function(x) prod(q_list[x]))
  
  # Index of the last term of each site-combination.
  pos <- cumsum(c(0L, n_terms))
  
  # Sorting the effects for each site-combination into an element in a list.
  e_parsed <- lapply(1L:length(site_combn), function(i) e[(pos[i] + 1L):pos[i + 1L]])
  
  if(as_array) { 
    
    # Within a site-combination, state-combinations are ordered by the last
    # site moving fastest.
    
    for(i in which(ord > 1L)) {
      
      sites <- site_combn[[i]]
      qs <- q_list[sites] # Number of states in each site
      E <- array(e_parsed[[i]], rev(qs))
      e_parsed[[i]] <- aperm(E, length(sites):1L)
    }
  }
  
  e_parsed
}

# Calculate the fraction of variance due to each reference-free effect.
partition_variance <- function(site_combn, e, as_array = TRUE) {
  
  # site_combn : list
  #   Site-combinations modeled.
  # e : numeric vector
  #   Reference-free effects.
  # as_array : logical
  #   If TRUE, state-combinations within each site-combination are arranged in
  #   an array. Otherwise, they are arranged in a vector.
  
  e_parsed <- parse_effects_by_site_combn(site_combn, e, as_array)
  v <- lapply(e_parsed[-1L], function(x) x ^ 2L / length(x))
  lapply(v, function(x) x / sum(unlist(v)))
}

# Enforce the zero-mean property.
normalize_effects <- function(e, site_combn, genotypes) {
  
  # e : numeric vector
  #   Reference-free effects to be normalized.
  # site_combn : list
  #   Site-combinations corresponding to e.
  # genotypes : int matrix
  #   Genotypes from which e was inferred.
  #
  # The algorithm used for normalization is explained in the supplementary
  # text, Appendix 6, of the paper.
  
  # Normalize all k-th-order effects without changing the other effects.
  # But if k is 1, the intercept is also normalized.
  normalize_effects_k <- function(e, site_combn, k) {
    
    ord <- sapply(site_combn, length) # Order of each site-combination
    
    if(sum(ord == k) == 0L) return (e) # No k-th-order effects provided
    
    e_parsed <- parse_effects_by_site_combn(site_combn, e)
    
    if(k == 1L) {
      
      # Normalizing the intercept.
      e_parsed[[1L]] <- e_parsed[[1L]] + sum(sapply(which(ord == 1L),
                                                    function(i) mean(e_parsed[[i]])))
      
      # Enforcing the constraint on all first-order effects.
      for(i in which(ord == 1L))
        e_parsed[[i]] <- e_parsed[[i]] - mean(e_parsed[[i]])
      
    } else {
      
      # Implementing the normalization formula in the supplementary text.
      
      for(i in which(ord == k)) {
        
        sites <- site_combn[[i]]
        qi <- q_list[sites]
        rev_qi <- rev(qi)
        
        E <- array(e_parsed[[i]], rev_qi)
        D <- E
        sgn <- 1
        
        for(low_ord in (k - 1L):1L) {
          
          correction_matrices <- lapply(combn(k, low_ord, simplify = FALSE),
                                        function(margin) {
            
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
  
  ord <- sapply(site_combn, length) # Order of each site combination
  
  if(max(ord) == 1L) return(normalize_effects_k(e, site_combn, 1L))
  
  X <- construct_phenotype_operator(genotypes, site_combn)
  s <- as.vector(X %*% e) # Genetic score
  
  # Number of terms in each site-combination.
  n_terms <- sapply(site_combn, function(x) prod(q_list[x]))
  
  # Order of each effect
  e_ord <- unlist(lapply(1L:length(site_combn), function(i) rep(ord[i], n_terms[i])))
  
  for(k in max(ord):2L) {
    
    e <- normalize_effects_k(e, site_combn, k) # Normalizing order k
    
    # Altering the lower-order effects to keep the predicted phenotypes same.
    
    X_l_k <- X[, e_ord < k]
    X_ge_k <- X[, e_ord >= k]
    e_l_k <- e[e_ord < k]
    e_ge_k <- e[e_ord >= k]
    
    # Objective function for correcting the lower-order effects.
    objective <- function(e_l_k) {
      
      z0 <<- s - as.vector(X_l_k %*% e_l_k + X_ge_k %*% e_ge_k)
      sum(z0 ^ 2L)
    }
    
    # Gradient function for correcting the lower-order effects.
    gradient <- function(e_l_k) as.vector(-2 * z0 %*% X_l_k)
    
    e[e_ord < k] <- lbfgs(objective, gradient, e_l_k, invisible = 1)$par
  }

  normalize_effects_k(e, site_combn, 1L)
}

# Objective function to minimize for the joint inference of reference-free
# effects and nonspecific epistasis.
objective <- function(param) {
  
  # param : numeric vector
  #   Parameters to evaluate; reference-free effects followed by the nonspecific
  #   epistasis parameters.
  #
  # This function computes the sum of the squared residuals. It relies on a
  # globally stored phenotype vector and phenotype operator (`y_active` and
  # `X_active`), which is created by the optimizer that calls this function.
  # The lasso penalty is calculated by the optimizer.
  
  if(nonspec_type == 1L) {
    
    e <- param
    
    # Residuals; stored as a global variable to be used in the gradient
    # calculation by the function `gradient`.
    z0 <<- y_active - as.vector(X_active %*% e)
    
    sum(z0 ^ 2L)
    
  } else if(nonspec_type == 2L) {
    
    nonspec_param <- param[(length(param) - 1L):length(param)]
    e <- param[1L:(length(param) - 2L)]
    
    # Partial residuals; stored as a global variable to be used in the gradient
    # calculation by the function `gradient`.
    z0 <<- nonspec_param[2L] / (1 + exp(-as.vector(X_active %*% e)))
    
    sum((y_active - nonspec_param[1L] - z0) ^ 2L)
  }
}

# Gradient function to use in the joint inference of reference-free effects and
# the nonspecific epistasis parameters.
gradient <- function(param) {
  
  # See the description for the function `objective`.
  
  if(nonspec_type == 1L) {
    
    -2 * as.vector(z0 %*% X_active)
    
  } else if(nonspec_type == 2) {
    
    nonspec_param <- param[(length(param) - 1L):length(param)]
    e <- param[1L:(length(param) - 2L)]
    
    z1 <- -2 * (y_active - nonspec_param[1L] - z0)
    z2 <- z1 * z0 / nonspec_param[2L]
    z3 <- z2 * (nonspec_param[2L] - z0)
    c(as.vector(z3 %*% X_active), sum(z1), sum(z2))
  }
}

# Jointly infer reference-free effects and nonspecific epistasis.
infer_model <- function(genotypes, y, site_combn, lambda = 0, init_model = NULL,
                        norm_e = TRUE, silent = 0L, epsilon = 1e-5) {
  
  # genotypes : numeric matrix
  #   Genotype matrix.
  # y : numeric vector
  #   Phenotype vector.
  # site_combn : list
  #   Site-combinations to model.
  # lambda : numeric
  #   Lasso penalty.
  # init_model : list
  #   Initial model to use for optimization.
  # norm_e : logical
  #   If TRUE, the zero-mean constraint is enforced on the returned model.
  # silent : 0 or 1
  #   If 1, the console output of lbfgs::lbfgs is suppressed.
  # epsilon : numeric
  #   Optimization precision; see the documentation for lbfgs::lbfgs.
  #
  # Returns a list with the following elements.
  # nonspec_param : numeric vector
  #   The inferred nonspecific epistasis parameters.
  # e : numeric vector
  #   The inferred reference-free effects.
  # R2 : numeric
  #   Model fit.
  # residual : numeric vector
  #   Residual for each genotype.
  
  # Phenotye vector and operator for `objective` and `gradient`.
  X_active <<- construct_phenotype_operator(genotypes, site_combn)
  y_active <<- y
  
  if(var(y) < 1e-10) { # No phenotypic variation to model
    
    if(nonspec_type == 1L) {
      
      return(list(nonspec_param = NULL, e = c(mean(y), rep(0, ncol(X_active) - 1L)),
                  R2 = NA_real_, residual = rep(0, length(y))))
      
    } else if(nonspec_type == 2L) {
      
      return(list(nonspec_param = c(mean(y), 0), e = rep(0, ncol(X_active)),
                  R2 = NA_real_, residual = rep(0, length(y))))
    }
  }
  
  # Initial model for optimization.
  
  if(!is.null(init_model)) {
    
    init_param <- c(init_model$e, init_model$nonspec_param)
    
  } else {
    
    # Initial value for every reference-free effect is zero.
    init_param <- rep(0, ncol(X_active))
    
    # Initial value for the sigmoid parameters is the data range.
    if(nonspec_type == 2L)
      init_param <- c(init_param, min(y), max(y) - min(y))
  }
  
  # Nonspecific epistasis parameters are subject to regularization to avoid
  # non-identifiability for some cases.
  res <- lbfgs(objective, gradient, init_param, orthantwise_c = lambda,
               orthantwise_start = 1L, invisible = silent, epsilon = epsilon)$par
  
  if(norm_e)
    e <- normalize_effects(res[1:ncol(X_active)], site_combn, genotypes)
  else
    e <- res[1L:ncol(X_active)]
  
  if(nonspec_type == 1)
    nonspec_param <- NULL
  else
    nonspec_param <- res[(ncol(X_active) + 1):length(res)]
  
  # Pheontypes predicted by the model.
  p <- apply_nonspecific_epistasis(as.vector(X_active %*% e), nonspec_param)
  
  list(nonspec_param = nonspec_param,
       e = e,
       R2 = 1 - sum((y - p) ^ 2) / sum((y - mean(y)) ^ 2),
       residual = y - p)
}

# Transform genetic score by nonspecific epistasis.
apply_nonspecific_epistasis <- function(s, nonspec_param) {
  
  # s : numeric vector
  #   Genetic scores.
  # nonspec_param : numeric vector
  #   Parameters of the nonspecific epistasis function.
  #
  # The type of nonspecific epistasis is defined by the global variable
  # nonspec_type. Only two types of nonspecific epistasis is currently
  # supported: none or sigmoid.
  
  if(nonspec_type == 1L) # None
    s
  else if(nonspec_type == 2L) # Sigmoid
    nonspec_param[1L] + nonspec_param[2L] / (1 + exp(-s))
}

# Show nonspecific epistasis and the distribution of genetic score.
show_latent_space <- function(model, genotypes, y, site_combn) {
  
  # model : list
  #   Model to evaluate.
  # genotypes, y, site_combn
  #   Data used to infer the model.
  
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
export_model <- function(path, n, q_list, site_combn, nonspec_type, model) {
  
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
  nonspec_type <<- as.integer(file[which(file == 'nonspec_type') + 1L])
  
  # Nonspecific epistasis parameters.
  nonspec_param <- as.numeric(strsplit(file[which(file == 'nonspec_param') + 1L], split = ' ')[[1L]])
  
  # Reference-free effects.
  e <- as.numeric(file[(which(file == 'e') + 1L):length(file)])
  
  # Returning model.
  list(nonspec_param = nonspec_param, e = e)
}

# Perform k-fold cross-validation.
cross_validation <- function(genotypes, y, site_combn, k, lambda_list, niter,
                             silent = 0, epsilon = 1e-5) {
  
  # genotypes, y, site_combn, epsilon : As defined in `infer_model`
  # k : int
  #   Number of genotype partitions for cross-validation.
  # lambda_list : numeric vector
  #   Lasso penalty values to evaluate.
  # niter : int
  #   Number of times the k-fold cross-validation should be performed.
  # silent : 0 or 1
  #   If 1, cross-validation progress is not reported.
  
  X <- construct_phenotype_operator(genotypes, site_combn)
  cvR2 <- matrix(NA_real_, k * niter, length(lambda_list))
  
  for(iter in 1L:niter) {
    
    # Random partition of genotypes into k categories.
    membership <- sample(1L:nrow(genotypes) %% k + 1L)
    
    for(j in seq_along(lambda_list)) { # j indexes lasso penalties
      
      for(i in 1L:k) { # i indexes genotype partitions
        
        # If the test set has no phenotypic variance to explain.
        if(var(y[membership == i]) < 1e-10) next
        
        model <- infer_model(genotypes[membership != i, ], y[membership != i],
                             site_combn, lambda_list[j], norm_e = FALSE,
                             silent = 1L, epsilon = epsilon)
        
        # Predicted phenotype for test set.
        p <- apply_nonspecific_epistasis(as.vector(X[membership == i, ] %*% model$e),
                                         model$nonspec_param)
        
        # Out-of-sample R2.
        cvR2[k * (iter - 1L) + i, j] <-
          1 - sum((y[membership == i] - p) ^ 2L) / sum((y[membership == i] - mean(y[membership == i])) ^ 2L)
        
        if(!silent) cat('CV iteration ', iter, ', lambda = ', lambda_list[j],
                        ', partition ', i, '; ', 'R2 = ',
                        round(cvR2[k * (iter - 1L) + i, j], 3L), '\n', sep = '')
      }
    }
  }
  
  cvR2
}

# Perform cross-validation across measurement replicates. This is useful when
# there are too few genotypes to partition them into training and test sets.
cv_across_replicates <- function(genotypes, replicates, site_combn, lambda_list,
                                 silent = 0L, epsilon = 1e-5) {
  
  # genotypes, site_combn, lambda_list, silent, epsilon :
  #   As defined in `infer_model` and `cross_validation`.
  # replicates : numeric matrix
  #   Replicate measurements; rows correspond to genotypes, columns to replicates.
  
  # Excluding genotypes missing any replicate measurement.
  exclude <- rowSums(is.finite(replicates)) < ncol(replicates)
  genotypes <- genotypes[!exclude, ]
  replicates <- replicates[!exclude, ]
  
  X <- construct_phenotype_operator(genotypes, site_combn)
  cvR2 <- matrix(NA_real_, ncol(replicates), length(lambda_list))
  
  for(i in 1L:ncol(replicates)) { # i indexes replicates
    
    # Mean phenotype excluding the test replicate set.
    y_minus_i <- apply(replicates[, -i, drop = FALSE], 1L, mean, na.rm = TRUE)
    
    for(j in seq_along(lambda_list)) { # j indexes lasso penalties
      
      model <- infer_model(genotypes, y_minus_i, site_combn, lambda_list[j],
                           norm_e = FALSE, silent = 0L, epsilon = epsilon)
      
      # Predicted phenotype.
      p <- apply_nonspecific_epistasis(as.vector(X %*% model$e), model$nonspec_param)
      
      # Out-of-sample R2.
      cvR2[i, j] <-
        1 - sum((replicates[, i] - p) ^ 2L) / sum((replicates[, i] - mean(replicates[, i])) ^ 2L)
      
      if(!silent) cat('CV leaving out replicate ', i, ', lambda = ', lambda_list[j],
                      '; R2 = ', round(cvR2[i, j], 3L), '\n', sep = '')
    }
  }
  
  cvR2
}

# Perform cross-validation to evaluate a series of nested models constructed by
# including reference-free effects in the order of their phenotypic contribution.
minimal_models <- function(genotypes, y, site_combn, model, max_size, k, lambda_list) {
  
  # genotypes, y, site_combn :
  #   As defined in `infer_model`.
  # model : list
  #   Model inferred from genotypes, y, and site_combn; used to calculate and
  #   rank the phenotypic contribution of each reference-free effect.
  # max_size : int
  #   Maximum number of effects to evaluate.
  # k, lambda_list :
  #   As defined in `cross_validation`.
  #
  # Returns the out-of-sample R2 of the minimal models.
  
  # Cross-validation function customized for this function. It takes the
  # phenotype operator as input instead of a list of site-combinations, which
  # allows the phenotype operator to be calcualted just once.
  cv_internal <- function(X, y, k, lambda_list, epsilon = 1e-3) {
    
    # X : Phenotype operator matrix
    # k : Number of genotype partitions for cross-validation.
    
    cvR2 <- matrix(NA_real_, k, length(lambda_list))
    
    # Genotype partitions for cross-validation.
    membership <- sample(1L:length(y) %% k + 1L) 
    
    for(i in 1L:k) { # i indexes genotype partitions
      
      # No phenotypic variance to explain in the test set.
      if(var(y[membership == i]) < 1e-10) next
      
      X_active <<- X[membership != i, ]
      y_active <<- y[membership != i]
      
      for(j in seq_along(lambda_list)) { # j indexes lasso penalties
        
        # Initial values for optimization.
        init_param <- rep(0, ncol(X_active))
        
        if(nonspec_type == 2L)
          init_param <- c(init_param, min(y_active), max(y_active) - min(y_active))
        
        res <- lbfgs(objective, gradient, init_param,
                     orthantwise_c = lambda_list[j], orthantwise_start = 1L,
                     invisible = 1L, epsilon = epsilon)$par
        
        e <- res[1L:ncol(X_active)] # Inferred reference-free effects
        
        if(nonspec_type == 1L)
          nonspec_param <- NULL
        else
          nonspec_param <- res[(ncol(X_active) + 1L):length(res)]
        
        # Predicted phenotype for the test set.
        p <- apply_nonspecific_epistasis(as.vector(X[membership == i, ] %*% e),
                                         nonspec_param)
        
        # Out-of-sample R2.
        cvR2[i, j] <-
          1 - sum((y[membership == i] - p) ^ 2L) / sum((y[membership == i] - mean(y[membership == i])) ^ 2L)
      }
    }
    
    cvR2
  }
  
  # Ordering reference-free effects by phenotypic contribution.
  
  # Fraction of phenotypic variance due to each effect.
  v <- unlist(partition_variance(site_combn, model$e))
  
  # Order of effects to include; intercept is always included.
  ord <- c(1L, order(v, decreasing = TRUE) + 1L)
  
  # Full phenotype operator.
  X <- construct_phenotype_operator(genotypes, site_combn)
  
  # Lasso is activated only when the number of parameters is sufficiently large.
  perform_lasso <- FALSE
  
  # Intercept-only model explains no variance.
  cat('Intercept: df = 1; R2 = '); cat(rep(0, k), sep = ', '); cat('\n')
  minR2 <- matrix(rep(0, k), nrow = 1L)
  
  # Degree of freedom of the current model.
  df <- 1L
  
  for(i in 2L:min(length(model$e), max_size)) {
    
    cat('Effect number ', i - 1, ':', sep = '')
    
    # Calculating the degree of freedom of the model.
    rank_i <- rankMatrix(X[, ord[1L:i]], method = 'qr')
    
    if(!(rank_i > df)) { # The effect is redundant with already included ones
      
      cat(' skipped because it is redundant\n')
      next
      
    } else {
      
      df <- df + 1L
      cat(' df = ', df, ';', sep = '')
    }
    
    if(perform_lasso) {
      
      cvR2 <- cv_internal(X[, ord[1L:i]], y, k, lambda_list)
      cvR2 <- cvR2[, which.max(colSums(cvR2))] # Best out-of-sample R2
      
    } else {
      
      if(i %% 20L == 0L) { # Perform lasso regression with each 20 additional parameters.
        
        cvR2 <- cv_internal(X[, ord[1L:i]], y, k, c(0, lambda_list))
        cvR2_mean <- apply(cvR2, 2L, mean)
        
        # If better model fit is found for a nonzero lambda, begin lasso.
        if(max(cvR2_mean) - cvR2_mean[1L] > 0.002) {
          
          cat('  Performing lasso regression from now on\n')
          perform_lasso <- TRUE
        }
        
        cvR2 <- cvR2[, which.max(cvR2_mean)] # Best out-of-sample R2
        
      } else {
        
        cvR2 <- as.vector(cv_internal(X[, ord[1L:i]], y, k, 0))
      }
    }
    
    cat(' R2 = '); cat(round(cvR2, 3L), sep = ', '); cat('\n')
    minR2 <- rbind(minR2, cvR2)
  }
  
  unname(minR2)
}

# Evaluate minimal models by performing cross-validation across measurement
# replicates.
minimal_models_replicates <- function(genotypes, replicates, site_combn, model,
                                      max_size, lambda_list) {
  
  # replicates : numeric matrix
  #   Replicate measurements; rows correspond to genotypes, columns to replicates.
  # All other arguments are defined as in `evaluate_minimal_models`.
  
  # Cross-validation function customized for this function. It takes the
  # phenotype operator as input instead of a list of site-combinations, which
  # allows the phenotype operator to be calcualted just once.
  cv_internal <- function(X, replicates, lambda_list, epsilon = 1e-3) {
    
    X_active <<- X
    cvR2 <- matrix(NA_real_, ncol(replicates), length(lambda_list))
    
    for(i in 1L:ncol(replicates)) {
      
      # Mean phenotype excluding the test replicate set.
      y_active <<- apply(replicates[, -i, drop = FALSE], 1L, mean, na.rm = TRUE)
      
      for(j in seq_along(lambda_list)) {
        
        # Initial parameter values.
        
        init_param <- rep(0, ncol(X))
        
        if(nonspec_type == 2L)
          init_param <- c(init_param, min(y_active), max(y_active) - min(y_active))
        
        res <- lbfgs(objective, gradient, init_param, orthantwise_c = lambda_list[j],
                     orthantwise_start = 1L, invisible = 1L, epsilon = epsilon)$par
        
        e <- res[1L:ncol(X)]
        
        if(nonspec_type == 1L)
          nonspec_param <- NULL
        else
          nonspec_param <- res[(ncol(X) + 1L):length(res)]
        
        # Predicted phenotype.
        p <- apply_nonspecific_epistasis(as.vector(X %*% e), nonspec_param)
        
        # Out-of-sample R2.
        cvR2[i, j] <-
          1 - sum((replicates[, i] - p) ^ 2L) / sum((replicates[, i] - mean(replicates[, i])) ^ 2L)
      }
    }
    
    cvR2
  }
 
  # Excluding genotypes missing any replicate measurement.
  exclude <- rowSums(is.finite(replicates)) < ncol(replicates)
  genotypes <- genotypes[!exclude, ]
  replicates <- replicates[!exclude, ]

  # Ordering reference-free effects by phenotypic contribution.
  
  # Fraction of phenotypic variance due to each effect.
  v <- unlist(partition_variance(site_combn, model$e))
  
  # Order of effects to include; intercept is always included.
  ord <- c(1L, order(v, decreasing = TRUE) + 1L)
  
  # Full phenotype operator.
  X <- construct_phenotype_operator(genotypes, site_combn)
  
  # Lasso is activated only when the number of parameters is sufficiently large.
  perform_lasso <- FALSE
  
  # Intercept-only model explains no variance.
  cat('Intercept: df = 1; R2 = '); cat(rep(0, k), sep = ', '); cat('\n')
  minR2 <- matrix(rep(0, k), nrow = 1L)
  
  # Degree of freedom of the current model.
  df <- 1L
  
  for(i in 2L:min(length(model$e), max_size)) {
    
    cat('Effect number ', i - 1, ':', sep = '')
    
    # Calculating the degree of freedom of the model.
    rank_i <- rankMatrix(X[, ord[1L:i]], method = 'qr')
    
    if(!(rank_i > df)) { # The effect is redundant with already included ones
      
      cat(' skipped because it is redundant\n')
      next
      
    } else {
      
      df <- df + 1L
      cat(' df = ', df, ';', sep = '')
    }
    
    if(perform_lasso) {
      
      cvR2 <- cv_internal(X[, ord[1L:i]], replicates, lambda_list)
      cvR2 <- cvR2[, which.max(colSums(cvR2))] # Best out-of-sample R2
      
    } else {
      
      if(i %% 20L == 0L) { # Perform lasso regression with each 20 additional parameters.
        
        cvR2 <- cv_internal(X[, ord[1L:i]], replicates, c(0, lambda_list))
        cvR2_mean <- apply(cvR2, 2L, mean)
        
        # If better model fit is found for a nonzero lambda, begin lasso.
        if(max(cvR2_mean) - cvR2_mean[1L] > 0.002) {
          
          cat('  Performing lasso regression from now on\n')
          perform_lasso <- TRUE
        }
        
        cvR2 <- cvR2[, which.max(cvR2_mean)] # Best out-of-sample R2
        
      } else {
        
        cvR2 <- as.vector(cv_internal(X[, ord[1L:i]], replicates, 0))
      }
    }
    
    cat(' R2 = '); cat(round(cvR2, 3L), sep = ', '); cat('\n')
    minR2 <- rbind(minR2, cvR2)
  }
  
  unname(minR2)
}

