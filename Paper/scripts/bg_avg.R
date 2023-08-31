
# This script implements functions for background-averaged analysis.

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

# Considering state 1 as the reference state, identify the mutation order of genotypes in their canonical arrangement.
identify_genotype_order <- function(n, s) {
  
  O1 <- diag(s, nrow = s, ncol = s)
  O1[1L, 1L] <- 1
  
  O <- O1
  
  if(n > 1L) {
    
    for(i in 2L:n) {
      O <- kronecker(O1, O)
    }
  }
  
  round(diag(log(O, s)))
}

# Generate the background-averaged epistasis operator.
construct_epistasis_operator <- function(n, s) {
  
  VH1 <- diag(1, nrow = s, ncol = s)
  VH1[, 1L] <- -1
  VH1[1L, ] <- 1 / s
  
  VH <- VH1
  
  if(n > 1L) {
    
    for(i in 2L:n) VH <- kronecker(VH1, VH)
  }
  
  VH
}

# Construct the phenotype operator matrix.
# This is the matrix that, when multiplied to a vector of background-averaged terms, yields the phenotypes.
construct_phenotype_operator <- function(n, s) {
  
  A1 <- diag(1, nrow = s, ncol = s) - 1 / s * matrix(1, nrow = s, ncol = s)
  A1[, 1L] <- 1
  
  A <- A1
  
  if(n > 1L) {
    
    for(i in 2L:n) A <- kronecker(A1, A)
  }
  
  A
}

# Simulate a GP map.
simulate_GP_map <- function(n, s, var_spec, spar_spec, noise_var, distr = 'normal') {
  
  # var_spec : numeric
  #   Variance spectrum, or the absolute genetic variance due to each epistatic order.
  # spar_spec : numeric
  #   Sparsity spectrum, or the fraction of terms that equal zero at each epistatic order.
  # noise_var : numeric
  #   Phenotypic variance due to noise.
  # distr : character ('normal', 'uniform')
  #   Distribution from which non-zero effects are drawn.
  
  ord <- identify_genotype_order(n, s)
  genotypes <- enumerate_genotypes(n, s)
  P <- construct_phenotype_operator(n, s)
  b <- rep(0, s ^ n)
  
  # Simulating terms without regards to variance spectrum.
  b_list <- lapply(1L:n, function(model_ord) {
    
    if(var_spec[model_ord] == 0) return(rep(0, sum(ord == model_ord)))
    
    if(distr == 'normal')
      b_model_ord <- rnorm(sum(ord == model_ord), 0, 1)
    else if(distr == 'uniform')
      b_model_ord <- runif(sum(ord == model_ord), -1, 1)
    
    b_model_ord[sample(sum(ord == model_ord), sum(ord == model_ord) * spar_spec[model_ord])] <- 0
    b_model_ord
  })
  
  # Adjusting variance spectrum.
  
  for(model_ord in 1L:n) {
    
    if(var_spec[model_ord] == 0) next
    
    y_model_ord <- as.vector(P[, ord == model_ord, drop = FALSE] %*% b_list[[model_ord]])
    var_model_ord <- mean((y_model_ord - mean(y_model_ord)) ^ 2L)
    b[ord == model_ord] <- b_list[[model_ord]] * sqrt(var_spec[model_ord] / var_model_ord)
  }
  
  # Adding measurement noise.
  y <- as.vector(P %*% b) + rnorm(s ^ n, 0, sqrt(noise_var))
  
  list(genotypes = genotypes, ord = ord, y = y, b = b)
}

# Identify the best model for a subsample of specified size via lasso regression and
# use it to predict the phenotype of every genotype.
predict_from_subsample <- function(n, s, genotypes, y, model_ord, frac_list, niter, CV_k) {
  
  prediction <- matrix(NA_real_, length(frac_list), niter, dimnames = list(frac_list))
  
  ord <- identify_genotype_order(n, s)
  A <- construct_phenotype_operator(n, s) # Indicator matrix (y = Ab)
  
  for(i in seq_along(frac_list)) {
    
    cat('Prediction with sample size', frac_list[i], '\n')
    
    for(j in 1L:niter) {
      
      train_set <- sample(s ^ n, s ^ n * frac_list[i]) # Sampling training/test set.
      
      if(var(y[-train_set]) < 1e-10) next # No phenotypic variation in the test set to explain
      
      # Cross-validation to determine the best-fit model for the training set.
      
      A_train <- A[train_set, ord <= model_ord & ord != 0L]
      cv <- cv.glmnet(A_train, y[train_set], family = 'gaussian', nfolds = CV_k, alpha = 1)
      res <- glmnet(A_train, y[train_set], family = 'gaussian', alpha = 1, lambda = cv$lambda.min)
      b_inferred <- rep(0, s ^ n)
      b_inferred[ord <= model_ord] <- coef(res)
      
      # Evaluation on the test set.
      
      y_inferred <- as.vector(A %*% b_inferred)
      
      prediction[i, j] <-
        1 - sum((y[-train_set] - y_inferred[-train_set]) ^ 2L) / sum((y[-train_set] - mean(y[-train_set])) ^ 2L)
    }
  }
  
  prediction
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

