library(Matrix)
library(lbfgs)


# Import data.
import_data <- function(genotype_dir, phenotype_dir) {
  
  # genotype_dir : Path to a csv file listing the genotypes.
  # phenotype_dir : Path to a text file listing the phenotypes (line by line for each genotype).
  #
  # The following global variables are created:
  #
  # genotypes : Genotype matrix.
  # y : Phenotype vector.
  # n : Number of sites.
  # q_list : Number of states in each site.
  #
  # States are numbered from 1, 2, 3, ...

  genotypes <<- unname(as.matrix(read.csv(genotype_dir, header = FALSE))) # Genotypes
  if(!all(is.finite(genotypes))) stop('Genotype matrix contains missing entries.\n')
  
  y <<- as.numeric(readLines(phenotype_dir)) # Phenotypes
  if(!all(is.finite(y))) stop('Phenotype vector contains missing entries\n')
  if(nrow(genotypes) != length(y)) stop('Genotype matrix and phenotype vector of different lengths.\n')
  
  n <<- ncol(genotypes) # Number of sites
  q_list <<- apply(genotypes, 2L, function(x) length(unique(x))) # Number of observed states in each site
  
  # Renumbering states if necessary.
  
  ori_q_list <- lapply(1L:ncol(genotypes), function(j) sort(unique(genotypes[, j]))) # States in each site
  
  reducible <- vapply(ori_q_list, function(x) length(x) != max(x), TRUE) # Sites to renumber states.
  
  if(any(reducible)) {
    cat('State space reduced for sites', which(reducible), '\n')
    genotypes <<- vapply(1L:ncol(genotypes), function(j) match(genotypes[, j], ori_q_list[[j]]), genotypes[, 1L])
  }
}

# Infer a specified model.
infer_model <- function(ref, model_ord, nonspec_type, calc_R2_higher_order = FALSE) {
  
  # ref : integer
  #   Row number of the reference genotype.
  # model_ord : integer
  #   Model order.
  # nonspec_type : integer
  #   Type of nonspecific epistasis (1 = none, 2 = logistic).
  # calc_R2_higher_order : logical
  #   If TRUE, calculates the r2 between predicted and observed phenotypes of higher-order mutants.
  #
  # Returns a list containing the following elements:
  #
  # ref_g : Reference genotype.
  # ref_y : Phenotype of the reference genotype.
  # model_ord : Model order.
  # nonspec_type : Type of nonspecific epistasis.
  # nonspec_param : Inferred parameters of the nonspecific epistasis model.
  # mut_ord : Order of each mutant.
  # calculable : Calculable effects and genotypes involved in their calculation.
  # e : Value of the calculable effects.
  
  # Reference genotype.
  
  if(ref > nrow(genotypes)) stop('Invalid reference genotype.\n')
  ref_g <- genotypes[ref, ]
  ref_y <- y[ref]
  
  # Mutant order.
  mut_ord <- rowSums(vapply(1L:n, function(i) genotypes[, i] != ref_g[i], logical(nrow(genotypes))))
  
  # Identifying the calculable effects and the genotypes involved in their calculation.
  calculable <- identify_calculable_effects(ref_g, model_ord)
  
  # Deriving the phenotype operator.
  # This is the matrix that, when multiplied to the vector of effects, yields the phenotypes.
  G <- construct_phenotype_operator(genotypes[calculable$indices, ], calculable$effects)
  
  if(nonspec_type == 1L) { # Not modeling nonspecific epistasis.
    
    nonspec_param <- NULL
    
    e <- as.vector(solve(G, y[calculable$indices]))
    
    # Fit betwen the predicted and observed phenotype (all genotypes).
    
    G_all <- construct_phenotype_operator(genotypes, calculable$effects)
    p <- as.vector(G_all %*% e)
    R2 <- 1 - sum((y - p) ^ 2L) / sum((y - mean(y)) ^ 2L)
    
    # Predicting the phenotypes of higher-order mutants.
    
    R2_higher_ord <- NA_real_
    
    if(calc_R2_higher_order) {
      
      y_higher_ord <- y[mut_ord > model_ord]
      
      # Operator for calculating the phenotypes of higher-order mutants from the calculable effects.
      G_higher_ord <- construct_phenotype_operator(genotypes[mut_ord > model_ord, , drop = FALSE], calculable$effects)
      
      # Total number of effects involved in the calculation of each high-order phenotype.
      n_terms_required <- sapply(mut_ord[mut_ord > model_ord], function(ord) sum(sapply(0L:model_ord, function(i) choose(ord, i))))
      
      # Number of calculable effects for each high-order phenotype.
      n_terms_calculable <- rowSums(G_higher_ord)
      
      # Examining only the high-order genotypes for which all relevant effects are calculable.
      
      if(sum(n_terms_required == n_terms_calculable) > 0.1 * length(n_terms_required)) {
        
        p_valid <- as.vector(G_higher_ord[n_terms_required == n_terms_calculable, ] %*% e)
        y_valid <- y_higher_ord[n_terms_required == n_terms_calculable]
        
        if(var(y_valid) > 1e-10) # Only when there is phenotypic variation among higher-order mutants
          R2_higher_ord <- 1 - sum((y_valid - p_valid) ^ 2L) / sum((y_valid - mean(y_valid)) ^ 2L)
        
      } else {
        
        cat('Too few higher-order genotypes available for checking extrapolability.\n')
      }
    }
      
    return(list(ref_g = ref_g, ref_y = ref_y, model_ord = model_ord, nonspec_type = nonspec_type, nonspec_param = nonspec_param,
                mut_ord = mut_ord, calculable = calculable, e = e, R2 = R2, R2_higher_ord = R2_higher_ord))
    
  } else { # Modeling nonspecific epistasis.
  
    # Identifying higher-order genotypes for which all relevant effects can be calculated.
    
    y_higher_ord <- y[mut_ord > model_ord]
    
    # Operator for calculating the phenotypes of higher-order mutants from the calculable effects.
    G_higher_ord <- construct_phenotype_operator(genotypes[mut_ord > model_ord, ], calculable$effects)
    
    # Total number of effects involved in the calculation of each high-order phenotype.
    n_terms_required <- sapply(mut_ord[mut_ord > model_ord], function(ord) sum(sapply(0L:model_ord, function(i) choose(ord, i))))
    
    # Number of calculable effects for each high-order phenotype.
    n_terms_calculable <- rowSums(G_higher_ord)
    
    if(sum(n_terms_required == n_terms_calculable) > 0.1 * length(n_terms_required)) {
      
      G_higher_ord <- G_higher_ord[n_terms_required == n_terms_calculable, ]
      y_higher_ord <- y_higher_ord[n_terms_required == n_terms_calculable]
      
      if(var(y_higher_ord) > 1e-10) # Only when there is phenotypic variation among higher-order mutants
        res <- infer_nonspecific_epistasis(calculable$effects, G, y[calculable$indices], G_higher_ord, y_higher_ord)
      else
        res <- list('nonspec_param' = NA_real_, 'e' = NA_real_, 'R2_higher_ord' = NA_real_)
        
      return(list(ref_g = ref_g, ref_y = ref_y, model_ord = model_ord, nonspec_type = nonspec_type, nonspec_param = res$nonspec_param,
                  mut_ord = mut_ord, calculable = calculable, e = res$e, R2 = NA_real_, R2_higher_ord = res$R2_higher_ord))
    
    } else {
      
      cat('Too few higher-order genotypes can be used for inferring nonspecific epistasis.\n')
      
      return(list(ref_g = ref_g, ref_y = ref_y, model_ord = model_ord, nonspec_type = nonspec_type, nonspec_param = NA_real_,
                  mut_ord = mut_ord, calculable = calculable, e = NA_real_, R2 = NA_real_, R2_higher_ord = NA_real_))
    }
  }
}


## Required specialized functions.

# Identify the effects that can be calculated and the genotypes required for their calculation.
identify_calculable_effects <- function(ref_g, model_ord, verbose = TRUE) {
  
  # ref_g : integer
  #   Index of the reference genotype.
  # model_ord : integer
  #   Model order.
  #
  # Returns a list with the following elements:
  #
  # effects : A matrix indicating the calculable effects (the identity of the effects, not their values).
  # indices : A vector of indices indicating the genotypes involved in the calculation of calculable effects.
  
  # Distance from the reference genotype.
  dist <- Reduce(`+`, lapply(1L:n, function(j) genotypes[, j] != ref_g[j]))
  
  # If all required genotypes are present:
  
  n_terms <- lapply(0L:model_ord, function(ord) apply(combn(n, ord), 2L, function(x) prod(q_list[x] - 1L)))
  
  if(sum(unlist(n_terms)) == sum(dist <= model_ord)) {
    
    if(verbose) cat('All effects calculable.\n')
    
    effects <- list(matrix(ref_g, nrow = 1L)) # List to be filled with calculable effects.
    indices <- list(which(dist == 0L)) # List to be filled with genotypes required for calculation.
    
    for(ord in 1L:model_ord) {
      
      indices_ord <- which(dist == ord)
      effects_ord <- vapply(1L:n, function(j) {g <- genotypes[indices_ord, j]; g[g == ref_g[j]] <- 0L; g}, indices_ord)
      if(!is.matrix(effects_ord)) effects_ord <- matrix(effects_ord, nrow = 1L)
      
      # Ordering by site-combination.
      sites <- apply(effects_ord, 1L, function(e) which(e != 0L))
      if(ord == 1L)
        site_ord <- order(sites)
      else
        site_ord <- do.call(order, lapply(1L:ord, function(j) sites[j, ]))
      
      indices[[ord + 1L]] <- indices_ord[site_ord]
      effects[[ord + 1L]] <- effects_ord[site_ord, , drop = FALSE]
    }
    
    return(list(effects = do.call(rbind, effects), indices = do.call(c, indices)))
  }
  
  # Identifying the calculable effects given the pattern of missing data.
  
  if(verbose) cat('Some effects are not calculable due to missing genotypes.\n')
  
  effects <- list(matrix(ref_g, nrow = 1L)) # List to be filled with calculable effects.
  indices <- list(which(dist == 0L)) # List to be filled with genotypes required for calculation.
  
  for(ord in 1L:model_ord) {
    
    if(ord == 1L) { # First-order effects.
      
      indices_ord <- which(dist == ord)
      if(length(indices_ord) == 0L) stop('All point mutants missing; no effect is calculable.\n')
      
      # Recording calculable effects.
      effects_ord <- vapply(1L:n, function(j) {g <- genotypes[indices_ord, j]; g[g == ref_g[j]] <- 0L; g}, indices_ord)
      if(!is.matrix(effects_ord)) effects_ord <- matrix(effects_ord, nrow = 1L)
      
      # Ordering by site.
      sites <- apply(effects_ord, 1L, function(e) which(e != 0L))
      effects[[ord + 1L]] <- effects_ord[order(sites), , drop = FALSE]
      indices[[ord + 1L]] <- indices_ord[order(sites)]
      next
    }
    
    # Second and higher-order effects.
    
    indices_ord <- which(dist == ord)
    if(length(indices_ord) == 0L) {
      if(verbose) cat('Maximum possible calculable order:', ord - 1L, '\n')
      break
    }
    
    # All potentially calculable effects (those for which the highest-order mutant is observed).
    effects_ord <- vapply(1L:n, function(j) {g <- genotypes[indices_ord, j]; g[g == ref_g[j]] <- 0L; g}, indices_ord)
    if(!is.matrix(effects_ord)) effects_ord <- matrix(effects_ord, nrow = 1L)
    
    # Checking calculability.
    
    # Hash table; element [[i]]][[s]] is a logical vector indicating whether each lower-order effect has state s at site i.
    hash <- lapply(1L:n, function(j) lapply(1L:q_list[j], function(s) effects[[ord]][, j] == s))
    
    # Checking if each effect in calculable, i.e., all lower-order terms are calculable.
    is_calculable <- apply(effects_ord, 1L, function(e) {
      
      sites <- which(e != 0L)
      states <- e[sites]
      low_ord_sites <- combn(sites, ord - 1L, simplify = FALSE) # Required lower-order terms; sites
      low_ord_states <- combn(states, ord - 1L, simplify = FALSE) # Required lower-order terms; states
      is_calculable <- TRUE
      
      # Checking if each required lower-order term is calculable.
      for(k in 1L:length(low_ord_sites)) {
        
        is_calculable_k <- mapply(function(i, s) hash[[i]][[s]], low_ord_sites[[k]], low_ord_states[[k]], SIMPLIFY = FALSE)
        is_calculable <- any(Reduce(`&`, is_calculable_k))
        if(!is_calculable) break
      }
      
      is_calculable
    })
    
    if(sum(is_calculable) == 0L) {
      if(verbose) cat('Maximum possible calculable order:', ord - 1L, '\n')
      break
    }
    
    indices_ord <- indices_ord[is_calculable]
    effects_ord <- effects_ord[is_calculable, , drop = FALSE]
    
    # Ordering effects by site-combination.
    
    sites <- apply(effects_ord, 1L, function(e) which(e != 0L))
    site_ord <- do.call(order, lapply(1L:ord, function(j) sites[j, ]))
    indices[[ord + 1L]] <- indices_ord[site_ord]
    effects[[ord + 1L]] <- effects_ord[site_ord, , drop = FALSE]
  }
  
  list(effects = do.call(rbind, effects), indices = do.call(c, indices))
}

# Construct the phenotype operator.
# This is the matrix that, when multiplied to the vector of reference-based effects, yields the phenotypes.
construct_phenotype_operator <- function(genotypes, effects, as_pattern_matrix = FALSE) {
  
  # genotypes : numeric matrix
  #   Genotypes for which the phenotype should be calculated.
  # effects : numeric matrix
  #   Reference-based effects to be used for phenotypic calculation; the first row must be the intercept.
  # as_pattern_matrix : logical
  #   If TRUE, a sparse pattern matrix is returned; otherwise a sparse double matrix.
  
  # Hash table; element [[i]]][[s]] is a logical vector indicating whether each genotype has state s at site i.
  hash <- lapply(1L:n, function(j) lapply(1L:q_list[j], function(s) genotypes[, j] == s))
  
  # Calculate the row and column indices of the sparse genotype indicator matrix.
  ab <- lapply(1L:nrow(effects), function(k) {
    
    if(k == 1L) return(list(1L:nrow(genotypes), rep(1L, nrow(genotypes)))) # Intercept
    
    sites <- which(effects[k, ] != 0L) # An epistatic effect; sites
    states <- effects[k, sites] # An epistatic effect; states
    
    # Identify all genotypes with `states` in `sites`.
    is_present <- Reduce(`&`, mapply(function(i, s) hash[[i]][[s]], sites, states, SIMPLIFY = FALSE))
    list(which(is_present), rep(k, sum(is_present)))
  })
  
  a <- unlist(lapply(ab, `[[`, 1L)) # Row indices
  b <- unlist(lapply(ab, `[[`, 2L)) # Column indices
  
  if(as_pattern_matrix)
    sparseMatrix(a, b, dims = c(nrow(genotypes), nrow(effects)))
  else
    sparseMatrix(a, b, x = 1, dims = c(nrow(genotypes), nrow(effects)))
}

# Infer nonspecific epistasis.
infer_nonspecific_epistasis <- function(effects, G_lower_ord, y_lower_ord, G_higher_ord, y_higher_ord) {
  
  # For a reference-based model of a given order,
  # nonspecific epistasis is inferred by maximizing the extrapolability to higher-order mutants.
  # By the definition of reference-based analysis, the model must exactly reproduce
  # the phenotypes of low-order mutants.
  #
  # effects : numeric matrix
  #   Reference-based effects included in the model.
  # G_lower_ord : numeric matrix
  #   Phenotype operator matrix for the mutants whose reference-based effects are calculated.
  #   The reference-based model must exactly reproduce their phenotypes.
  # y_lower_ord : numeric vector
  #   Phenotype vector for the low-order mutants.
  # G_higher_ord : numeric matrix
  #   Phenotype operator matrix for the higher-order mutants.
  # y_higher_ord : numeric vector
  #   Phenotype vector for the higher-order mutants.
  
  # Infer a logistic nonspecific epistasis function.
  fit <- function(param, return_e = FALSE) {
    
    # Here, param[1L] = lower bound and param[2L] = upper bound.
    # This is different from other scripts, in which param[2L] = range (upper bound - lower bound).
    # This is for use in L-BFGS-B box constraint.
    
    # If the logistic function cannot cover the phenotypic range of lower-order genotypes.
    if(param[1L] > min(y_lower_ord) || param[2L] <= max(y_lower_ord)) return(0) 
    
    s_lower_ord <- log((y_lower_ord - param[1L]) / (param[2L] - y_lower_ord)) # Genetic score of the lower-ord genotypes
    e <- as.vector(solve(G_lower_ord, s_lower_ord)) # Reference-based effects on genetic score
    s_higher_ord <- as.vector(G_higher_ord %*% e) # Genetic score of the higher-order mutants.
    y_higher_ord_pred <- param[1L] + (param[2L] - param[1L]) / (1 + exp(-s_higher_ord))
    
    R2 <- 1 - sum((y_higher_ord - y_higher_ord_pred) ^ 2L) / sum((y_higher_ord - mean(y_higher_ord)) ^ 2L)
    
    if(return_e) return(e)
    
    if(!is.finite(R2)) 0 else -R2
  }
  
  res <- optim(c(min(y_lower_ord) - 0.1, max(y_lower_ord) + 0.1), fit, method = 'L-BFGS-B',
               lower = c(-Inf, max(y_lower_ord) + 0.01), upper = c(min(y_lower_ord) - 0.01, Inf))
  
  list(nonspec_param = c(res$par[1L], res$par[2L] - res$par[1L]), e = fit(res$par, return_e = TRUE),
       R2_higher_ord = -res$val)
}


## Functions for simulating GP maps.

# Enumerate all genotypes
sim_enumerate_genotypes <- function(n, s) {
  
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
sim_identify_genotype_order <- function(n, s) {
  
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

# Generate the reference-based phenotype operator.
# This is the matrix that, when multiplied to a vector of reference-based terms (in canonical order), yields the phenotypes.
sim_construct_phenotype_operator <- function(n, s) {
  
  Q1 <- diag(1, nrow = s, ncol = s)
  Q1[-1L, 1L] <- 1
  
  Q <- Q1
  
  if(n > 1L) {
    
    for(i in 2L:n) {
      Q <- kronecker(Q1, Q)
    }
  }
  
  Q
}

# Generate the reference-based epistasis operator.
# This is the matrix that, when multiplied to a vector of phenotypes (in canonical order), yields the reference-based terms.
sim_construct_epistasis_operator <- function(n, s) {
  
  R1 <- diag(1, nrow = s, ncol = s)
  R1[-1L, 1L] <- -1
  
  R <- R1
  
  if(n > 1L) {
    
    for(i in 2L:n) {
      R <- kronecker(R1, R)
    }
  }
  
  R
}

