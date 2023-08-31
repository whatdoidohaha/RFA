
# This script implements functions for reference-free analysis of complete or near-complete datasets.

library(Matrix)
library(abind)


# Import model settings from a control file.
import_control_file <- function(dir) {
  
  # Control file must contain two lines:
  #   1. Path to a csv file listing the genotypes.
  #   2. Path to a text file listing the phenotypes (line by line for each genotype).
  
  input <- readLines(dir)
  lapply(strsplit(input, split = '#'), function(x) trimws(x[1L], 'right'))
}

# Set model as specified in a control file.
set_model <- function(control) {
  
  # The following global variables are created:
  #
  # `genotypes` : Genotype matrix.
  # `y` : Phenotype vector.
  # `n` : Number of sites.
  # `q_list` : Number of states in each site.
  # `site_combn` : List of site-combinations for which effects are inferred.
  # `global_type` : Type of globalific nonlinearity (1 = linear, 2 = logistic).
  #
  # States are renumbered to be 1, 2, 3, ...
  
  genotypes <<- unname(as.matrix(read.csv(control[[1L]], header = FALSE))) # Genotypes
  if(!all(is.finite(genotypes))) stop('Genotype matrix contains missing entries.\n')
  
  y <<- as.numeric(readLines(control[[2L]])) # Phenotypes
  if(!all(is.finite(y))) stop('Phenotype vector contains missing entries\n')
  if(nrow(genotypes) != length(y)) stop('Genotype matrix and phenotype vector of different lengths.\n')
  
  n <<- ncol(genotypes) # Number of sites
  q_list <<- apply(genotypes, 2L, function(x) length(unique(x))) # Number of observed states in each site
  
  # Renumbering states if necessary.
  
  ori_q_list <- lapply(1L:ncol(genotypes), function(j) sort(unique(genotypes[, j]))) # States in each site
  
  reducible <- vapply(ori_q_list, function(x) length(x) != max(x), TRUE) # Sites to renumber states.
  
  if(any(reducible)) {
    
    cat('States should be numbered 1, 2, 3, ...\n')
    cat('States renumbered for sites', which(reducible), '\n')
    genotypes <<- vapply(1L:ncol(genotypes), function(j) match(genotypes[, j], ori_q_list[[j]]), genotypes[, 1L])
  }
}

# Generate all permutations of elements in a list of vectors, right-most changing fastest.
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

# Calculate ensemble phenotypes, reducing the state space by omitting the last state in each site.
transform_to_ensemble_representation <- function() {
  
  # Hash table; element [[i]]][[s]] is a logical vector indicating whether each genotype has state s at site i.
  hash <- lapply(1L:n, function(j) lapply(1L:q_list[j], function(s) genotypes[, j] == s))
  
  # Ensemble genotypes, ordered as the right-most site fastest-changing.
  genotypes <<- permutations(lapply(q_list, function(qj) 0L:(qj - 1L)))
  
  # Ensemble phenotypes.
  y <<- apply(genotypes, 1L, function(g) {
    
    sites <- which(g != 0L)
    if(length(sites) == 0L) return(mean(y)) # Intercept
    states <- g[sites]
    
    mean(y[Reduce(`&`, mapply(function(i, s) hash[[i]][[s]], sites, states, SIMPLIFY = FALSE))])
  })
}

# Calculate the reference-based epistasis operator by recursion.
# This is the matrix that, when multiplied to a vector of phenotypes, yields the reference-based effects.
calculate_epistasis_operator <- function(G0 = sparseMatrix(1L, 1L, x = 1, dims = c(1L, 1L)), k = 0L) {
  
  if(k == n) return(G0) # Recursion complete.
  
  Z <- sparseMatrix(c(), c(), dims = dim(G0))
  
  G1 <- do.call(cbind, lapply(1L:q_list[n - k], function(j) {
    
    if(j == 1L) {
      
      do.call(rbind, c(list(G0), rep(list(-G0), q_list[n - k] - 1L)))
      
    } else {
      
      do.call(rbind, c(rep(list(Z), j - 1L), list(G0), rep(list(Z), q_list[n - k] - j)))
    }
  }))
  
  calculate_epistasis_operator(G1, k + 1L)
}

# Calculate reference-free effects.
calculate_effects <- function() {
  
  G <- calculate_epistasis_operator()
  
  if(all(is.finite(y))) { # No missing ensemble phenotypes
    
    as.vector(G %*% y)
    
  } else {
    
    missing_effects <- !is.finite(y)
    effects <- as.vector(G[, !missing_effects] %*% y[!missing_effects])
    effects[missing_effects] <- NA_real_
    effects
  }
}

# Parse computed effects and compute effects including the last state in each site.
# Sorted first by order, then by site-combination (same order as in `combn`), then by state-combination (in an array).
parse_effects <- function(effects) {

  ord <- Reduce(`+`, lapply(1L:n, function(j) genotypes[, j] != 0L))
  
  lapply(0L:n, function(k) {
    
    if(k == 0L) return(effects[1L])
    
    indices <- which(ord == k)
    
    if(k == 1L)
      sites <- matrix(apply(genotypes[indices, , drop = FALSE], 1L, function(g) which(g != 0L)), nrow = 1L)
    else
      sites <- apply(genotypes[indices, , drop = FALSE], 1L, function(g) which(g != 0L))
    
    site_ord <- do.call(order, lapply(1L:k, function(j) sites[j, ]))
    
    # Effects ordered by site-combination (same order as `combn`) and by state-combination (last site moving fastest).
    effects_k <- effects[indices[site_ord]]
    site_combn <- combn(n, k)
    n_terms <- Reduce(`*`, lapply(1L:k, function(j) q_list[site_combn[j, ]] - 1L))
    pos <- as.integer(cumsum(c(0L, n_terms)))
    
    # Parsed effects.
    parsed <- lapply(1L:ncol(site_combn), function(j) effects_k[(pos[j] + 1L):pos[j + 1L]])
    
    # Adding omitted effects.
    if(k == 1L) {
      
      for(j in seq_along(parsed)) parsed[[j]] <- c(parsed[[j]], -sum(parsed[[j]]))
      
    } else if(k == 2L) {
      
      for(j in seq_along(parsed)) {
        qj <- q_list[site_combn[, j]]
        E <- t(matrix(parsed[[j]], nrow = qj[2L] - 1L))
        parsed[[j]] <- cbind(rbind(E, -colSums(E)), c(-rowSums(E), sum(E)))
      }
      
    } else {
      
      for(j in seq_along(parsed)) {
        
        qj <- q_list[site_combn[, j]]
        E <- aperm(array(parsed[[j]], rev(qj - 1L)), k:1L)
        
        # Incrementally expanding the effect array.
        for(i in 1L:k) {
          
          margin <- setdiff(1L:k, i)
          marginal_E <- array(-apply(E, margin, sum), c(dim(E)[margin], 1L))
          marginal_E <- aperm(marginal_E, order(c(margin, i)))
          E <- abind(E, marginal_E, along = i)
        }
        
        parsed[[j]] <- E
      }
    }
    
    parsed
  })
}


# Run example.

#control <- import_control_file('scripts/control.txt')
#set_model(control)
#transform_to_ensemble_representation()
#e <- calculate_effects()
#e_parsed <- parse_effects(e) # Costly operation



