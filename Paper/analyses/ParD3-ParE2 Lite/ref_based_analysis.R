
# Reference-based analysis.


# Setting and data import -------------------------------------------------

dataset <- 'ParD3-ParE2 Lite'
WT_index <- 954L # Row number of the wild-type genotype.
nonspec_type <- 2L # 1: No nonspecific epistasis; 2: Logistic nonspecific epistasis.
max_ord <- 2L # Maximum model order to be tested.
niter <- 500L # Number of random references to test.

source('../../scripts/ref_based.R')
import_data(paste0('../../data/', dataset, '/genotypes.txt'), paste0('../../data/', dataset, '/phenotypes.txt'))


# Run ---------------------------------------------------------------------

# If niter is greater than the number of possible genotypes, every genotype is used as reference once.
if(niter > nrow(genotypes)) comprehensive <- TRUE else comprehensive <- FALSE

if(comprehensive) { # Using each genotype as reference.
  
  ref <- rep(WT_index, nrow(genotypes) + 1L)
  res <- matrix(NA_real_, nrow(genotypes) + 1L, max_ord)
  
  cat('  Using WT as reference:\n')
  res[1L, ] <- sapply(1L:max_ord, function(model_ord) infer_model(WT_index, model_ord, nonspec_type)$R2_higher_ord)
  
  for(i in 1L:nrow(genotypes)) {
    
    cat('  Reference', i, ':\n')
    ref[i + 1L] <- i
    res[i + 1L, ] <- sapply(1L:max_ord, function(model_ord) infer_model(i, model_ord, nonspec_type)$R2_higher_ord)
  }
  
} else { # Randomly sampling references.
  
  ref <- rep(WT_index, niter + 1L)
  res <- matrix(NA_real_, niter + 1L, max_ord)
  
  cat('  Using WT as reference:\n')
  res[1L, ] <- sapply(1L:max_ord, function(model_ord) infer_model(WT_index, model_ord, nonspec_type)$R2_higher_ord)
  
  for(i in 1L:niter) {
    
    cat('  Random reference iteration', i, ':\n')
    ref[i + 1L] <- sample(1L:nrow(genotypes), 1L)
    res[i + 1L, ] <- sapply(1L:max_ord, function(model_ord) infer_model(ref[i + 1L], model_ord, nonspec_type)$R2_higher_ord)
  }
}

res <- cbind(ref, round(res, 3L))
colnames(res) <- c('ref', 1L:max_ord)

write.csv(res, 'ref_based_analysis.txt', row.names = FALSE)
