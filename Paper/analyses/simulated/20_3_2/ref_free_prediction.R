
# This script simulates prediction from subsampling under the reference-free formalism.

library(glmnet)
source('../../../scripts/ref_free.R')
dataset <- '20_3_2'

set_model(list(paste0('../../../data/simulated/', dataset, '/genotypes.txt'),
               paste0('../../../data/simulated/', dataset, '/phenotypes.txt'), '1', 1))

# Modified version of the function predict_from_subsample in ref_free.R.
# Lasso regression is performed using the package glmnet.
predict_from_subsample <- function(genotypes, y, model_ord, frac_list, niter, CV_k) {
  
  prediction <- matrix(NA_real_, length(frac_list), niter, dimnames = list(frac_list))
  
  A <- construct_phenotype_operator(genotypes, generate_site_combn(model_ord))
  
  for(i in seq_along(frac_list)) {
    
    cat('Prediction with sample size', frac_list[i], '\n')
    
    for(j in 1L:niter) {
      
      train_set <- sample(nrow(genotypes), nrow(genotypes) * frac_list[i]) # Sampling training/test set.
      
      if(var(y[-train_set]) < 1e-10) next # No phenotypic variation to explain
      
      # Cross-validation to determine the best-fit model for the training set.
      
      A_train <- A[train_set, -1L] # Excluding the intercept
      cv <- cv.glmnet(A_train, y[train_set], family = 'gaussian', nfolds = CV_k, alpha = 1)
      res <- glmnet(A_train, y[train_set], family = 'gaussian', alpha = 1, lambda = cv$lambda.min)
      e_inferred <- coef(res)
      
      # Evaluation on the test set.
      
      y_inferred <- as.vector(A %*% e_inferred)
      
      prediction[i, j] <-
        1 - sum((y[-train_set] - y_inferred[-train_set]) ^ 2L) / sum((y[-train_set] - mean(y[-train_set])) ^ 2L)
    }
  }
  
  prediction
}


model_ord <- 2L
frac_list <- 1 / 2 ^ seq(7, 1, by = -1)
niter <- 200L
CV_k <- 5L

pred_ref_free <- predict_from_subsample(genotypes, y, model_ord, frac_list, niter, CV_k)

dotplot_with_sd(log(frac_list, 2), pred_ref_free, ylim = c(0, 1))

write.csv(round(t(pred_ref_free), 3L), 'ref_free_prediction.txt', row.names = FALSE)
