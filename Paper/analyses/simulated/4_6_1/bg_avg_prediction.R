
# This script simulates prediction from subsampling under the background-averaged formalism.

library(glmnet)
source('../../../scripts/bg_avg.R')
dataset <- '4_6_1'

genotypes <- read.csv(paste0('../../../data/simulated/', dataset, '/genotypes.txt'), header = FALSE)
y <- as.numeric(readLines(paste0('../../../data/simulated/', dataset, '/phenotypes.txt')))
n <- ncol(genotypes)
s <- max(genotypes)

model_ord <- 2L
frac_list <- 1 / 2 ^ seq(7, 1, by = -1)
niter <- 250L
CV_k <- 5L

pred_bg_avg <- predict_from_subsample(n, s, genotypes, y, model_ord, frac_list, niter, CV_k)

dotplot_with_sd(log(frac_list, 2), pred_bg_avg, ylim = c(-1, 1))

write.csv(round(t(pred_bg_avg), 3L), 'bg_avg_prediction.txt', row.names = FALSE)
