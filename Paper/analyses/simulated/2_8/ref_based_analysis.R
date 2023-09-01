
# Computing reference-based effects for a complete, simulated dataset.

source('../../../scripts/ref_based.R')

import_data('../../../data/Simulated/2_8/genotypes.txt', '../../../data/Simulated/2_8/phenotypes.txt')
model <- infer_model(1L, n, 1L) # Using genotype 1 as reference; no nonspecific epistasis
ord <- model$mut_ord

writeLines(as.character(model$e), 'ref_based_analysis.txt')
