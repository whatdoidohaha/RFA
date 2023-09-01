
# Computing reference-free effects for a complete, simulated dataset.

source('../../../scripts/ref_free_complete_data.R')
 
set_model(list('../../../data/Simulated/2_8/genotypes.txt', '../../../data/Simulated/2_8/phenotypes.txt'))
transform_to_ensemble_representation()
e <- calculate_effects()
e_ord <- rowSums(genotypes != 0L) # Order of each effect

writeLines(as.character(e), 'ref_free_analysis.txt')
