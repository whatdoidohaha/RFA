
# Computing background-averaged effects for a complete, simulated dataset.

source('../../../scripts/bg_avg.R')

genotypes <- read.csv('../../../data/Simulated/2_8/genotypes.txt', header = FALSE)
y <- as.numeric(readLines('../../../data/Simulated/2_8/phenotypes.txt'))
n <- ncol(genotypes)
s <- max(genotypes) # Assumes the same number of states for every site.

ord <- identify_genotype_order(n, s)
H <- construct_epistasis_operator(n, s)

b <- as.vector(H %*% y)
writeLines(as.character(b), 'bg_avg_analysis.txt')

