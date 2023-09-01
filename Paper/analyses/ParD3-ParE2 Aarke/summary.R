
# This script visualizes the analysis results.


# Main data import --------------------------------------------------------

source('../../scripts/ref_free.R')
dataset <- 'ParD3-ParE2 Aarke'

# Control file:
#   1. Directory to the genotype matrix.
#   2. Directory to the phenotype vector.
#   3. Model order (or directory to a site-combination list).
#   4. Type of nonspecific epistasis.
control <- list(paste0('../../data/', dataset, '/genotypes.txt'), paste0('../../data/', dataset, '/phenotypes.txt'), '1', 2)
set_model(control)


# 1. Distribution of phenotype ---------------------------------------------

h <- hist(y, breaks = 20)
h$counts <- h$counts / sum(h$counts)
plot(h, ylim = c(0, 1))


# 2. Cross-validation -----------------------------------------------------

cv_first_linear <- identify_lambda(dir = 'cv_R2_1st_linear.txt', show_plot = TRUE, also_return_R2 = TRUE, 
                                   xlim = c(-3, 1), ylim = c(0, 1), cex = 1.6, pch = 19L)
cv_second_linear <- identify_lambda(dir = 'cv_R2_2nd_linear.txt', show_plot = TRUE, also_return_R2 = TRUE, 
                                    xlim = c(-3, 1), ylim = c(0, 1), cex = 1.6, pch = 19L)
cv_third_linear <- identify_lambda(dir = 'cv_R2_3rd_linear.txt', show_plot = TRUE, also_return_R2 = TRUE, 
                                   xlim = c(-3, 1), ylim = c(0, 1), cex = 1.6, pch = 19L)

cv_first <- identify_lambda(dir = 'cv_R2_1st.txt', show_plot = TRUE, also_return_R2 = TRUE,
                            xlim = c(-3, 1), ylim = c(0, 1), cex = 1.6, pch = 19L)
cv_second <- identify_lambda(dir = 'cv_R2_2nd.txt', show_plot = TRUE, also_return_R2 = TRUE,
                             xlim = c(-3, 1), ylim = c(0, 1), cex = 1.6, pch = 19L)
cv_third <- identify_lambda(dir = 'cv_R2_3rd.txt', show_plot = TRUE, also_return_R2 = TRUE, p_cutoff = 0.5,
                            xlim = c(-3, 1), ylim = c(0, 1), cex = 1.6, pch = 19L)


# 3. Cross-validated model: Nonspecific epistasis -------------------------

M1 <- import_model('regularized_fit_1st.txt')
M2 <- import_model('regularized_fit_2nd.txt')
M3 <- import_model('regularized_fit_3rd.txt')

G1 <- construct_phenotype_operator(genotypes, generate_site_combn(1L))
G2 <- construct_phenotype_operator(genotypes, generate_site_combn(2L))
G3 <- construct_phenotype_operator(genotypes, generate_site_combn(3L))

# Genetic scores.
s1 <- as.vector(G1 %*% M1$e)
s2 <- as.vector(G2 %*% M2$e)
s3 <- as.vector(G3 %*% M3$e)

# Histogram of genetic score.
h1 <- hist(s1, breaks = 20, plot = FALSE)
h1$counts <- h1$counts / sum(h1$counts)
h2 <- hist(s2, breaks = 20, plot = FALSE)
h2$counts <- h2$counts / sum(h2$counts)
h3 <- hist(s3, breaks = 20, plot = FALSE)
h3$counts <- h3$counts / sum(h3$counts)

# Overlaid with nonspecific epistasis.
plot(h1)
plot(seq(min(h1$breaks), max(h1$breaks), length.out = 100L),
     apply_nonspecific_epistasis(seq(min(h1$breaks), max(h1$breaks), length.out = 100L), M1$nonspec_param),
     type = 'l', col = 'red')

plot(h2)
plot(seq(min(h2$breaks), max(h2$breaks), length.out = 100L),
     apply_nonspecific_epistasis(seq(min(h2$breaks), max(h2$breaks), length.out = 100L), M2$nonspec_param),
     type = 'l', col = 'red')

plot(h3)
plot(seq(min(h3$breaks), max(h3$breaks), length.out = 100L),
     apply_nonspecific_epistasis(seq(min(h3$breaks), max(h3$breaks), length.out = 100L), M3$nonspec_param),
     type = 'l', col = 'red')

# Fraction of genotypes in the lower 2.5% or upper 2.5% of the dynamic range.
frac_out <- (sum(y < M2$nonspec_param[1L] + 0.025 * M2$nonspec_param[2L]) + 
               sum(y > M2$nonspec_param[1L] + 0.975 * M2$nonspec_param[2L])) / length(y)


# 4. Cross-validated model: Model fit -------------------------------------

p1 <- apply_nonspecific_epistasis(s1, M1$nonspec_param)
p2 <- apply_nonspecific_epistasis(s2, M2$nonspec_param)
p3 <- apply_nonspecific_epistasis(s3, M3$nonspec_param)

plot(p1, y, pch = 16L, cex = 0.5, col = rgb(0, 0, 0, 0.2))
plot(p2, y, pch = 16L, cex = 0.5, col = rgb(0, 0, 0, 0.2))
plot(p3, y, pch = 16L, cex = 0.5, col = rgb(0, 0, 0, 0.2))

# Fraction of outliers in the residual of the second-order model.
frac_outlier_second <- sum(abs(y - p2) > 0.2 * diff(range(y))) / length(y)


# 5. Sparsity of each epistatic order -------------------------------------

v1 <- M1$e[-1L] ^ 2L / sum(M1$e[-1L] ^ 2L)
F90_1 <- min(which(cumsum(sort(v1, decreasing = TRUE)) > 0.9)) / length(v1)

v2 <- M2$e[-1L] ^ 2L / sum(M2$e[-1L] ^ 2L)
v2 <- unlist(v2[sapply(generate_site_combn(2L)[-1L], length) == 2L]) # Isolating second-order terms
F90_2 <- min(which(cumsum(sort(v2, decreasing = TRUE)) > 0.9 * sum(v2))) / length(v2)

v3 <- M3$e[-1L] ^ 2L / sum(M3$e[-1L] ^ 2L)
v3 <- unlist(v3[sapply(generate_site_combn(3L)[-1L], length) == 3L]) # Isolating site-triples
F90_3 <- min(which(cumsum(sort(v3, decreasing = TRUE)) > 0.9 * sum(v3))) / length(v3)


# 6. Sparsity of genetic architecture -------------------------------------

sparsity <- apply(read.csv('sparsity.txt', header = FALSE), 1L, mean)

F90_upper <- min(which(sparsity > 0.9)) / prod(q_list)
F90_lower <- min(which(sparsity > 0.9 * cv_third[2L])) / prod(q_list)

# Plotting.

# Model used for sparsity analysis.
model_ord <- 1L
model <- switch(as.character(model_ord), '1' = M1, '2' = M2, '3' = M3)

# Order of model terms.
e_ord <- unlist(lapply(generate_site_combn(model_ord), function(sites) if(length(sites) == 0L) 0 else rep(length(sites), prod(q_list[sites]))))

# Terms (site-combinations) ordered by contribution to variance; intercept is always included.
ord_by_contribution <- c(1L, order(model$e[-1L] ^ 2L, decreasing = TRUE) + 1L)

# Color scheme: intercept = black, first-order = magenta, second-order = cyan, third-order = green.
ord_col <- c('#000000', '#ec008b', '#00adef', '#fff100')

n_terms <- 40L
plot(1L:n_terms / prod(q_list), sparsity[1L:n_terms], type = 'l', xlim = c(0, n_terms / prod(q_list)), ylim = c(0, 1))
points(1L:n_terms / prod(q_list), sparsity[1L:n_terms], cex = 2)
points(1L:n_terms / prod(q_list), sparsity[1L:n_terms], pch = 16L, cex = 2,
       col = ord_col[e_ord[ord_by_contribution[1L:n_terms]] + 1L])
abline(h = 0.9)
abline(v = F90_upper)


# 7. Prediction from subsamples -------------------------------------------

pred <- as.matrix(read.csv('prediction.txt', header = TRUE))
frac <- vapply(colnames(pred), function(x) as.numeric(substr(x, 2L, nchar(x))), 1, USE.NAMES = FALSE)
S90_upper <- interpolate(frac, apply(pred, 2L, mean, na.rm = TRUE), 0.9)
S90_lower <- interpolate(frac, apply(pred, 2L, mean, na.rm = TRUE), 0.9 * cv_third[2L])

if(!is.finite(S90_upper)) S90_upper <- 1
if(!is.finite(S90_lower)) S90_lower <- 1

dotplot_with_sd(frac, t(pred), xlim = c(0, 0.5), ylim = c(0, 1), pch = 16L, cex = 1.5)
abline(h = 0.9)


# 8. Summary --------------------------------------------------------------

writeLines(as.character(c(prod(q_list),
                          cv_first_linear[2L],
                          cv_second_linear[2L],
                          cv_third_linear[2L],
                          cv_first[2L],
                          cv_second[2L],
                          cv_third[2L],
                          frac_out,
                          F90_upper,
                          F90_lower,
                          F90_1,
                          F90_2,
                          F90_3,
                          S90_upper, 
                          S90_lower,
                          frac_outlier_second)), 'summary.txt')
