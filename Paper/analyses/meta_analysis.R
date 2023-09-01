
# This script summarizes the analyses for individual datasets and generates plots.


# Libraries, sources, and imports -----------------------------------------

library(ggplot2)

setwd('~/Desktop/GP2/analyses/')

# Dataset name; in the order of the number of possible genotypes.
dataset_list <- c('MPH', 'beta-lactamase', 'DHFR', 'hemagglutinin', 'CR6261-H1', 'CR6261-H9',
                  'ParD3-ParE2 Lite', 'ParD3-ParE3 Lite', 'avGFP', 'ParD3-ParE2 Aarke',
                  'ParD3-ParE3 Aarke', 'spike', 'CH65-G189E', 'CH65-MA90', 'CH65-SI06',
                  'CR9114-B', 'CR9114-H1', 'CR9114-H3', 'ParB', 'GB1')

# Number of states per site; 1 = two, 2 = intermediate, 3 = 20
data_type <- c(1, 1, 2, 2, 1, 1, 3, 3, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 3, 3) 

# Importing summary files for individual datasets.
summary <- data.frame(t(sapply(dataset_list, function(dataset) 
  as.numeric(readLines(paste0(dataset, '/summary.txt'))))))

colnames(summary) <- c('N', 'cv_1_lin', 'cv_2_lin', 'cv_3_lin', 'cv_1', 'cv_2', 'cv_3',
                       'frac_out', 'F90_upper', 'F90_lower', 'F90_1', 'F90_2', 'F90_3',
                       'S90_upper', 'S90_lower', 'frac_outlier_second')

# N : Number of possible genotypes
# cv_x_lin : Out-of-sample R2 without modeling nonspecific epistasis; order x
# cv_x : Out-of-sample R2 when nonspecific epistasis is modeled; order x
# frac_out : Fraction of genotypes outside the phenotype bounds
# F90_upper, F90_lower : Upper and lower estimate of F90, which is T90 multiplied by N
# F90_x : F90 calculated for each order of terms, instead of all terms
# S90_upper, S90_lower : Upper and lower estimate of S90, which is N90 dividied by N
# frac_outlier_second : Fraction of outliers in the second-order model.


# Reference-free analysis -------------------------------------------------

# Fig. 2A.

data <- data.frame(order = as.factor(c(rep(c(1L, 2L, 3L), each = 20L))),
                   R2 = c(summary$cv_1, summary$cv_2, summary$cv_3))

ggplot(data, aes(x = order, y = R2)) +
  geom_boxplot(coef = 100) +
  ylim(0, 1) +
  geom_jitter(color="black", size = 4, alpha = 0.5, width = 0.1) +
  geom_hline(yintercept = 0.9)

# Fig. 2B, right panel.

frac_outlier_second <- log(summary$frac_outlier_second, 10)
frac_outlier_second[frac_outlier_second < -4] <- -5
plot(rep(1, nrow(summary)), frac_outlier_second,
     ylim = c(-5.25, 0), xlim = c(0.75, 1.25), cex = 2)


# Reference-based analysis ------------------------------------------------

# Importing data.
dir_list <- vapply(dataset_list, function(dataset) paste0(dataset, '/ref_based_analysis.txt'), character(1L))

# R2 for higher-order mutants, based on the wild-type genotype as reference.
res_WT <- lapply(dir_list, function(dir) {res <- read.csv(dir); res[1L, -1L]})

# R2 for higher-order mutants, based on randomly chosen references.
res_random <- lapply(dir_list, function(dir) {res <- read.csv(dir); res[-1L, ]})

# Fig. 2C.

data <- data.frame(order = as.factor(unlist(lapply(res_WT, seq_along))), R2 = unlist(res_WT))
data$R2[data$R2 < 0] <- 0

ggplot(data, aes(x = order, y = R2)) +
  geom_boxplot(coef = 100) +
  ylim(0, 1) +
  geom_jitter(color="black", size = 3, alpha = 0.5, width = 0.05) +
  geom_hline(yintercept = 0.9)

# Fig. S2.

plot(NA_real_, xlim = c(1, 3 * length(dataset_list)), ylim = c(0, 1))
for(i in seq_along(dataset_list)) {
  
  R2_random <- as.matrix(res_random[[i]][, -1L]) # R2 based on random references
  R2_random[R2_random < 0] <- 0 # Negative values plotted as 0
  
  # First-order model.
  points(rep(3L * (i - 1L) + 1L, 5L),
         quantile(R2_random[, 1L], c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE), pch = 3L) # Random
  points(3L * (i - 1L) + 1L, max(res_WT[[i]][1L], 0), pch = 16L) # WT
  
  # Second-order model.
  points(rep(3L * (i - 1L) + 2L, 5L),
         quantile(R2_random[, 2L], c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE), pch = 3L) # Random
  points(3L * (i - 1L) + 2L, max(res_WT[[i]][2L], 0), pch = 16L) # WT
  
  # Third-order model, if applicable.
  if(ncol(R2_random) == 3L) {
    points(rep(3L * (i - 1L) + 3L, 5L),
           quantile(R2_random[, 3L], c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE), pch = 3L) # Random
    points(3L * (i - 1L) + 3L, max(res_WT[[i]][3L], 0), pch = 16L) # WT
  }
}


# Impact of modeling nonspecific epistasis --------------------------------

# Fig. 3A.

data <- data.frame(order = as.factor(c(rep(c(1L, 2L, 3L), each = 20L))),
                   R2 = c(summary$cv_1_lin, summary$cv_2_lin, summary$cv_3_lin))

ggplot(data, aes(x = order, y = R2)) +
  geom_boxplot(coef = 100) +
  ylim(0, 1) +
  geom_jitter(color="black", size = 4, alpha = 0.5, width = 0.05) +
  geom_hline(yintercept = 0.9)

# Fig. 3B.

plot(summary$cv_1_lin, summary$cv_1, xlim = c(0, 1), ylim = c(0, 1), cex = 2.7) # Main effects
abline(a = 0, b = 1)

plot(summary$cv_2_lin - summary$cv_1_lin, summary$cv_2 - summary$cv_1, # Pairwise interactions
     xlim = c(0, 1), ylim = c(0, 1), cex = 2.7)
abline(a = 0, b = 1)

plot(summary$cv_3_lin - summary$cv_2_lin, summary$cv_3 - summary$cv_2, # Three-way interactions
     xlim = c(0, 1), ylim = c(0, 1), cex = 2.7)
abline(a = 0, b = 1)

# Fig. 3C.

plot(log(1 - summary$frac_out, 10), summary$cv_1 - summary$cv_1_lin,
     xlim = c(-3, 0), ylim = c(0, 1), cex = 3, pch = 16L) # First-order
cor(log(1 - summary$frac_out, 10), summary$cv_1 - summary$cv_1_lin, method = 'spearman')

plot(log(1 - summary$frac_out, 10), summary$cv_2 - summary$cv_2_lin,
     xlim = c(-3, 0), ylim = c(0, 1), cex = 3, pch = 16L) # Second-order
cor(log(1 - summary$frac_out, 10), summary$cv_2 - summary$cv_2_lin, method = 'spearman')

plot(log(1 - summary$frac_out, 10), summary$cv_3 - summary$cv_3_lin,
     xlim = c(-3, 0), ylim = c(0, 1), cex = 3, pch = 16L) # Thrd-order
cor(log(1 - summary$frac_out, 10), summary$cv_3 - summary$cv_3_lin,
    method = 'spearman', use = 'pairwise.complete')


# Sparsity (T90) ----------------------------------------------------------

T90_upper <- summary$N * summary$F90_upper
T90_lower <- summary$N * summary$F90_lower
T90 <- (T90_upper + T90_lower) / 2

# Fig. 4B.

pal <- c('#33a02c', '#b15928', '#ee82b4') # Coloring by number of states.

plot(log(summary$N, 10), log(T90, 10), xlim = c(0.5, 5.5), ylim = c(0.5, 5.5),
     cex = 2.8, pch = 16L, col = pal[data_type])
points(log(summary$N, 10), log(T90_upper, 10), cex = 3, pch = 3L, col = pal[data_type])
points(log(summary$N, 10), log(T90_lower, 10), cex = 3, pch = 3L, col = pal[data_type])
abline(a = 0, b = 1)

l <- lm(log(T90, 10) ~ log(summary$N, 10))
abline(a = l$coefficients[1L], b = l$coefficients[2L])

# Fig. 4C.

pal <- c('#33a02c', '#b15928', '#ee82b4') # Coloring by number of states.

plot(log(summary$N, 10), log((summary$F90_lower + summary$F90_upper) / 2, 10),
     xlim = c(0.5, 5.5), ylim = c(-4.5, 0), cex = 2.8, pch = 16L, col = pal[data_type])
points(log(summary$N, 10), log(summary$F90_upper, 10), cex = 3, pch = 3L, col = pal[data_type])
points(log(summary$N, 10), log(summary$F90_lower, 10), cex = 3, pch = 3L, col = pal[data_type])

cor(log((summary$F90_lower + summary$F90_upper) / 2, 10), log(summary$N, 10), method = 'spearman')


# Prediction from subsample -----------------------------------------------

# Fig. 5B.

N90_upper <- summary$N * summary$S90_upper
N90_lower <- summary$N * summary$S90_lower

# Excluding the three small datasets
N90_upper <- N90_upper[-c(1, 2, 3)]
N90_lower <- N90_lower[-c(1, 2, 3)]
N90 <- (N90_upper + N90_lower) / 2
N <- summary$N[-c(1, 2, 3)]

plot(log(N, 10), log(N90, 10), xlim = c(0.5, 5.5), ylim = c(0.5, 5.5), cex = 2.4)
points(log(N, 10), log(N90_upper, 10), cex = 3, pch = 3L)
points(log(N, 10), log(N90_lower, 10), cex = 3, pch = 3L)
abline(a = 0, b = 1)

plot(rep(1, length(N90)), log(N90 / N, 10), ylim = c(-3, 0))

# Fig. 5C.

plot(log(T90[-c(1, 2, 3)], 10), log(N90, 10), xlim = c(0.5, 5.5), ylim = c(0.5, 5.5), cex = 2.4)
abline(a = 0, b = 1)

# Fig. 5D.

plot(log(1 - summary$frac_out[-c(1, 2, 3)], 10), log(N90, 10),
     xlim = c(-3, 0), ylim = c(0.5, 5.5), cex = 2)
cor(log(1 - summary$frac_out[-c(1, 2, 3)], 10), log(N90, 10), method = 'spearman')

# Fig. 5E.

l <- lm(log(N90, 10) ~ log(N, 10) + log(1 - summary$frac_out[-c(1, 2, 3)], 10))

plot(l$fitted.values, log(N90, 10), xlim = c(0.5, 5.5), ylim = c(0.5, 5.5), cex = 2)
abline(a = 0, b = 1)


# Simulated data ----------------------------------------------------------

# Fig. 1D.

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

# Simulated map
n <- 8L
s <- 2L
ord <- identify_genotype_order(n, s)

# Reference-based effects.
lambda <- as.numeric(readLines('simulated/2_8/ref_based_analysis.txt'))

# Background-averaged effects.
b <- as.numeric(readLines('simulated/2_8/bg_avg_analysis.txt'))

# Reference-free effects.
e <- as.numeric(readLines('simulated/2_8/ref_free_analysis.txt'))

# Plotting.

spacing <- 2

plot(NA, xlim = c(-0.5, n * spacing + 0.5), ylim = c(-max(c(lambda, b, e)), max(c(lambda, b, e))))
points(spacing * ord - 0.5, e, pch = 16L, col = rgb(0, 0, 0, 0.5), cex = 1.5)
points(spacing * ord, b, pch = 16L, col = rgb(0, 1, 0, 0.5), cex = 1.5)
points(spacing * ord + 0.5, lambda, pch = 16L, col = rgb(1, 0, 0, 0.5), cex = 1.5)
abline(h = c(-1, 0, 1))
abline(v = spacing * 1L:n - 1)

# Standard deviations.
for(i in 1L:(n - 1L)) {
  
  points(spacing * i - 0.5, mean(e[ord == i]), pch = 3L)  
  arrows(x0 = spacing * i - 0.5, y0 = mean(e[ord == i]) - sd(e[ord == i]),
         x1 = spacing * i - 0.5, y1 = mean(e[ord == i]) + sd(e[ord == i]), angle = 90, code = 3, length = 0.05)
  
  points(spacing * i, mean(b[ord == i]), pch = 3L)  
  arrows(x0 = spacing * i, y0 = mean(b[ord == i]) - sd(b[ord == i]),
         x1 = spacing * i, y1 = mean(b[ord == i]) + sd(b[ord == i]),
         angle = 90, code = 3, length = 0.05)
  
  points(spacing * i + 0.5, mean(lambda[ord == i]), pch = 3L)  
  arrows(x0 = spacing * i + 0.5, y0 = mean(lambda[ord == i]) - sd(lambda[ord == i]),
         x1 = spacing * i + 0.5, y1 = mean(lambda[ord == i]) + sd(lambda[ord == i]),
         angle = 90, code = 3, length = 0.05)
}


# Fig. 1E; shown for the 20-state simulations.

frac <- seq(-7, -1, by = 1)

bg_avg <- read.csv('simulated/20_3_1/bg_avg_prediction.txt', header = TRUE)
ref_free <- read.csv('simulated/20_3_1/ref_free_prediction.txt', header = TRUE)

d <- data.frame(frac = as.factor(c(rep(frac, each = 250L), rep(frac, each = 200L))),
                form = c(rep('b', 250L * length(frac)), rep('r', 200L * length(frac))),
                R2 = unname(c(unlist(bg_avg), unlist(ref_free))))

d$R2[d$R2 < 0] <- 0 # Negative R2 plotted as 0

# Boxplot with whiskers extending to the data limits.
ggplot(d, aes(x = frac, y = R2, fill = form)) + geom_boxplot(coef = 100)


# For Fig. 1F, see the /analysis/simulated/10_4_only_3rd_nonspec


# Miscellaneous plots -----------------------------------------------------

# Fig. 1B.

t <- seq(-6.25, 6.25, by = 0.01)
plot(t, 0.5 + 2.5 / (1 + exp(-t)), type = 'l', ylim = c(0.3, 3.2))
abline(h = c(0.5, 3))

# Fig. 6A.

t <- seq(-2.7, 2.7, by = 0.01)
plot(t, 1 / (1 + exp(-t)), type = 'l', ylim = c(0, 1))
abline(v = c(-2.3, -1.15, 0, 1.15, 2.3))
abline(h = c(0, 1))
