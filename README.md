Functional Data Subgroup Clustering Algorithm via Minimum Spanning Tree
================

## Overview

`FDSC` implements two functional data clustering algorithms based on
minimum spanning tree (MST).

**FDSC-MST** addresses subgroup clustering by constructing an MST using
both geographical distance and functional curve similarity, and applying
subgroup pairwise penalties to MST edges. The optimization problem is
solved via the alternating direction method of multipliers (ADMM).

**FDSCOD-MST** extends this framework to contaminated data by
simultaneously performing clustering and outlier detection.

The package provides the following main functions:

- `fdsc_mst()`: Functional Data Subgroup Clustering via MST (FDSC-MST)
- `fdscod_mst()`: Functional Data Subgroup Clustering and Outlier
  Detection via MST (FDSCOD-MST)

## Installation

``` r
# install.packages("devtools")
devtools::install_github("LC8736/FDSC")
```

## Examples

### Example 1: Clustering without Outliers (FDSC-MST)

We simulate $n = 200$ functional observations over $t \in [1, 21]$ with
step size 0.2, under two clusters $G_1$ and $G_2$ with mean functions:

- For $i \in G_1$, $\mu_i(t) = \max(6 - |t - 7|, 0)$
- For $i \in G_2$, $\mu_i(t) = \max(6 - |t - 15|, 0)$

Spatial locations $(h_{xi},h_{yi})$ are sampled uniformly from
$[0, 1]^2$, which satisfy:

$$G_1 = \lbrace i \mid h_{xi} > h_{yi}\rbrace, \quad G_2 = \lbrace i \mid h_{xi} \leq h_{yi}\rbrace$$

Each observation is a smooth curve with Gaussian noise of standard
deviation 0.1.

``` r
library(FDSC)
library(flexclust) ## Compute ARI
set.seed(123)

# Simulation setup
n     <- 200          # number of observations
K     <- 2            # number of true clusters
grid  <- seq(1, 21, length = 101)  # time grid of length 101
sig_x <- 0.1          # noise level

# Generate data: two clusters separated by the diagonal y = x in [0,1]^2
dat1 <- datageneration(n = n, grid = grid, sig_x = sig_x,
                        K = K, spatial_dist_type = "uniform",
                        is_outlier = FALSE)
# dat1$Y        : 200 x 101 matrix of functional observations
# dat1$coord    : 200 x 2 matrix of spatial coordinates
# dat1$membership : list of true cluster memberships

# Run FDSC-MST
res_fdscmst <- fdsc_mst(Y = dat1$Y, grid = grid, coord = dat1$coord,
                    lambda = 0.1, tau = 20, spline_order = 3,
                    n_basis = 6, sp_weight = 0, tol = 0.1)
summary(res_fdscmst)

# Evaluate clustering performance
# ARI: Adjusted Rand Index (1 = perfect clustering, 0 = random)
label_true <- integer(n)
for (k in 1:K) label_true[dat1$membership[[k]]] <- k
cat(flexclust::comPart(label_true, res_fdscmst$label, type = "ARI")) ## ARI calculation, close to 1 indicates accurate clustering
```

### Example 2: Clustering with Outlier Detection (FDSCOD-MST)

We introduce $n_{\text{out}} = 10$ outliers into the dataset from
Example 1. The outlier mean functions follow the contamination scheme:

$$\mu_i(t) = \theta\lbrace\delta \min(|t - 7| - 6, 0) + (1 - \delta)\min(|t - 15| - 6, 0)\rbrace, \quad i \in \mathcal{O},$$

where $\delta \sim \text{Bernoulli}(0.5)$ and $\theta = 0.5$.

``` r
set.seed(123)
n <- 200; K <- 2; grid <- seq(1, 21, length = 101); sig_x <- 0.1
theta <- 0.5   # outlier contamination intensity

# Generate contaminated data: 200 observations, 10 of which are outliers
dat2 <- datageneration(n = n, grid = grid, sig_x = sig_x,
                        K = K, spatial_dist_type = "uniform",
                        is_outlier = TRUE, outliercase = "case1",
                        theta = theta)
# dat2$outlier : indices of the 10 outliers

# Run FDSCOD-MST
res_fdscodmst <- fdscod_mst(Y = dat2$Y, coord = dat2$coord, grid = grid,
                      lambda = 0.1, tau = 40, spline_order = 3, n_basis = 6,
                      sp_weight = 0, tol = 0.1, alpha = 0.05)
summary(res_fdscodmst)

# Evaluate clustering and outlier detection performance
# ARI : clustering accuracy
# TPR : True Positive Rate (proportion of true outliers correctly detected)
# FPR : False Positive Rate (proportion of non-outliers falsely detected)
label_true2 <- integer(n)
for (k in 1:K) label_true2[dat2$membership[[k]]] <- k
cat(flexclust::comPart(label_true2, res_fdscodmst$label, type = "ARI")) ## ARI calculation
TPR <- length(intersect(res_fdscodmst$outlier, dat2$outlier)) / length(dat2$outlier) ## TPR calculation
cat(TPR)
FPR <- length(setdiff(res_fdscodmst$outlier, dat2$outlier)) / (n - length(dat2$outlier)) ## FPR calculation
cat(FPR)
```

## Reference

Liu, C., Lin, S., and Li, J. (2026). Functional data clustering and
outlier identification via subgroup analysis with minimum spanning tree.
*Information Processing & Management*, 63(5), 104646.
<https://doi.org/10.1016/j.ipm.2026.104646>
