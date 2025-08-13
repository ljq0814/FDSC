###  Example:
#########################################################################################
###  Simulation setup
###  datageneration 
library(MASS)
library(fda)
library(fdapace)
library(igraph)
library(flexclust)
source("./R/datageneration.R")
source("./R/FDSC-MST.R")
source("./R/FDSCOD-MST.R")
source("./R/utils.R")
n_ = 200
K_ = 2
timerange_ = seq(1,21,length = 101)
sig_x_ = 0.1
theta_ = 0.5
data1 <- datageneration(n = n_, timerange = timerange_, sig_x = sig_x_, K = K_, spatial_dist_type = "uniform", is_outlier = FALSE, outliercase = NULL)
data2 <- datageneration(n = n_, timerange = timerange_, sig_x = sig_x_, K = K_, spatial_dist_type = "uniform", is_outlier = TRUE, outliercase = "case1", theta = theta_)

### FDSC-MST
result1 <- fdsc_mst(fdy = data1$fdy, dx = data1$dx, dy = data1$dy, timerange = timerange_, lambda = 0.1, tau = 20, n_order = 3, el = 6, w = 0, ep = 0.1)

####### Evaluation : by RI for clustering accuracy; MSE for estimation accuracy
evaluation1 <- evaluation_fun(group_id = data1$group_id, group_id_hat = result1$group_id_hat, K = K_, K_hat = result1$K_hat, fdy = data1$fdy,
                              mu_true = data1$mu_true, is_outlier = FALSE, out_true = NULL, out_set_final = NULL)
print(evaluation1)

### FDSCOD-MST
result2 <- fdscod_mst(fdy = data2$fdy, dx = data2$dx, dy = data2$dy, timerange = timerange_, rho = 1, lambda = 0.1, tau = 40, n_order = 3, 
                      el = 6, w = 0, ep = 0.1, alpha0 = 0.05)
evaluation2 <- evaluation_fun(group_id = data2$group_id, group_id_hat=result2$group_id_final, K=K_, K_hat=result2$K_hat_final,
                              fdy = data2$fdy, mu_true = data2$mu_true, is_outlier = TRUE, out_true = data2$out_true, out_set_final=result2$out_set_final)
print(evaluation2)
############ End Evaluation for simulation studies ######################################
