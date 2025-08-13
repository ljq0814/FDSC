# main function: "fun_group_mst" and "fun_est_ls_mst" 
## --------- FDSC-MST algorithm -----------

########################### For estimation by ADMM algorithm #########################
#' @title ADMM algorithm
#' Input:
#' @param fdy a numeric matrix of dimension n * T, representing the functional data for n observations evaluated at T time points
#' @param timerange numeric vector of time points
#' @param B B=kronecker(diag(1,n),B.basis) B-spline bases, as the role of covariates in regression-based methods
#' @param A A=kronecker(mst_mat,diag(p))
#' @param diff_mat The second order difference operator
#' @param mst_mat (n-1)*n matrix obtained based on the constructed MST
#' @param gamma_old Initial value of gamma (gamma in our manuscript)
#' @param vl_old Initial value of vl (v_l in our manuscript)
#' @param zeta_old Initial value of zeta (zeta in our manuscript)
#' @param rho ADMM learning rate parameter rho
#' @param aa concavity parameter a of the MCP penalty
#' @param tau the regularization parameter tau refers to MCP penalty
#' @param lambda the smoothing parameter lambda
#' @param ep ADMM stopping tolerance (default = 1e-1)
#' Output: 
#' @return gamma_tmp estimated gamma by ADMM iteration (gamma in our manuscript)
#' @return vl_tmp estimated vl by ADMM iteration (v_l in our manuscript)
#####################################################################################
fun_est_ls_mst <- function(fdy, timerange, B, A, diff_mat, mst_mat, gamma_old, vl_old, zeta_old, rho=1, aa=3, lambda = 0.1, tau = 20, ep=0.1) {  
  n <- nrow(fdy)
  T <- length(timerange)
  y <- matrix(c(t(fdy)))
  tau1 = tau / rho 
  rr = 10
  p <- 6
  C <- crossprod(B) + lambda * diff_mat + rho * crossprod(A)
  chol_C <- chol(C)
  while(rr > ep) {
    gamma_tmp <- backsolve(chol_C, forwardsolve(chol_C, crossprod(B, y) + rho * t(A) %*% matrix(as.vector((t(vl_old) - 1 / rho * t(zeta_old))))))
    gamma_matrix = matrix(gamma_tmp, nrow = n, byrow = TRUE)
    delta_gamma = mst_mat %*% gamma_matrix
    ksi = delta_gamma + rho^(-1) * zeta_old
    ksinorm = apply(ksi^2, 1, sum)
    ksinorm = sqrt(ksinorm)
    SS = 1 - tau1 / ksinorm
    S = (((SS > 0) * SS) %*% t(seq(1, 1, length = p))) * ksi
    thr1 = (ksinorm > (aa * tau)) %*% t(seq(1, 1, length = p))
    thr2 = (ksinorm <= (aa * tau)) %*% t(seq(1, 1, length = p))
    tau2 = tau1 * aa / (aa - 1)
    SS2_ = 1 - tau2 / ksinorm
    S2_ = (((SS2_ > 0) * SS2_) %*% t(seq(1, 1, length = p))) * ksi
    thr3 = (ksinorm <= (tau + tau1)) %*% t(seq(1, 1, length = p))
    thr4 = ((ksinorm > (tau + tau1)) * (ksinorm <= (aa * tau))) %*% t(seq(1, 1, length = p))
      
    vl_tmp = ksi * thr1 + (S / (1 - (rho * aa)^(-1))) * thr2

    zeta_tmp = zeta_old + rho * (delta_gamma - vl_tmp)
    r = delta_gamma - vl_tmp
    r = as.vector(r)
    rr= sqrt(t(r) %*% r / length(r))
    gamma_old = gamma_tmp
    vl_old = vl_tmp
    zeta_old = zeta_tmp
  }
  list(gamma_tmp = gamma_tmp, vl_tmp = vl_tmp)
}     

##########################################################################################
########################### End  For estimation by ADMM algorithm #########################


####################################### for grouping #################################
#####################################################################################
#' Input:
#' @param fdy a numeric matrix of dimension n * T, representing the functional data for n observations evaluated at T time points
#' @param vl_new_ls Output from function "fun_est_ls_mst" $vl.tmp
#' @param edges edges of construced MST
#' Output:
#' @return group_id list of true cluster memberships
#' @return K_hat estimated number of groups
#####################################################################################
fun_group_mst = function(fdy, vl_new_ls, edges) {
  n <- nrow(fdy)
  group_w = matrix(0, n - 1, n)
  for (i in 1:(n - 1)) {
    if (sum(vl_new_ls[i, ]^2) == 0) { group_w[i, edges[i, ]] = 1 }
  }
  group_id_hat = list()
  add = list()
  i = 1
  k = 1
  group_done = setdiff(c(1:2), c(1:2))
  group_id_hat[[k]] = edges[1, 1]
  id_new = edges[1, 1]
  while (!is.na(id_new)) {
    initial_id = group_id_hat[[k]]
    ini_group = which(edges[, 1] == id_new)
    for (tt in 1:length(ini_group)) {
      group_id_hat[[k]] = union(group_id_hat[[k]], which(group_w[ini_group[tt], ] == 1))
    }
    add[[i]] = setdiff(group_id_hat[[k]], initial_id)
    while (length(add[[i]]) != 0) {
      initial_id = group_id_hat[[k]]
      group_id_temp = group_id_hat[[k]]
      for (j in 1:length(add[[i]])) {
        n_row = which(group_w[, add[[i]][j]] == 1)
        if (length(n_row) != 0) {
          for (m in 1:length(n_row)) {
            group_id_temp = union(group_id_temp, which(group_w[n_row[m], ] == 1))
          }
        }
      }
      group_id_hat[[k]] = group_id_temp
      i = i + 1
      add[[i]] = setdiff(group_id_hat[[k]], initial_id)
    }
    group_done = union(group_id_hat[[k]], group_done)
    id_new = as.vector(edges)[which(c(!edges[, 1] %in% group_done, !edges[, 2] %in% group_done))[1]]
    k = k + 1
    group_id_hat[[k]] = id_new
  }
  K_hat = k - 1
  group_id = list()
  for (kk in 1:K_hat) {
    group_id[[kk]] = group_id_hat[[kk]]
  }
  list(group_id = group_id, K_hat = K_hat)
}

########
##########################################################################################
####################################### end grouping #################################


##################### to choose the optimal taubda###############
##############################################################################
#' Input:
#' @param fdy a numeric matrix of dimension n * T, representing the functional data for n observations evaluated at T time points
#' @param B B=kronecker(diag(1,n),B.basis) B-spline bases, as the role of covariates in regression-based methods
#' @param A A=kronecker(mst_mat,diag(p))
#' @param diff_mat The second order difference operator
#' @param mst_mat (n-1)*n matrix obtained based on the constructed MST
#' @param gamma_old Initial value of gamma (gamma in our manuscript)
#' @param vl_old Initial value of vl (v_l in our manuscript)
#' @param zeta_old Initial value of zeta (zeta in our manuscript)
#' @param edges edges of construced MST
#' @param rho ADMM learning rate parameter rho
#' @param aa concavity parameter a of the MCP penalty
#' @param ep ADMM stopping tolerance (default = 1e-1)
#' Output:
#' @return tau.opt optimal regularization parameter taubda selected by BIC
#########################################################################
cv_tau_ls <- function(fdy, B, A, diff_mat,  mst_mat, gamma_old, vl_old, zeta_old, edges , rho=1, aa=3, ep = 0.1) {
  n <- nrow(fdy)
  y <- matrix(c(t(fdy)))
  tauU = 30
  tauL = 20
  taus = seq(tauL, tauU, by = 1)
  ss = length(taus)
  BIC = seq(-1, -1, length = ss)
  Qnn = seq(-1, -1, length = ss)
  for (i in 1:ss) {
    tau = taus[i]
    estimation = fun_est_ls_mst(fdy = fdy, timerange = timerange, B = B, A = A, 
                                diff_mat = diff_mat, rho = rho, aa = aa, mst_mat = mst_mat, 
                                gamma_old = gamma_old, vl_old = vl_old, zeta_old = zeta_old, tau = tau, ep = ep)
    gamma_new = estimation$gamma_tmp
    vl_new = estimation$vl_tmp
    result_group = fun_group_mst(fdy, vl_new, edges)
    group_id_hat = result_group$group_id
    K_hat = result_group$K_hat
    Error = rep(0, length(group_id_hat))
    for (k in 1:length(group_id_hat)) {
      y2_ = matrix(t(fdy[group_id_hat[[k]], ]))  
      alpha2 = solve(t(kronecker(diag(1, length(group_id_hat[[k]])), B.basis)) %*% kronecker(diag(1, length(group_id_hat[[k]])), B.basis)) %*% t(kronecker(diag(1, length(group_id_hat[[k]])), B.basis)) %*% y2_
      alpha2_matrix_ = matrix(alpha2, nrow = length(group_id_hat[[k]]), byrow = TRUE)
      alpha2_hat_ = apply(alpha2_matrix_, 2, mean)
      Error[k] = sum((fdy[group_id_hat[[k]], ] - matrix(rep(1, length(group_id_hat[[k]])), ncol = 1) %*% matrix(B.basis %*% as.matrix(alpha2_hat_), nrow = 1))^2) 
    }
    dff = K_hat * p
    Qnn[i] = sum(Error) / (n * T)
    BIC[i] = log(Qnn[i]) + log(n * p) * log(n) * dff / n
  }
  tau_opt = taus[(max(which(BIC == min(BIC))))]
  return(tau_opt)
}
#########################  End of choose the optimal taubda########################################





###########################################################
#' @title fdsc_mst
#' @description
#' Performs Functional Data Spectral Clustering via Minimum Spanning Tree (FDSC-MST).
#' Given a matrix of functional observations and optional spatial coordinates, it
#' computes a combined distance matrix (spatial + functional), constructs the MST,
#' and then solves the ADMM-based LS-MST optimization to partition the data into
#' subgroups. The output is the estimated cluster assignments and number of clusters.
#' Input:
#' @param fdy a numeric matrix of dimension n * T, representing the functional data for n observations evaluated at T time points
#' @param dx numeric vector representing the spatial x-coordinate of each observation
#' @param dy numeric vector representing the spatial y-coordinate of each observation
#' @param timerange numeric vector of time points
#' @param rho ADMM learning rate parameter rho (rho in the manuscript, typically set to 1), used in augmented Lagrangian and threshold scaling
#' @param aa concavity parameter a of the MCP penalty (default = 3), controlling shrinkage segments such as a * tau
#' @param lambda smoothing penalty parameter lambda for B-spline roughness penalty, controlling smoothness of estimated curves
#' @param tau fusion penalty parameter tau, controlling the strength of edge difference shrinkage and thus the number of clusters, usually selected by BIC
#' @param n_order B-spline order (default = 3)
#' @param el the number of B-spline basis functions
#' @param w spatial weight in [0,1] (default = 0, i.e. no spatial effect)
#' @param ep ADMM stopping tolerance (default = 1e-1)
#' Output:
#' @return group_id_hat list of estimated cluster memberships
#' @return K_hat estimated number of groups
###########################################################
fdsc_mst <- function(fdy,
                     dx = NULL, dy = NULL,
                     timerange, rho = 1, aa = 3,
                     lambda = 0.1,
                     tau = 20,
                     n_order = 3,
                     el       = 6,
                     w       = 0,
                     ep      = 10^(-1)) {
  n <- nrow(fdy)
  T <- length(timerange)
  y <- matrix(c(t(fdy)))
  n_knots <- el - n_order
  breaks <- seq(min(timerange), max(timerange), length = n_knots + 2)
  B.basis <- bsplineS(timerange, breaks, n_order)       
  B <- kronecker(diag(1, n), B.basis)
  D_matrix = matrix(0, el - 2, el)
  for (i in 1:(el - 2)) {
    D_matrix[i, i] = 1
    D_matrix[i, i + 2] = 1
    D_matrix[i, i + 1] = -2
  }
  D2 = t(D_matrix) %*% D_matrix
  diff_mat = kronecker(diag(1, n), D2)
  distance_matrix <- matrix(0, n, n)
  for(i in 1:n) {
    for(j in 1:n) {
      if(!is.null(dx) && !is.null(dy)) {
        sp_dist <- sqrt((dx[i] - dx[j])^2 / 2 + (dy[i] - dy[j])^2 / 2)
      } else {
        sp_dist <- 0
      }
      func_dist <- sqrt(mean((fdy[i, ] - fdy[j, ])^2))
      distance_matrix[i, j] <- w * sp_dist + (1 - w) * func_dist
    }
  }
  graph_cut <- graph_from_adjacency_matrix(distance_matrix, weighted = TRUE, mode = "undirected")
  mst_cut <- mst(graph_cut)
  edges_cut <- as_edgelist(mst_cut)
  mst_mat_mst <- matrix(0, n - 1, n)
  for (i in seq_len(n - 1)) {
    mst_mat_mst[i, edges_cut[i, 1]] <-  1
    mst_mat_mst[i, edges_cut[i, 2]] <- -1
  }
  A_mst <- kronecker(mst_mat_mst, diag(el))
  C_init <- crossprod(B) + 0.1 * kronecker(diag(n), D2) + 0.001 * crossprod(A_mst)
  chol_C_init <- chol(C_init)
  gamma_initial <- backsolve(chol_C_init, forwardsolve(chol_C_init, crossprod(B, y)))
  vl_initial <- mst_mat_mst %*% matrix(gamma_initial, nrow = n, byrow = TRUE)
  zeta_initial <- matrix(0, n - 1, el)
  est <- fun_est_ls_mst(fdy = fdy, timerange = timerange, B = B, A = A_mst, diff_mat = diff_mat, 
                         rho = rho, aa = aa, mst_mat = mst_mat_mst,
                         gamma_old = gamma_initial, vl_old = vl_initial, zeta_old = zeta_initial,
                         tau = tau, lambda = lambda, ep = 0.1)
  grp <- fun_group_mst(fdy = fdy, vl_new_ls = est$vl_tmp, edges = edges_cut)
  list(
    group_id_hat = grp$group_id,
    K_hat    = grp$K_hat
  )
}

################################################################################
########## End FDSC-MST algorithm ##################################################
