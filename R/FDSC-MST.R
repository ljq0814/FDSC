# main function: "recover_clusters" and "admm_mst" 
## --------- FDSC-MST algorithm -----------

#' ADMM Estimation for FDSC-MST
#'
#' Solves the ADMM optimization problem for functional data subgroup clustering
#' via minimum spanning tree, iteratively updating the spline coefficient vector
#' \eqn{\gamma}, the auxiliary variable \eqn{v_l}, and the dual variable
#' \eqn{\zeta} until convergence or the maximum number of iterations is reached.
#'
#' @param Y A numeric matrix of dimension \code{n x T}, representing
#'   functional data for \code{n} observations evaluated at \code{T} time
#'   points.
#' @param grid A numeric vector of time points of length \code{T}.
#' @param B A numeric matrix equal to \code{kronecker(diag(1, n), B_basis)},
#'   representing the B-spline basis expanded for all observations.
#' @param A A numeric matrix equal to \code{kronecker(mst_mat, diag(p))}.
#' @param diff_mat The second-order difference penalty matrix.
#' @param mst_mat An \code{(n-1) x n} matrix derived from the constructed MST.
#' @param gamma_old Initial value of the spline coefficient vector
#'   \eqn{\gamma}.
#' @param vl_old Initial value of the auxiliary variable \eqn{v_l}.
#' @param zeta_old Initial value of the dual variable \eqn{\zeta}.
#' @param rho A positive numeric scalar specifying the ADMM learning rate.
#'   Defaults to \code{1}.
#' @param a A positive numeric scalar specifying the concavity parameter of
#'   the MCP penalty. Defaults to \code{3}.
#' @param lambda A positive numeric scalar specifying the smoothing penalty
#'   parameter for the B-spline roughness penalty. Defaults to \code{0.1}.
#' @param tau A positive numeric scalar specifying the regularization parameter
#'   of the MCP penalty. Defaults to \code{20}.
#' @param tol A positive numeric scalar specifying the ADMM stopping tolerance.
#'   Defaults to \code{0.1}.
#' @param max_iter A positive integer specifying the maximum number of ADMM
#'   iterations. Defaults to \code{10000}.
#'
#' @return A list with two components:
#'   \describe{
#'     \item{gamma}{Estimated spline coefficient vector \eqn{\gamma} at
#'       convergence.}
#'     \item{vl}{Estimated auxiliary variable \eqn{v_l} at convergence.}
#'   }
#'
#' @keywords internal
admm_mst <- function(Y, grid, B, A, diff_mat, mst_mat, gamma_old, vl_old, zeta_old, rho=1, a=3, lambda = 0.1, tau = 20, tol=0.1, max_iter = 10000) {  
  n <- nrow(Y)
  n_time <- length(grid)
  y <- matrix(c(t(Y)))
  tau1 <- tau / rho 
  rr <- 10
  iter <- 0
  p <- ncol(B)/n
  C <- crossprod(B) + lambda * diff_mat + rho * crossprod(A)
  chol_C <- chol(C)
  while(rr > tol && iter < max_iter) {
    gamma_tmp <- backsolve(chol_C, forwardsolve(chol_C, crossprod(B, y) + rho * t(A) %*% matrix(as.vector((t(vl_old) - 1 / rho * t(zeta_old))))))
    gamma_matrix <- matrix(gamma_tmp, nrow = n, byrow = TRUE)
    delta_gamma <- mst_mat %*% gamma_matrix
    ksi <- delta_gamma + rho^(-1) * zeta_old
    ksinorm <- apply(ksi^2, 1, sum)
    ksinorm <- sqrt(ksinorm)
    SS <- 1 - tau1 / ksinorm
    S <- (((SS > 0) * SS) %*% t(seq(1, 1, length = p))) * ksi
    thr1 <- (ksinorm > (a * tau)) %*% t(seq(1, 1, length = p))
    thr2 <- (ksinorm <= (a * tau)) %*% t(seq(1, 1, length = p))
    tau2 <- tau1 * a / (a - 1)
    SS2_ <- 1 - tau2 / ksinorm
    S2_ <- (((SS2_ > 0) * SS2_) %*% t(seq(1, 1, length = p))) * ksi
    thr3 <- (ksinorm <= (tau + tau1)) %*% t(seq(1, 1, length = p))
    thr4 <- ((ksinorm > (tau + tau1)) * (ksinorm <= (a * tau))) %*% t(seq(1, 1, length = p))
      
    vl_tmp <- ksi * thr1 + (S / (1 - (rho * a)^(-1))) * thr2

    zeta_tmp <- zeta_old + rho * (delta_gamma - vl_tmp)
    r <- delta_gamma - vl_tmp
    r <- as.vector(r)
    rr <- sqrt(t(r) %*% r / length(r))
    gamma_old <- gamma_tmp
    vl_old <- vl_tmp
    zeta_old <- zeta_tmp
    iter <- iter + 1
  }
  if (iter == max_iter)
    warning("ADMM did not converge within ", max_iter, " iterations")
  list(gamma = gamma_tmp, vl = vl_tmp)
}     

##########################################################################################
########################### End  For estimation by ADMM algorithm #########################


#' Recover Cluster Memberships from MST Edge Penalties
#'
#' Recovers cluster memberships from the estimated auxiliary variable \code{vl}
#' by identifying zero-penalty edges in the MST and propagating connected
#' components to assign observations to clusters.
#'
#' @param Y A numeric matrix of dimension \code{n x T}, representing
#'   functional data for \code{n} observations evaluated at \code{T} time
#'   points.
#' @param vl The estimated auxiliary variable \eqn{v_l}, as returned by
#'   \code{admm_mst()} in the \code{vl} component. An \code{(n-1) x p}
#'   matrix where a zero row indicates the corresponding MST edge is not
#'   cut.
#' @param edges An \code{(n-1) x 2} integer matrix giving the edge list of
#'   the constructed MST, as returned by \code{igraph::as_edgelist()}.
#'
#' @return A list with two components:
#'   \describe{
#'     \item{membership}{A list of integer vectors, each giving the observation
#'       indices belonging to one estimated cluster.}
#'     \item{cluster}{A positive integer giving the estimated number of
#'       clusters.}
#'   }
#'
#' @keywords internal
recover_clusters <- function(Y, vl, edges) {
  n <- nrow(Y)
  group_w = matrix(0, n - 1, n)
  for (i in 1:(n - 1)) {
    if (sum(vl[i, ]^2) == 0) { group_w[i, edges[i, ]] = 1 }
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
  list(membership = group_id, cluster = K_hat)
}

########
##########################################################################################
####################################### end grouping #################################


##################### to choose the optimal tau ###############
#' BIC-Based Selection of Regularization Parameter
#'
#' Selects the optimal regularization parameter \eqn{\tau} for FDSC-MST by
#' minimizing a BIC criterion over a grid of candidate values, with
#' \eqn{\lambda} fixed at the provided value.
#'
#' @param Y A numeric matrix of dimension \code{n x T}, representing
#'   functional data for \code{n} observations evaluated at \code{T} time
#'   points.
#' @param grid A numeric vector of time points of length \code{T}.
#' @param B_basis A numeric matrix of dimension \code{T x n_basis} giving the
#'   B-spline basis evaluated at the time points in \code{grid}.
#' @param mst_mat An \code{(n-1) x n} matrix derived from the constructed MST.
#' @param gamma_old Initial value of the spline coefficient vector \eqn{\gamma}.
#' @param vl_old Initial value of the auxiliary variable \eqn{v_l}.
#' @param zeta_old Initial value of the dual variable \eqn{\zeta}.
#' @param edges An \code{(n-1) x 2} integer matrix giving the edge list of
#'   the constructed MST.
#' @param rho A positive numeric scalar specifying the ADMM learning rate.
#'   Defaults to \code{1}.
#' @param a A positive numeric scalar specifying the concavity parameter of
#'   the MCP penalty. Defaults to \code{3}.
#' @param lambda A positive numeric scalar specifying the smoothing penalty
#'   parameter, fixed during the search over \eqn{\tau}. Defaults to
#'   \code{0.1}.
#' @param tol A positive numeric scalar specifying the ADMM stopping tolerance.
#'   Defaults to \code{0.1}.
#' @param max_iter A positive integer specifying the maximum number of ADMM
#'   iterations. Defaults to \code{10000}.
#' @param tau_lower A positive numeric scalar specifying the lower bound of the
#'   search range for \eqn{\tau}. Defaults to \code{20}.
#' @param tau_upper A positive numeric scalar specifying the upper bound of the
#'   search range for \eqn{\tau}. Defaults to \code{30}.
#' @param tau_step A positive numeric scalar specifying the step size of the
#'   search grid for \eqn{\tau}. Defaults to \code{1}.
#'
#' @return A numeric scalar giving the optimal regularization parameter
#'   \eqn{\tau} selected by BIC.
#'
#' @export
bic_tau <- function(Y, grid, B_basis, mst_mat, gamma_old, vl_old, zeta_old, edges, 
                      rho=1, a=3, lambda = 0.1, tol = 0.1, max_iter = 10000, 
                      tau_lower = 20, tau_upper = 30, tau_step = 1) {
  n <- nrow(Y)
  y <- matrix(c(t(Y)))
  B <- kronecker(diag(1, n), B_basis)
  p <- ncol(B) / n
  A <- kronecker(mst_mat,diag(p))
  
  D_matrix = matrix(0, p - 2, p)
  for (i in 1:(p - 2)) {
    D_matrix[i, i] = 1
    D_matrix[i, i + 2] = 1
    D_matrix[i, i + 1] = -2
  }
  D2 = t(D_matrix) %*% D_matrix
  diff_mat = kronecker(diag(1, n), D2)

  n_time <- length(grid)
  taus <- seq(tau_lower, tau_upper, by = tau_step)
  ss = length(taus)
  BIC = seq(-1, -1, length = ss)
  Qnn = seq(-1, -1, length = ss)
  for (i in 1:ss) {
    tau = taus[i]
    estimation = admm_mst(Y = Y, grid = grid, B = B, A = A, 
                          diff_mat = diff_mat, rho = rho, a = a, mst_mat = mst_mat, 
                          gamma_old = gamma_old, vl_old = vl_old, zeta_old = zeta_old, 
                          lambda = lambda, tau = tau, tol = tol, max_iter = max_iter)
    gamma_new = estimation$gamma
    vl_new = estimation$vl
    result_group = recover_clusters(Y, vl_new, edges)
    group_id_hat = result_group$membership
    K_hat = result_group$cluster
    Error = rep(0, length(group_id_hat))
    for (k in 1:length(group_id_hat)) {
      y2_ = matrix(t(Y[group_id_hat[[k]], ]))  
      alpha2 = solve(t(kronecker(diag(1, length(group_id_hat[[k]])), B_basis)) %*% kronecker(diag(1, length(group_id_hat[[k]])), B_basis)) %*% t(kronecker(diag(1, length(group_id_hat[[k]])), B_basis)) %*% y2_
      alpha2_matrix_ = matrix(alpha2, nrow = length(group_id_hat[[k]]), byrow = TRUE)
      alpha2_hat_ = apply(alpha2_matrix_, 2, mean)
      Error[k] = sum((Y[group_id_hat[[k]], ] - matrix(rep(1, length(group_id_hat[[k]])), ncol = 1) %*% matrix(B_basis %*% as.matrix(alpha2_hat_), nrow = 1))^2) 
    }
    dff = K_hat * p
    Qnn[i] = sum(Error) / (n * n_time)
    BIC[i] = log(Qnn[i]) + log(n * p) * log(n) * dff / n
  }
  tau_opt = taus[(max(which(BIC == min(BIC))))]
  return(tau_opt)
}
#########################  End of choose the optimal taubda########################################


#' Functional Data Subgroup Clustering via Minimum Spanning Tree (FDSC-MST)
#'
#' Performs functional data subgroup clustering by constructing a minimum
#' spanning tree (MST) from a combined spatial and functional distance matrix,
#' and solving an ADMM-based optimization problem with MCP subgroup penalties
#' applied to MST edges.
#'
#' @param Y A numeric matrix of dimension \code{n x T}, representing
#'   functional data for \code{n} observations evaluated at \code{T} time
#'   points.
#' @param grid A numeric vector of time points of length \code{T}.
#' @param coord An optional \code{n x 2} numeric matrix of spatial coordinates,
#'   where each row gives the \code{(x, y)} coordinates of one observation.
#'   Set to \code{NULL} (default) if no spatial information is available.
#' @param rho A positive numeric scalar specifying the ADMM learning rate.
#'   Defaults to \code{1}.
#' @param a A positive numeric scalar specifying the nonconvexity of
#'   the MCP penalty. Defaults to \code{3}.
#' @param lambda A positive numeric scalar specifying the smoothing penalty
#'   parameter for the B-spline roughness penalty. Defaults to \code{0.1}.
#' @param tau A positive numeric scalar specifying the regularization parameter
#'   of the MCP penalty, controlling the shrinkage and thus the number of clusters. 
#'   Usually selected by BIC. Defaults to \code{20}.
#' @param spline_order A positive integer specifying the B-spline order.
#'   Defaults to \code{3}.
#' @param n_basis A positive integer specifying the number of B-spline basis
#'   functions. Defaults to \code{6}.
#' @param sp_weight A numeric scalar in \eqn{[0, 1]} specifying the spatial
#'   weight in the combined distance matrix. Defaults to \code{0} (no spatial
#'   effect).
#' @param tol A positive numeric scalar specifying the ADMM stopping tolerance.
#'   Defaults to \code{0.1}.
#' @param max_iter A positive integer specifying the maximum number of ADMM
#'   iterations. Defaults to \code{10000}.
#'
#' @return An object of class \code{"fdsc"} with the following components:
#'   \describe{
#'     \item{membership}{A list of integer vectors, each giving the observation
#'       indices belonging to one estimated cluster.}
#'     \item{cluster}{A positive integer giving the estimated number of
#'       clusters.}
#'     \item{label}{An integer vector of length \code{n} giving the cluster
#'       label of each observation.}
#'     \item{method}{Character string \code{"FDSC-MST"}.}
#'     \item{n}{Number of observations.}
#'     \item{tau}{Regularization parameter used.}
#'     \item{lambda}{Smoothing penalty parameter used.}
#'   }
#'
#' @seealso \code{\link{fdscod_mst}}
#'
#' @export
fdsc_mst <- function(Y, grid, coord = NULL, 
                      rho = 1, a = 3,
                      lambda = 0.1, tau = 20,
                      spline_order = 3,
                      n_basis = 6, sp_weight = 0, tol = 0.1,max_iter = 10000) {
  n <- nrow(Y)
  n_time <- length(grid)
  if (ncol(Y) != n_time)
    warning("ncol(Y) and length(grid) differ, please check your input")
  y <- matrix(c(t(Y)))
  n_knots <- n_basis - spline_order
  breaks <- seq(min(grid), max(grid), length = n_knots + 2)
  B_basis <- fda::bsplineS(grid, breaks, norder = spline_order)       
  B <- kronecker(diag(1, n), B_basis)
  D_matrix = matrix(0, n_basis - 2, n_basis)
  for (i in 1:(n_basis - 2)) {
    D_matrix[i, i] = 1
    D_matrix[i, i + 2] = 1
    D_matrix[i, i + 1] = -2
  }
  D2 = t(D_matrix) %*% D_matrix
  diff_mat = kronecker(diag(1, n), D2)
  distance_matrix <- matrix(0, n, n)
  for(i in 1:n) {
    for(j in 1:n) {
      if(!is.null(coord)) {
        sp_dist <- sqrt(sum((coord[i, ] - coord[j, ])^2) / 2)
      } else {
        sp_dist <- 0
      }
      func_dist <- sqrt(mean((Y[i, ] - Y[j, ])^2))
      distance_matrix[i, j] <- sp_weight * sp_dist + (1 - sp_weight) * func_dist
    }
  }
  graph_cut <- igraph::graph_from_adjacency_matrix(distance_matrix, weighted = TRUE, mode = "undirected")
  mst_cut <- igraph::mst(graph_cut)
  edges_cut <- igraph::as_edgelist(mst_cut)
  mst_mat_mst <- matrix(0, n - 1, n)
  for (i in seq_len(n - 1)) {
    mst_mat_mst[i, edges_cut[i, 1]] <-  1
    mst_mat_mst[i, edges_cut[i, 2]] <- -1
  }
  A_mst <- kronecker(mst_mat_mst, diag(n_basis))
  C_init <- crossprod(B) + 0.1 * kronecker(diag(n), D2) + 0.001 * crossprod(A_mst)
  chol_C_init <- chol(C_init)
  gamma_initial <- backsolve(chol_C_init, forwardsolve(chol_C_init, crossprod(B, y)))
  vl_initial <- mst_mat_mst %*% matrix(gamma_initial, nrow = n, byrow = TRUE)
  zeta_initial <- matrix(0, n - 1, n_basis)
  est <- admm_mst(Y = Y, grid = grid, B = B, A = A_mst, diff_mat = diff_mat, 
                         rho = rho, a = a, mst_mat = mst_mat_mst,
                         gamma_old = gamma_initial, vl_old = vl_initial, zeta_old = zeta_initial,
                         tau = tau, lambda = lambda, tol = tol, max_iter = max_iter)
  grp <- recover_clusters(Y = Y, vl = est$vl, edges = edges_cut)
  new_fdsc(
    membership = grp$membership,
    cluster  = grp$cluster,
    n        = n,
    tau      = tau,
    lambda   = lambda,
    method   = "FDSC-MST"
  )
}

################################################################################
########## End FDSC-MST algorithm ##################################################
