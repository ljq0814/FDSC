#' Compute Test Statistics for Functional Outlier Detection
#'
#' Computes test statistics for potential outliers using functional principal
#' component analysis (FPCA). The statistic for each potential outlier is the
#' weighted sum of squared projections onto the leading functional principal
#' components.
#'
#' @param Y A numeric matrix of dimension \code{n x T}, representing
#'   functional data for \code{n} observations evaluated at \code{T} time
#'   points.
#' @param grid A numeric vector of time points of length \code{T}.
#' @param out_pot An integer vector of indices of potential outliers, obtained
#'   from the preliminary outlier screening step.
#' @param group_id_hat A list of integer vectors giving the initial cluster
#'   memberships, obtained from the initial clustering step.
#'
#' @return A list with two components:
#'   \describe{
#'     \item{test_stat}{A numeric vector of test statistics for each potential
#'       outlier.}
#'     \item{D}{A positive integer giving the number of retained functional
#'       principal components.}
#'   }
#'
#' @keywords internal
fun_statistic <- function(Y, grid, out_pot, group_id_hat) 
{
  n_time <- length(grid)
  n = length(out_pot)
  test_stat <- rep(0,n)
  lt=list()
  ly=list()
  for (i in 1:length(group_id_hat))
  {
    lt[[i]]=grid
    ly[[i]] <- Y[group_id_hat[[i]],]
  }
  R <- fdapace::FPCA(ly,lt,list(nRegGrid=n_time,FVEthreshold=0.999))   
  mu <- R$mu
  lambda <- R$lambda
  phi <- R$phi 
  D <- dim(phi)[2]
  ksi=matrix(0,n,D)
  for (i in 1:n) 
  {   
    for (k in 1:D) {   
      ksi[i,k]<- t(matrix((Y[out_pot[i],]-mu)))%*%matrix(phi[,k])}
    test_stat[i] <- sum(ksi[i,]^2/lambda)
  }
  
  list(test_stat=test_stat,D=D)
}
##############################End of compute test statistics#####################################


#########################Multiple Hypothesis Testing for Outlier Detection########################
##############################################################################
#' Multiple Hypothesis Testing for Functional Outlier Detection
#'
#' Performs FDR-controlled multiple hypothesis testing to identify functional
#' outliers among a set of potential outliers. Test statistics are computed
#' via \code{fun_statistic()} and thresholded using an iterative FDR procedure.
#'
#' @param Y A numeric matrix of dimension \code{n x T}, representing
#'   functional data for \code{n} observations evaluated at \code{T} time
#'   points.
#' @param grid A numeric vector of time points of length \code{T}.
#' @param out_pot An integer vector of indices of potential outliers, obtained
#'   from the preliminary outlier screening step.
#' @param group_id_hat A list of integer vectors giving the initial cluster
#'   memberships, obtained from the initial clustering step.
#' @param alpha A numeric scalar in \eqn{(0, 1)} specifying the significance
#'   level for FDR-controlled outlier detection. Defaults to \code{0.05}.
#'
#' @return An integer vector of indices (within \code{out_pot}) of detected
#'   outliers.
#'
#' @keywords internal
fun_test <- function(Y, grid, out_pot, group_id_hat, alpha)
{  
  n_time <- length(grid)
  n = length(out_pot)
  test_stat=rep(0,n)
  Q_T=rep(0,n)  
  
  AA <- fun_statistic(Y = Y, grid = grid, out_pot = out_pot, group_id_hat = group_id_hat) 
  test_stat=AA$test_stat
  D=AA$D
  #test_stat=test_stat
  th_median=median(test_stat)/qchisq(0.5,D)  
  
  test_stat=test_stat/th_median
  for (i in 1:n)
  {Q_T[i]=qnorm((pchisq(test_stat[i],D)+1)/2)}
  
  ############  find the threshold t0 by iteration 
  len.tt=100
  U_thr=sqrt(2*log(n))
  t_it=seq(0,U_thr,length=len.tt)
  t_fdr=rep(0,len.tt)
  for (it in 1:len.tt)
  {
    t_fdr[it]=2*n*(1-pnorm(t_it[it]))/max(length(which(abs(Q_T)>=t_it[it])),1)
  }
  
  if (length(which(t_fdr<=alpha))>0){
    t0=t_it[min(which((t_fdr-alpha)<=0))]}
  else
  {t0=sqrt(2*log(n))}
  ###############  perform outlier detection
  out_set=which(Q_T>t0)
  return(out_set)
}
##############################End of compute test statistics#####################################fd
#############################################End of the test########################################


#' Functional Data Subgroup Clustering and Outlier Detection via MST (FDSCOD-MST)
#'
#' Performs functional data subgroup clustering with simultaneous outlier
#' detection based on minimum spanning tree (MST). The algorithm first applies
#' preliminary outlier screening, then runs \code{fdsc_mst()} on the cleaned
#' data to obtain initial cluster assignments, and finally performs
#' FDR-controlled hypothesis testing to identify outliers among the potential
#' outlying set.
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
#' @param a A positive numeric scalar specifying the concavity parameter of
#'   the MCP penalty. Defaults to \code{3}.
#' @param lambda A positive numeric scalar specifying the smoothing penalty
#'   parameter for the B-spline roughness penalty. Defaults to \code{0.1}.
#' @param tau A positive numeric scalar specifying the regularization parameter
#'   of the MCP penalty, controlling the strength of edge difference shrinkage
#'   and thus the number of clusters. Defaults to \code{20}.
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
#' @param alpha A numeric scalar in \eqn{(0, 1)} specifying the significance
#'   level for FDR-controlled outlier detection. Defaults to \code{0.05}.
#'
#' @return An object of class \code{"fdsc"} with the following components:
#'   \describe{
#'     \item{membership}{A list of integer vectors, each giving the observation
#'       indices belonging to one estimated cluster.}
#'     \item{cluster}{A positive integer giving the estimated number of
#'       clusters.}
#'     \item{label}{An integer vector of length \code{n} giving the cluster
#'       label of each observation.}
#'     \item{method}{Character string \code{"FDSCOD-MST"}.}
#'     \item{n}{Number of observations.}
#'     \item{tau}{Regularization parameter used.}
#'     \item{lambda}{Smoothing penalty parameter used.}
#'     \item{outliers}{An integer vector of detected outlier indices, or an
#'       empty integer vector if no outliers are detected.}
#'   }
#'
#' @seealso \code{\link{fdsc_mst}}
#'
#' @export
fdscod_mst<- function(Y, grid, coord = NULL, rho = 1, a = 3,
                      lambda = 0.1, tau = 20, spline_order = 3, 
                      n_basis = 6, sp_weight = 0, tol = 0.1, max_iter = 10000, alpha = 0.05) {
  n <- nrow(Y)
  n_time <- length(grid)
  n.knots <- n_basis - spline_order
  breaks <- seq(min(grid), max(grid), length = n.knots + 2)
  
  # Step 1: Preliminary Outlier Screening
  dif <- matrix(0, n, n)
  h_len <- numeric(n)
  for (i in 1:n) {
    for (j in 1:n) {
      dif[i, j] <- mean((Y[i, ] - Y[j, ])^2)
    }
    h_len[i] <- length(which(dif[i, ] < 4 * n^(-0.01)))
  }
  h.sort <- sort.int(h_len, index.return = TRUE)
  out_pot <- h.sort$ix[1:floor(0.2 * n)]
  clean_pot <- setdiff(1:n, out_pot)
  fdy_cut <- Y[clean_pot, ]
  dx_cut <- coord[clean_pot,1]
  dy_cut <- coord[clean_pot,2]
  n_cut=dim(fdy_cut)[1]
  grp <- fdsc_mst(Y = fdy_cut, coord = coord[clean_pot,],
                  grid = grid,
                  lambda = lambda,
                  rho = rho,
                  a = a,
                  tau = tau,
                  spline_order = spline_order,
                  n_basis = n_basis,
                  sp_weight = sp_weight,
                  tol = tol,
                  max_iter = max_iter)
  
  group_id_mst <- grp$membership
  K_hat_mst <- grp$cluster
  
  group_id_mst_final=list()
  for (k in 1:K_hat_mst)
  {
    group_id_mst_final[[k]] <- sort(clean_pot[group_id_mst[[k]]])
  }
  
  out_index=list()
  for (k in 1:K_hat_mst)
  {
    out_index[[k]] <- fun_test(Y = Y, grid = grid, out_pot = out_pot,
                               group_id_hat = group_id_mst_final[[k]], alpha = alpha)
  }
  
  out_final=out_index[[1]]
  for (k in 1:K_hat_mst)
  {
    out_final <- intersect(out_final,out_index[[k]])
  }
  out_set_final <- out_pot[out_final]
  
  set_m <- setdiff(out_pot,out_set_final)
  
  for (i in 1:length(set_m)){
    mses <- rep(0,K_hat_mst)
    for (kk in 1:K_hat_mst){  
      mses[kk] <- mean((apply(matrix(Y[group_id_mst_final[[kk]],],ncol=n_time),2,mean)-Y[set_m[i],])^2)
    }
    k <- which.min(mses)
    group_id_mst_final[[k]] <- union(group_id_mst_final[[k]],set_m[i])
  }

  return(new_fdsc(
    membership = group_id_mst_final,
    cluster    = K_hat_mst,
    n          = n,
    tau        = tau,
    lambda     = lambda,
    method     = "FDSCOD-MST",
    outlier    = out_set_final
  ))
}

