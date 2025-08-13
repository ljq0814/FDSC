fun_PHI <- function(s) {   
  phi_1=-cos(pi*s/10)/sqrt(5)
  phi_2=sin(pi*s/10)/sqrt(5)
  PHI=cbind(phi_1,phi_2)
  list(PHI=PHI)
}

fun_mu=function(s) {
  mu=s+sin(s)
  list(mu=mu)
}

fun_mu_cluster1 <- function(s) {
  h1 <- ifelse(s >= 1 & s <= 13, 6 - abs(s - 7), 0)
  mu <- h1 
  list(mu = mu)
}

fun_mu_cluster2 <- function(s) {
  h2 <- ifelse(s >= 9 & s <= 21, 6 - abs(s - 15), 0)
  mu <- h2  
  list(mu = mu)
}

fun_mu_cluster3 <- function(s) {
  h3 <- -(ifelse(s >= 1 & s <= 13, 6 - abs(s - 7), 0))
  list(mu = h3)
}

fun_mu_cluster4 <- function(s) {
  h4 <- -(ifelse(s >= 9 & s <= 21, 6 - abs(s - 15), 0))
  list(mu = h4)
}

#################################################################
fun_data <- function(n, timerange, mu, sig_x) {
  PHI=fun_PHI(timerange)$PHI
  tau1=1
  tau2=0.5
  r1=rnorm(n,0, tau1)
  r2=rnorm(n,0, tau2)
  r=cbind(r1,r2)
  U=r%*%t(PHI)
  E_X=mvrnorm(n,rep(0,length(timerange)),diag(sig_x^2,length(timerange)))
  E_X=matrix(c(E_X),nrow=n,ncol= length(timerange))
  X=mu+U+E_X
  return(X)
}

fun_W <- function(group_id, n, K) {
  W=matrix(0,n,K)
  for (k in 1:K) {
    W[group_id[[k]],k]=1
  }
  list(W=W)
}


################################Simulation Data Generation##############################################
#################################################################
# Title: datageneration
# Description: Generates synthetic functional data for clustering analyses under various settings in our simulation.
# Supports K=2 or K=4 clusters with two spatial distribution cases (case1 and case2).
# Optionally introduces outliers (is_outlier=TRUE) for K=2 using two contamination patterns parameterized by theta.
## Input:
#' @param n Number of functional observations to simulate.
#' @param timerange numeric vector of time points
#' @param sig_x Standard deviation of the Gaussian noise added to each observation.
#' @param K Number of clusters to simulate. Must be either 2 or 4.
#' @param spatial_dist_type Specifies the spatial partitioning scheme. Must be either "case1" or "case2".
#' @param is_outlier If TRUE, outlier observations will be introduced (only supported when K = 2).
#' @param outliercase Specifies the outlier generation mechanism. 
#'   - "case1": generates sinusoidal outliers.
#'   - "case2": generates inverted weighted mixtures of cluster means.
#'   Only used when is_outlier = TRUE.
#' @param theta Outlier contamination intensity parameter. Default is 0.
## Output:
#' @return fdy a numeric matrix of dimension n * T, representing the functional data for n observations evaluated at T time points
#' @return dx numeric vector representing the spatial x-coordinate of each observation
#' @return dy numeric vector representing the spatial y-coordinate of each observation
#' @return group_id list of true cluster memberships
#' @return mu_true n x T matrix of true mean functions
##########################################################################
################################Simulation Data Generation##############################################
datageneration <- function(n, timerange, sig_x,
                           K,
                           spatial_dist_type = "case1",
                           is_outlier = FALSE,
                           outliercase = NULL,
                           theta = 0) {
  T_len <- length(timerange)
  
  # generate spatial coordinates
  dx <- runif(n, 0, 1)
  dy <- runif(n, 0, 1)
  group_id <- list()
  
  # handle outliers
  if (is_outlier) {
    if (K != 2) stop("Outlier generation only supports K=2.")
    if (is.null(outliercase)) stop("Specify outliercase='case1' or 'case2'.")
    
    if (tolower(spatial_dist_type) == "uniform") {
      group_id[[1]] <- which(dy < dx)
      group_id[[2]] <- setdiff(seq_len(n), group_id[[1]])
    } else if (tolower(spatial_dist_type) == "nonuniform") {
      group_id[[1]] <- which(dy < dx^2)
      group_id[[2]] <- setdiff(seq_len(n), group_id[[1]])
    } else {
      stop("Invalid spatial_dist_type: choose 'uniform' or 'nonuniform'.")
    }
    
    out_true <- seq(n-9, n)
    group_id[[1]] <- setdiff(group_id[[1]], out_true)
    group_id[[2]] <- setdiff(group_id[[2]], out_true)
    
    mu_true <- matrix(0, n, T_len)
    for (i in group_id[[1]]) mu_true[i, ] <- fun_mu_cluster1(timerange)$mu
    for (i in group_id[[2]]) mu_true[i, ] <- fun_mu_cluster2(timerange)$mu
    
    if (tolower(outliercase) == "case2") {
      for (i in out_true) mu_true[i, ] <- theta * sin(2 * timerange)
    } else if (tolower(outliercase) == "case1") {
      for (i in out_true) {
        s <- rbinom(1, 1, 0.5)
        mu1 <- fun_mu_cluster1(timerange)$mu
        mu2 <- fun_mu_cluster2(timerange)$mu
        mu_true[i, ] <- -theta * (s * mu1 + (1 - s) * mu2)
      }
    } else {
      stop("Invalid outliercase: choose 'case1' or 'case2'.")
    }
    
  } else {
    if (K == 2) {
      if (tolower(spatial_dist_type) == "uniform") {
        group_id[[1]] <- which(dy < dx)
        group_id[[2]] <- setdiff(seq_len(n), group_id[[1]])
      } else if (tolower(spatial_dist_type) == "nonuniform") {
        group_id[[1]] <- which(dy < dx^2)
        group_id[[2]] <- setdiff(seq_len(n), group_id[[1]])
      } else stop("Invalid spatial_dist_type for K=2.")
    } else if (K == 4) {
      if (tolower(spatial_dist_type) == "uniform") {
        group_id[[1]] <- intersect(which(dy >= dx), which(dy >= -dx + 1))
        group_id[[2]] <- intersect(which(dy >= dx), which(dy <  -dx + 1))
        group_id[[3]] <- intersect(which(dy <  dx), which(dy >= -dx + 1))
        group_id[[4]] <- setdiff(seq_len(n),
                                 c(group_id[[1]], group_id[[2]], group_id[[3]]))
      } else if (tolower(spatial_dist_type) == "nonuniform") {
        group_id[[1]] <- intersect(which(dy < 4 * (dx - 0.5)^2), which(dx <  0.5))
        group_id[[2]] <- intersect(which(dy >=4 * (dx - 0.5)^2), which(dx <  0.5))
        group_id[[3]] <- intersect(which(dy >=4 * (dx - 0.5)^2), which(dx >= 0.5))
        group_id[[4]] <- setdiff(seq_len(n),
                                 c(group_id[[1]], group_id[[2]], group_id[[3]]))
      } else stop("Invalid spatial_dist_type for K=4.")
    } else stop("Unsupported K: only 2 or 4.")
    
    mu_true <- matrix(0, n, T_len)
    for (k in seq_len(K)) {
      mu_k <- switch(k,
                     fun_mu_cluster1(timerange)$mu,
                     fun_mu_cluster2(timerange)$mu,
                     fun_mu_cluster3(timerange)$mu,
                     fun_mu_cluster4(timerange)$mu)
      mu_true[group_id[[k]], ] <- matrix(mu_k,
                                         nrow = length(group_id[[k]]),
                                         ncol = T_len,
                                         byrow = TRUE)
    }
  }
  
  fdy <- fun_data(n, timerange, mu_true, sig_x)
  
  out <- list(
    fdy = fdy,
    dx = dx,
    dy = dy,
    group_id = group_id,
    mu_true = mu_true
  )
  if (is_outlier) out$out_true <- out_true
  return(out)
}


##########################################################################################
####################################### end genarate data #################################
##########################################################################################

