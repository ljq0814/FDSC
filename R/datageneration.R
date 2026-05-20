#' Eigenfunctions for Functional Data Generation
#'
#' Computes two orthonormal basis functions used to generate the random
#' variation component in the functional data simulation.
#'
#' @param s A numeric vector of time points.
#'
#' @return A numeric matrix of dimension \code{length(s) x 2} giving the
#'   two eigenfunctions evaluated at \code{s}.
#'
#' @noRd
fun_PHI <- function(s) {   
  phi_1=-cos(pi*s/10)/sqrt(5)
  phi_2=sin(pi*s/10)/sqrt(5)
  PHI=cbind(phi_1,phi_2)
  return(PHI)
}

#' Mean Function for Functional Data Generation
#'
#' Computes the mean function \eqn{\mu(s) = s + \sin(s)}.
#'
#' @param s A numeric vector of time points.
#'
#' @return A numeric vector of the same length as \code{s}.
#'
#' @noRd
fun_mu <- function(s) {
  mu <- s+sin(s)
  return(mu)
}

#' Mean Function for Cluster 1
#'
#' Computes the mean function for cluster 1, defined as a triangular pulse
#' centered at \eqn{s = 7} over the interval \eqn{[1, 13]}.
#'
#' @param s A numeric vector of time points.
#'
#' @return A numeric vector of the same length as \code{s}.
#'
#' @noRd
fun_mu_cluster1 <- function(s) {
  h1 <- ifelse(s >= 1 & s <= 13, 6 - abs(s - 7), 0)
  mu <- h1 
  return(mu)
}

#' Mean Function for Cluster 2
#'
#' Computes the mean function for cluster 2, defined as a triangular pulse
#' centered at \eqn{s = 15} over the interval \eqn{[9, 21]}.
#'
#' @param s A numeric vector of time points.
#'
#' @return A numeric vector of the same length as \code{s}.
#'
#' @noRd
fun_mu_cluster2 <- function(s) {
  h2 <- ifelse(s >= 9 & s <= 21, 6 - abs(s - 15), 0)
  mu <- h2  
  return(mu)
}

#' Mean Function for Cluster 2
#'
#' Computes the mean function for cluster 2, defined as a triangular pulse
#' centered at \eqn{s = 15} over the interval \eqn{[9, 21]}.
#'
#' @param s A numeric vector of time points.
#'
#' @return A numeric vector of the same length as \code{s}.
#'
#' @noRd
fun_mu_cluster3 <- function(s) {
  h3 <- -(ifelse(s >= 1 & s <= 13, 6 - abs(s - 7), 0))
  return(h3)
}

#' Mean Function for Cluster 4
#'
#' Computes the mean function for cluster 4, defined as the negation of
#' the cluster 2 mean function.
#'
#' @param s A numeric vector of time points.
#'
#' @return A numeric vector of the same length as \code{s}.
#'
#' @noRd
fun_mu_cluster4 <- function(s) {
  h4 <- -(ifelse(s >= 9 & s <= 21, 6 - abs(s - 15), 0))
  return(h4)
}

#' Generate Functional Data from a Given Mean Function
#'
#' Generates functional observations by adding random variation and Gaussian
#' noise to a given mean function. The random variation is generated from two
#' eigenfunctions via \code{fun_PHI()}.
#'
#' @param n A positive integer specifying the number of observations.
#' @param grid A numeric vector of time points.
#' @param mu An \code{n x T} numeric matrix of mean functions, where \code{T}
#'   is the length of \code{grid}.
#' @param sig_x A positive numeric scalar specifying the standard deviation
#'   of the Gaussian noise.
#'
#' @return An \code{n x T} numeric matrix of generated functional observations.
#'
#' @noRd
fun_data <- function(n, grid, mu, sig_x) {
  PHI=fun_PHI(grid)
  tau1=1
  tau2=0.5
  r1=rnorm(n,0, tau1)
  r2=rnorm(n,0, tau2)
  r=cbind(r1,r2)
  U=r%*%t(PHI)
  E_X=MASS::mvrnorm(n,rep(0,length(grid)),diag(sig_x^2,length(grid)))
  E_X=matrix(c(E_X),nrow=n,ncol= length(grid))
  X=mu+U+E_X
  return(X)
}

fun_W <- function(group_id, n, K) {
  W=matrix(0,n,K)
  for (k in 1:K) {
    W[group_id[[k]],k]=1
  }
  return(W)
}


#' Generate Synthetic Functional Data for Clustering Simulation
#'
#' Generates synthetic functional data for clustering analyses. Supports
#' \code{K = 2} or \code{K = 4} clusters with uniform or nonuniform spatial
#' distributions. Optionally introduces outliers for \code{K = 2} using two
#' contamination patterns.
#'
#' @param n A positive integer specifying the number of observations.
#' @param grid A numeric vector of time points.
#' @param sig_x A positive numeric scalar specifying the standard deviation
#'   of the Gaussian noise added to each observation.
#' @param K A positive integer specifying the number of clusters. Must be
#'   either \code{2} or \code{4}.
#' @param spatial_dist_type A character string specifying the spatial
#'   partitioning scheme, either \code{"uniform"} (default) or
#'   \code{"nonuniform"}.
#' @param is_outlier A logical value indicating whether to introduce outlier
#'   observations. Only supported when \code{K = 2}. Defaults to \code{FALSE}.
#' @param outliercase A character string specifying the outlier generation
#'   mechanism, either \code{"case1"} (inverted weighted mixtures of cluster
#'   means) or \code{"case2"} (sinusoidal outliers). Only used when
#'   \code{is_outlier = TRUE}.
#' @param theta A numeric scalar specifying the outlier contamination
#'   intensity. Defaults to \code{0}.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{Y}{An \code{n x T} numeric matrix of functional observations.}
#'     \item{coord}{An \code{n x 2} numeric matrix of spatial coordinates,
#'       where columns correspond to x and y coordinates respectively.}
#'     \item{group_id}{A list of integer vectors giving the true cluster
#'       memberships.}
#'     \item{mu_true}{An \code{n x T} numeric matrix of true mean functions.}
#'     \item{out_true}{An integer vector of true outlier indices. Only
#'       included when \code{is_outlier = TRUE}.}
#'   }
#'
#' @export
datageneration <- function(n, grid, sig_x, K,
                           spatial_dist_type = "uniform",
                           is_outlier = FALSE,
                           outliercase = NULL,
                           theta = 0) {
  T_len <- length(grid)
  
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
    for (i in group_id[[1]]) mu_true[i, ] <- fun_mu_cluster1(grid)
    for (i in group_id[[2]]) mu_true[i, ] <- fun_mu_cluster2(grid)
    
    if (tolower(outliercase) == "case2") {
      for (i in out_true) mu_true[i, ] <- theta * sin(2 * grid)
    } else if (tolower(outliercase) == "case1") {
      for (i in out_true) {
        s <- rbinom(1, 1, 0.5)
        mu1 <- fun_mu_cluster1(grid)
        mu2 <- fun_mu_cluster2(grid)
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
                     fun_mu_cluster1(grid),
                     fun_mu_cluster2(grid),
                     fun_mu_cluster3(grid),
                     fun_mu_cluster4(grid))
      mu_true[group_id[[k]], ] <- matrix(mu_k,
                                         nrow = length(group_id[[k]]),
                                         ncol = T_len,
                                         byrow = TRUE)
    }
  }
  
  fdy <- fun_data(n, grid, mu_true, sig_x)
  coord <- cbind(dx, dy)
  out <- list(
    Y = fdy,
    coord = coord,
    membership = group_id,
    mu = mu_true
  )
  if (is_outlier) out$outlier <- out_true
  return(out)
}


##########################################################################################
####################################### end genarate data #################################
##########################################################################################

