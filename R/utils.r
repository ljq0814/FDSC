### This function is for evalation: to compute RI
### This function is designed to align group index for true group and estimated group, then compute RI
fun_align=function(group_id,group_id_hat,K,K_hat)
{
  group_id_final=list()
  group_id_new=list()
  
  if (K<=K_hat)
  {
    group_id_new=group_id
    k_choose=Error=rep(0,K)
    for (k in 1:K)
    {
      same_len=rep(0,K_hat)
      for (kk in 1:K_hat )
      {
        same_len[kk]=length(intersect(group_id[[k]],group_id_hat[[kk]]))
      }
      k_choose[k]=which.max(same_len)
      group_id_final[[k]]=group_id_hat[[k_choose[k]]]
    }
    if (K<K_hat)
    {
      k_cut=setdiff(c(1:K_hat),k_choose)
      for (i in 1:length(k_cut))
      {group_id_final[[K+i]]=group_id_hat[[k_cut[i]]]}
    }
  }
  
  if (K>K_hat){
    group_id_final=group_id_hat
    k_choose=rep(0,K_hat)
    for (kk in 1:K_hat )
    {same_len=rep(0,K)
    for (k in 1:K)
    {
      same_len[k]=length(intersect(group_id[[k]],group_id_hat[[kk]]))
    }
    k_choose[kk]=which.max(same_len)
    group_id_new[[kk]]=group_id[[k_choose[kk]]]
    }
    
    k_cut=setdiff(c(1:K),k_choose)
    for (i in 1:length(k_cut))
    {group_id_new[[K_hat+i]]=group_id[[k_cut[i]]]}
  }
  
  
  list(group_id_new=group_id_new,group_id_final=group_id_final)
}


## This function is for evaluation: to compute MSE
fun_error=function(group_id,group_id_hat,K,K_hat,fdy,mu_true) 
{
  ######################################################################
  n <- nrow(fdy)
  n_order=3
  el=6
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
  
  ### initail value
  tau_in=0.001
  lambda=0.1
  
  if (K<=K_hat){
    Error=rep(0,K)
    
    for (k in 1:K)
    {
      y1=matrix(c(t(fdy[group_id_hat[[k]],]))) 
      C <- crossprod(kronecker(diag(1,length(group_id_hat[[k]])),B.basis))+ lambda*kronecker(diag(1,length(group_id_hat[[k]])),D2)  
      chol_C <- chol(C)
      alpha1=backsolve(chol_C, forwardsolve(chol_C,crossprod(kronecker(diag(1,length(group_id_hat[[k]])),B.basis), y1)))
      
      alpha1.matrix=matrix(alpha1,length(group_id_hat[[k]]),el,byrow=TRUE)
      alpha1_hat=apply(alpha1.matrix,2,mean)
      
      Error[k]=mean((matrix(mu_true[group_id[[k]][1],])-B.basis%*%matrix(alpha1_hat,ncol=1))^2)
    }
  }
  
  if (K>K_hat){
    Error=rep(0,K_hat)
    
    for (k in 1:K_hat)
    {
      y1=matrix(c(t(fdy[group_id_hat[[k]],]))) 
      C <- crossprod(kronecker(diag(1,length(group_id_hat[[k]])),B.basis))+ lambda*kronecker(diag(1,length(group_id_hat[[k]])),D2)  
      chol_C <- chol(C)
      alpha1=backsolve(chol_C, forwardsolve(chol_C,crossprod(kronecker(diag(1,length(group_id_hat[[k]])),B.basis), y1)))
      
      alpha1.matrix=matrix(alpha1,length(group_id_hat[[k]]),el,byrow=TRUE)
      alpha1_hat=apply(alpha1.matrix,2,mean)
      
      Error[k]=mean((matrix(mu_true[group_id[[k]][1],])-B.basis%*%matrix(alpha1_hat,ncol=1))^2)
    }
  }
  list(Error=Error)
}
####






####################################Evaluation#######################################
#' @title  evaluation_fun
#' @description
#' Evaluates clustering performance under two scenarios: with and without outlier detection.
#' When \code{is_outlier = FALSE}, the function computes the Adjusted Rand Index (ARI) for clustering accuracy 
#' and the Mean Squared Error (MSE) for estimation accuracy of the mean functions.  
#' When \code{is_outlier = TRUE}, the function additionally computes the True Positive Rate (TPR) 
#' and False Positive Rate (FPR) to assess the accuracy of outlier identification, 
#' where TPR is the proportion of correctly detected outliers among all true outliers, 
#' and FPR is the proportion of falsely detected outliers among all non-outliers.
#' Input:
#' @param group_id list of true cluster memberships
#' @param group_id_hat list of estimated cluster memberships
#' @param K number of true clusters
#' @param K_hat estimated number of groups
#' @param fdy n x T matrix of functional data
#' @param mu_true n x T matrix of true mean functions
#' @param is_outlier logical; if TRUE, compute outlier-related metrics (TPR, FPR) in addition to RI and MSE
#' @param out_true integer vector of indices of true outliers; required if is_outlier = TRUE
#' @param out_set_final integer vector of indices of detected outliers; required if is_outlier = TRUE
#' Output:
#' @return RI adjusted Rand index
#' @return MSE mean squared error
#' @return TPR true positive rate for outlier detection (only if is_outlier = TRUE)
#' @return FPR false positive rate for outlier detection (only if is_outlier = TRUE)
################################################################################
evaluation_fun <- function(group_id, group_id_hat, K, K_hat, fdy, mu_true,
                           is_outlier = FALSE, out_true = NULL, out_set_final = NULL) {
  library(flexclust) 
  
  n <- nrow(fdy)
  
  # Align true and estimated cluster labels
  align_res <- fun_align(group_id, group_id_hat,K, K_hat)
  group_id_mst_final <- align_res$group_id_final
  group_id_new       <- align_res$group_id_new
  
  id_true <- integer(n)
  for (k in seq_len(K)) {
    id_true[group_id_new[[k]]] <- k
  }
  id_pred <- integer(n)
  for (k in seq_len(K_hat)) {
    id_pred[group_id_mst_final[[k]]] <- k
  }
  
  RI_val <- comPart(id_true, id_pred, type = "ARI")  ### ARI
  
  # Compute Mean Squared Error (MSE)
  err_res <- fun_error(group_id_new, group_id_mst_final, K,K_hat, fdy, mu_true)
  MSE_val <- mean(err_res$Error)
  
  if (!is_outlier) {
    return(list(RI = RI_val, MSE = MSE_val))
  } else {
    # TPR/FPR for outlier setting
    if (is.null(out_true) || is.null(out_set_final)) {
      stop("When is_outlier = TRUE, please supply both out_true and out_set_final.")
    }
    TPR <- length(which(out_set_final %in% out_true)) / length(out_true)
    FPR <- length(which(!(setdiff(1:n, out_true) %in% setdiff(1:n, out_set_final)))) / (n - length(out_true))
    
    return(list(RI = RI_val, MSE = MSE_val, TPR = TPR, FPR = FPR))
  }
}