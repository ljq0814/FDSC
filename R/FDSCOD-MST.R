################to compute test statistics for identifying functional outliers##################
##############################################################################
#' @title fun_statistic
#' Input:
#' @param timerange numeric vector of time points
#' @param fdy a numeric matrix of dimension n * T, representing the functional data for n observations evaluated at T time points
#' @param out_pot the potential outlying set obtained from the Preliminary Outlier Screening step in Algorithm 2
#' @param group_id_hat initial cluster memberships obtained from the Initial Clustering step in Algorithm 2
#' Output:
#' @return T_stat: Test statistics for potential outliers (see equation (7) in the paper)
#' @return D: number of retained functional principal components
##############################################################################
fun_statistic=function(timerange,fdy,out_pot,group_id_hat) 
{
  T <- length(timerange)
  n = length(out_pot)
  T_stat <- rep(0,n)
  lt=list()
  ly=list()
  for (i in 1:length(group_id_hat))
  {
    lt[[i]]=timerange
    ly[[i]] <- fdy[group_id_hat[[i]],]
  }
  R <- FPCA(ly,lt,list(nRegGrid=T,FVEthreshold=0.999))   
  mu <- R$mu
  lambda <- R$lambda
  phi <- R$phi 
  D <- dim(phi)[2]
  ksi=matrix(0,n,D)
  for (i in 1:n) 
  {   
    for (k in 1:D) {   
      ksi[i,k]<- t(matrix((fdy[out_pot[i],]-mu)))%*%matrix(phi[,k])}
    T_stat[i] <- sum(ksi[i,]^2/lambda)
  }
  
  list(T_stat=T_stat,D=D)
}
##############################End of compute test statistics#####################################


#########################Multiple Hypothesis Testing for Outlier Detection########################
##############################################################################
#' Input:
#' @param timerange numeric vector of time points
#' @param fdy a numeric matrix of dimension n * T, representing the functional data for n observations evaluated at T time points
#' @param out_pot the potential outlying set obtained from the Preliminary Outlier Screening step in Algorithm 2
#' @param group_id_hat initial cluster memberships obtained from the Initial Clustering step in Algorithm 2
#' @param alpha0 default 0.05 the significance level for FDR-controlled outlier detection
#' Output:
#' @return out_set: the final detected outlier set
############################################################################
fun_test <- function(timerange,fdy,out_pot,group_id_hat, alpha0)
{  
  T <- length(timerange)
  n = length(out_pot)
  T.test=rep(0,n)
  Q_T=rep(0,n)  
  
  AA=fun_statistic(timerange,fdy,out_pot,group_id_hat) 
  T.test=AA$T_stat
  D=AA$D
  #T.test=T_stat
  th_median=median(T.test)/qchisq(0.5,D)  
  
  T.test=T.test/th_median
  for (i in 1:n)
  {Q_T[i]=qnorm((pchisq(T.test[i],D)+1)/2)}
  
  ############  find the threshold t0 by iteration 
  len.tt=100
  U_thr=sqrt(2*log(n))
  t_it=seq(0,U_thr,length=len.tt)
  t_fdr=rep(0,len.tt)
  for (it in 1:len.tt)
  {
    t_fdr[it]=2*n*(1-pnorm(t_it[it]))/max(length(which(abs(Q_T)>=t_it[it])),1)
  }
  # print(T.test)
  # print(Q_T)
  
  if (length(which(t_fdr<=alpha0))>0){
    t0=t_it[min(which((t_fdr-alpha0)<=0))]}
  else
  {t0=sqrt(2*log(n))}
  # print(t_fdr)
  #print(t0)  
  ###############  perform outlier detection
  out_set=which(Q_T>t0)
  list(out_set=out_set)
}
##############################End of compute test statistics#####################################fd
#############################################End of the test########################################


##########################################
#' @title fdscod_mst
#' @description
#' Implements the Functional Data Spatial Clustering with Outlier Detection 
#' and Minimum Spanning Tree optimization (FDSCOD-MST) algorithm.
#' This function performs clustering on functional spatial data while detecting 
#' potential outliers. It first applies preliminary outlier screening based on 
#' statistical testing, then fits B-spline basis functions to the functional data, 
#' optimizes the regularization parameter via BIC, and constructs a Minimum 
#' Spanning Tree to refine clustering assignments. The output includes the 
#' estimated cluster memberships, the number of detected groups, and the indices 
#' of identified outliers.
##############################################################################
## Input:
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
#' @param alpha0 The significance level (default 0.05) used for statistical testing 
## Output:
#' @return group_id_final list of estimated cluster memberships
#' @return K_hat_final estimated number of groups
#' @return out_set_final the indices of observations identified as outliers 

#####################################################################################
fdscod_mst<- function(fdy, dx, dy, timerange, rho = 1,aa = 3,
                      lambda = 0.1, tau = 20, n_order = 3, 
                      el = 6, w = 0, ep = 0.1, alpha0 = 0.05) {
  n <- nrow(fdy)
  T <- length(timerange)
  n.knots <- el - n_order
  breaks <- seq(min(timerange), max(timerange), length = n.knots + 2)
  
  # Step 1: Preliminary Outlier Screening
  dif <- matrix(0, n, n)
  h_len <- numeric(n)
  for (i in 1:n) {
    for (j in 1:n) {
      dif[i, j] <- mean((fdy[i, ] - fdy[j, ])^2)
    }
    h_len[i] <- length(which(dif[i, ] < 4 * n^(-0.01)))
  }
  h.sort <- sort.int(h_len, index.return = TRUE)
  out_pot <- h.sort$ix[1:floor(0.2 * n)]
  clean_pot <- setdiff(1:n, out_pot)
  fdy_cut <- fdy[clean_pot, ]
  dx_cut <- dx[clean_pot]
  dy_cut <- dy[clean_pot]
  n_cut=dim(fdy_cut)[1]
  grp <- fdsc_mst(fdy = fdy_cut, dx = dx_cut, dy = dy_cut,
                  timerange = timerange,
                  lambda = lambda,
                  rho = rho,
                  aa = aa,
                  tau = tau,
                  n_order = n_order,
                  el = el,
                  w = w,
                  ep = ep)
  
  group_id_mst <- grp$group_id_hat
  K_hat_mst <- grp$K_hat
  
  group_id_mst_final=list()
  for (k in 1:K_hat_mst)
  {
    group_id_mst_final[[k]] <- sort(clean_pot[group_id_mst[[k]]])
  }
  
  out_index=list()
  for (k in 1:K_hat_mst)
  {
    out_index[[k]] <- fun_test(timerange = timerange,fdy = fdy,out_pot = out_pot,
                            group_id_hat = group_id_mst_final[[k]], alpha0 = alpha0)$out_set
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
      mses[kk] <- mean((apply(matrix(fdy[group_id_mst_final[[kk]],],ncol=T),2,mean)-fdy[set_m[i],])^2)
    }
    k <- which.min(mses)
    group_id_mst_final[[k]] <- union(group_id_mst_final[[k]],set_m[i])
  }
  
  
  return(list(group_id_final = group_id_mst_final,
              K_hat_final = K_hat_mst,
              out_set_final = out_set_final))
}

