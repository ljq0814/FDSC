#' Create an fdsc object
#'
#' @param membership A list of integer vectors giving cluster memberships.
#' @param cluster A positive integer giving the estimated number of clusters.
#' @param n A positive integer giving the number of observations.
#' @param tau A numeric scalar giving the regularization parameter of MCP penalty used.
#' @param lambda A numeric scalar giving the smoothing penalty parameter used.
#' @param method A character string, either \code{"FDSC-MST"} or
#'   \code{"FDSCOD-MST"}. Defaults to \code{"FDSC-MST"}.
#' @param outlier An integer vector of detected outlier indices, or \code{NULL}
#'   (default). Only used when \code{method = "FDSCOD-MST"}.
#'
#' @return An object of class \code{"fdsc"}.
#'
#' @noRd
new_fdsc <- function(membership, cluster, n, tau, lambda, method = "FDSC-MST",
                        outlier = NULL){
    label <- integer(n)
    for (k in seq_len(cluster)){
        label[membership[[k]]] <- k
    }
    structure(
        list(
            membership = membership,
            cluster    = cluster,
            label   = label,
            method   = method,
            n        = n,
            tau      = tau,
            lambda   = lambda,
            outlier  = outlier   # NULL if FDSC-MST, integer vector if FDSCOD-MST
        ),
        class = "fdsc"
    )
}


#' Print method for fdsc objects
#'
#' @param x An object of class \code{"fdsc"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @export
print.fdsc <- function(x, ...){
    cat("\n\tFunctional Data Subgroup Clustering\n")
    cat("\tMethod:", x$method, "\n\n")
    cat("Number of observations:", x$n, "\n")
    cat("Number of clusters:    ", x$cluster, "\n")
    sizes <- sapply(x$membership, length)
    cat("Cluster sizes:         ", sizes, "\n")
    if (x$method == "FDSCOD-MST"){
        if (!is.null(x$outlier) && length(x$outlier) > 0){
            cat("Number of outliers:    ", length(x$outlier), "\n")
        } else {
            cat("Number of outliers:     0 (no outlier detected)\n")
        }
    }
    invisible(x)
}


#' Summary method for fdsc objects
#'
#' @param object An object of class \code{"fdsc"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @export
summary.fdsc <- function(object, ...){
    print(object)
    cat("\nCluster membership:\n")
    for (k in seq_len(object$cluster)){
        cat("  Cluster", k, ":", sort(object$membership[[k]]), "\n")
    }
    if (object$method == "FDSCOD-MST"){
        if (!is.null(object$outlier) && length(object$outlier) > 0){
            cat("\nOutliers:", sort(object$outlier), "\n")
        } else {
            cat("\nOutliers: none detected\n")
        }
    }
    invisible(object)
}