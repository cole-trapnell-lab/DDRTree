#' Compute the PCA projection
#'
#' @param C data matrix used for PCA projection
#' @param L number for the top principal components
#' @import irlba irlba
#' @export
pca_projection_R <- function(C, L) {
    if (L >= min(dim(C))){
        eigen_res <- eigen(C)

        U <- eigen_res$vector
        V <- eigen_res$value
        eig_sort <- sort(V, decreasing = T, index.return = T)
        eig_idx <- eig_sort$ix

        W <- U[, eig_idx[1:L]]
        return (W)
    } else{
        eigen_res <- irlba::irlba(C, nv = L)
        U <- eigen_res$u
        V <- eigen_res$v
        return (V)
    }
}

#' Get the top L eigenvalues
#' @param C data matrix used for eigendecomposition
#' @param L number for the top eigenvalues
#' @import irlba irlba
#' @export
get_major_eigenvalue <- function(C, L) {
    if (L >= min(dim(C))){
        return (base::norm(C, '2')^2);
    }else{
        #message("using irlba")
        eigen_res <- irlba(C, nv = L)
        return (max(abs(eigen_res$v)))
    }
    #     eig_sort <- sort(V, decreasing = T, index.return = T)
    #     eig_idx <- eig_sort$ix
    #
    #     W <- U[, eig_idx[1:L]]
}

#' calculate the square distance between a, b
#' @param a a matrix with \eqn{D \times N} dimension
#' @param b a matrix with \eqn{D \times N} dimension
#' @return a numeric value for the different between a and b
#' @export
sqdist_R <- function(a, b) {
    aa <- colSums(a^2)
    bb <- colSums(b^2)
    ab <- t(a) %*% b

    aa_repmat <- matrix(rep(aa, times = ncol(b)), ncol = ncol(b), byrow = F)
    bb_repmat <- matrix(rep(bb, times = ncol(a)), nrow = ncol(a), byrow = T)
    dist <- abs(aa_repmat + bb_repmat - 2 * ab)
}

#' Perform DDRTree construction
#' @param \eqn{\mathbf{X}} a matrix with \eqn{\mathbf{D \times N}} dimension which is needed to perform DDRTree construction
#' @param dimensions reduced dimension
#' @param maxIter maximum iterations
#' @param sigma bandwidth parameter
#' @param lambda regularization parameter for inverse graph embedding
#' @param ncenter number of nodes allowed in the regularization graph
#' @param param.gamma regularization parameter for k-means (the prefix of 'param' is used to avoid name collision with gamma)
#' @param tol relative objective difference
#' @param verbose emit extensive debug output
#' @return a list with W, Z, stree, Y, history
#' @export
DDRTree <- function(X,
                        dimensions = 2,
                        maxIter = 20,
                        sigma = 1e-3,
                        lambda = NULL,
                        ncenter = NULL,
                        param.gamma = 10,
                        tol = 1e-3,
                        verbose = F) {

    D <- nrow(X)
    N <- ncol(X)

    #initialization
    W <- pca_projection_R(X %*% t(X), dimensions)
    Z <- t(W) %*% X

    if(is.null(ncenter)) {
        K <- N
        Y <- Z[, 1:K]
    }
    else {
        K <- ncenter
        kmean_res <- kmeans(t(Z), K)
        Y <- kmean_res$centers
        Y <- t(Y)
    }

    if (is.null(lambda)){
        lambda = 5 * ncol(X)
    }
    ddrtree_res <- DDRTree_reduce_dim(X, Z, Y, W, dimensions, maxIter, K,  sigma,  lambda,  param.gamma, tol, verbose)

    return(list(W = ddrtree_res$W, Z = ddrtree_res$Z, stree = ddrtree_res$stree, Y = ddrtree_res$Y, history = NULL))
}
