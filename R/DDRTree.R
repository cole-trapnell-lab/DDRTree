#' perform PCA projection
#' solve the problem size(C) = NxN, size(W) = NxL
#' max_W trace( W' C W ) : W' W = I
#' @param C a matrix with N x N dimension
#' @param L a dataframe used to generate new data for interpolation of time points
#' @return a matrix (W) with N x L dimension
#' @export
#'
# pca_projection_R <- function(C, L) {
#     eigen_res <- eigen(C)
#
#     U <- eigen_res$vector
#     V <- eigen_res$value
#     eig_sort <- sort(V, decreasing = T, index.return = T)
#     eig_idx <- eig_sort$ix
#
#     W <- U[, eig_idx[1:L]]
# }

pca_projection_R <- function(C, L) {
    #message("using irlba")
    eigen_res <- irlba(C, nv = L)

    U <- eigen_res$u
    V <- eigen_res$v
    #     eig_sort <- sort(V, decreasing = T, index.return = T)
    #     eig_idx <- eig_sort$ix
    #
    #     W <- U[, eig_idx[1:L]]
}

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

#' perform PCA projection
#' solve the problem size(C) = NxN, size(W) = NxL
#' max_W trace( W' C W ) : W' W = I
#' @param a a matrix with D x N dimension
#' @param b a matrix with D x N dimension
#' @return a numeric value for the different between a and b
#' @export
#'
sqdist_R <- function(a, b) {
    aa <- colSums(a^2)
    bb <- colSums(b^2)
    ab <- t(a) %*% b

    aa_repmat <- matrix(rep(aa, times = ncol(b)), ncol = ncol(b), byrow = F)
    bb_repmat <- matrix(rep(bb, times = ncol(a)), nrow = ncol(a), byrow = T)
    dist <- abs(aa_repmat + bb_repmat - 2 * ab)
}

#' Perform DDRTree construction
#' @param X a matrix with D x N dimension which is needed to perform DDRTree construction
#' @param params a list with the following parameters:
#' maxIter : maximum iterations
#' eps     : relative objective difference
#' dim     : reduced dimension
#' lambda  : regularization parameter for inverse graph embedding
#' sigma   : bandwidth parameter
#' param.gamma   : regularization parameter for k-means (the prefix of 'param' is used to avoid name collision with param)
#' @return a list with W, Z, stree, Y, history
#' @export
#' gamma   : regularization parameter for k-means
#'
DDRTree_R <- function(X,  dimensions = 2,
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
        message("running k-means clustering")
        K <- ncenter
        kmean_res <- kmeans(t(Z), K)
        Y <- kmean_res$centers
        Y <- t(Y)
    }
    if (is.null(lambda)){
        lambda = 5 * ncol(X)
    }

    #main loop:
    objs <- c()
    history <- list()
    for(iter in 1:maxIter) {

        # 		#Kruskal method to find optimal B (use RBGL algorithm: http://stackoverflow.com/questions/16605825/minimum-spaning-tree-with-kruskal-algorithm)
        distsqMU <- sqdist_R(Y, Y)
        # 		#convert with graph packagege to BAM class of graph an calculate mst
        # 		mstKruskalBAM <- mstree.kruskal(graphBAM(as.data.frame(distsqMU)))
        # 		#build new data frame with resut
        # 		stree <- data.frame(cbind(t(mstKruskalBAM$edgeList),
        # 		                                 t(mstKruskalBAM$weight)))

        ##########################use mst from igraph: ##########################
        g <- graph.adjacency(distsqMU, mode = 'lower', diag = T, weighted = T)
        g_mst <- mst(g)
        stree <- get.adjacency(g_mst, attr = 'weight', type = 'lower')
        stree_ori <- stree

        #convert to matrix:
        stree <- as.matrix(stree)
        stree <- stree + t(stree)
        B_tmp <- stree != 0
        B <- B_tmp
        B[B_tmp == FALSE] <- 0
        B[B_tmp == TRUE] <- 1
        L <- diag(colSums(B)) - B

        # #convert back to igraph package
        # stree <- graph.data.frame(mstKruskalDF, directed=FALSE)

        #compute R usingmean-shift update rule
        distZY <- sqdist_R(Z, Y)
        min_dist <- matrix(rep(apply(distZY, 1, min), times = K), ncol = K, byrow = F)
        tmp_distZY <- distZY - min_dist
        tmp_R <- exp(-tmp_distZY / sigma)
        #print(tmp_R)
        R <- tmp_R / matrix(rep(rowSums(tmp_R), times = K), byrow = F, ncol = K)
        #print(R)
        Gamma_mat <- matrix(rep(0, ncol(R) ^ 2), nrow = ncol(R))
        diag(Gamma_mat) <- colSums(R)

        #termination condition
        obj1 <- - sigma * sum(log(rowSums(exp(-tmp_distZY / sigma)))
                                     - min_dist[, 1] /sigma)
        objs[iter] <- (base::norm(X - W %*% Z, '2'))^2 + lambda * sum(diag(Y %*% L %*% t(Y))) + param.gamma * obj1 #sum(diag(A))

        if(verbose)
            message('iter = ', iter, ' ', objs[iter])

        history$W[iter] <- W
        history$Z[iter] <- Z
        history$Y[iter] <- Y
        history$stree[iter] <- stree
        history$R[iter] <- R

        if(iter > 1) {
            if(abs(objs[iter] - objs[iter - 1]) / abs(objs[iter - 1]) < tol) {
                break
            }

        }

        #compute low dimension projection matrix
        tmp <- t(solve((((param.gamma + 1) / param.gamma) * ((lambda / param.gamma) * L + Gamma_mat) - t(R) %*% R), t(R)))
        Q <- 1 / (param.gamma + 1) * (diag(1, N) + tmp %*% t(R))
        C <- X %*% Q
        tmp1 <- C %*% t(X)
        W <- pca_projection_R((tmp1 + t(tmp1)) / 2, dimensions)
        Z <- t(W) %*% C
        Y <- t(solve((lambda / param.gamma * L + Gamma_mat), t(Z %*% R)))
        #print (Y)
    }

    history$objs <- objs

    return(list(W = W, Z = Z, stree = stree_ori, Y = Y, history = history))
}

#' Perform DDRTree construction
#' @param X a matrix with D x N dimension which is needed to perform DDRTree construction
#' @param params a list with the following parameters:
#' maxIter : maximum iterations
#' eps     : relative objective difference
#' dim     : reduced dimension
#' lambda  : regularization parameter for inverse graph embedding
#' sigma   : bandwidth parameter
#' param.gamma   : regularization parameter for k-means (the prefix of 'param' is used to avoid name collision with param)
#' @return a list with W, Z, stree, Y, history
#' @export
#' gamma   : regularization parameter for k-means
#'
DDRTree_cpp <- function(X,
                        dimensions = 2,
                        maxIter = 20,
                        sigma = 1e-3,
                        lambda = NULL,
                        ncenter = NULL,
                        gamma = 10,
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
    ddrtree_res <- DDRTree_reduce_dim(X, Z, Y, W, dimensions, maxIter, K,  sigma,  lambda,  gamma, tol, verbose)

    return(list(W = ddrtree_res$W, Z = ddrtree_res$Z, stree = ddrtree_res$stree, Y = ddrtree_res$Y, history = NULL))
}
