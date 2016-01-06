#include "DDRTree.h"

#include <boost/functional/hash.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

//using namespace boost;
//using boost::functional;
using namespace Rcpp;
using namespace Eigen;

typedef Eigen::SparseMatrix<double> SpMat;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
SEXP pca_projection(SEXP R_C, int dimensions){
    NumericMatrix Rcpp_C(R_C);
    const int n = Rcpp_C.nrow(), p = Rcpp_C.ncol();
    Map<MatrixXd> C(Rcpp_C.begin(), n, p);

    MatrixXd W;
    pca_projection_cpp(C, dimensions, W);
    return wrap(W);
}

void pca_projection_cpp(const MatrixXd& C, int dimensions,  MatrixXd& W){
    EigenSolver<MatrixXd> es(C, true);

    MatrixXd eVecs = es.eigenvectors().real();
    VectorXd eVals = es.eigenvalues().real();

    // Sort by ascending eigenvalues:
    std::vector<std::pair<double,MatrixXd::Index> > D;
    D.reserve(eVals.size());
    for (MatrixXd::Index i=0;i<eVals.size();i++)
        D.push_back(std::make_pair<double,MatrixXd::Index>(eVals.coeff(i,0),i));
    std::sort(D.rbegin(),D.rend());
    MatrixXd sortedEigs;
    sortedEigs.resize(eVecs.rows(), dimensions);
    for (int i=0; i < eVals.size() && i < dimensions; i++)
    {
        eVals.coeffRef(i,0)=D[i].first;
        sortedEigs.col(i)=eVecs.col(D[i].second);
    }
    W = sortedEigs;
}


// [[Rcpp::export]]
SEXP sqdist(SEXP R_a, SEXP R_b){
    NumericMatrix Rcpp_a(R_a);
    const int a_n = Rcpp_a.nrow(), a_p = Rcpp_a.ncol();
    Map<MatrixXd> a(Rcpp_a.begin(), a_n, a_p);

    NumericMatrix Rcpp_b(R_b);
    const int b_n = Rcpp_b.nrow(), b_p = Rcpp_b.ncol();
    Map<MatrixXd> b(Rcpp_b.begin(), b_n, b_p);

    MatrixXd W;
    sq_dist_cpp(a, b, W);
    return wrap(W);
}

void sq_dist_cpp(const MatrixXd& a, const MatrixXd& b,  MatrixXd& W){
//     aa <- colSums(a^2)
//     bb <- colSums(b^2)
//     ab <- t(a) %*% b
//
//     aa_repmat <- matrix(rep(aa, times = ncol(b)), ncol = ncol(b), byrow = F)
//     bb_repmat <- matrix(rep(bb, times = ncol(a)), nrow = ncol(a), byrow = T)
//     dist <- abs(aa_repmat + bb_repmat - 2 * ab)

    VectorXd aa = (a.array() * a.array()).colwise().sum();
    VectorXd bb = (b.array() * b.array()).colwise().sum();
    MatrixXd ab = a.transpose() * b;

    MatrixXd aa_repmat;
    aa_repmat.resize(b.cols(), b.cols());
    for (int i=0; i < aa_repmat.cols(); i++)
    {
        aa_repmat.col(i) = aa;
    }

    MatrixXd bb_repmat;
    bb_repmat.resize(a.cols(), a.cols());
    for (int i=0; i < bb_repmat.cols(); i++)
    {
        bb_repmat.row(i) = bb;
    }

    W = aa_repmat + bb_repmat - 2 * ab;
    W = W.array().abs().matrix();
}

void DDRTree_reduce_dim_cpp(const MatrixXd& X_in,
                            const MatrixXd& Z_in,
                            const MatrixXd& Y_in,
                            const MatrixXd& W_in,
                            int dimensions,
                            int maxIter,
                            int num_clusters,
                            double sigma,
                            double lambda,
                            double gamma,
                            double eps,
                            bool verbose,
                            MatrixXd& Y_out,
                            SpMat& stree,
                            MatrixXd& Z_out){

    Y_out = Y_in;
    MatrixXd W_out = W_in;
    Z_out = Z_in;

    int N_cells = X_in.cols();

    typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
    typedef boost::adjacency_matrix<
                                  boost::undirectedS, boost::no_property,
                                  EdgeWeightProperty> Graph;
    typedef boost::graph_traits < Graph >::edge_descriptor Edge;
    typedef boost::graph_traits < Graph >::vertex_descriptor Vertex;

    if (verbose)
        Rcpp::Rcout << "setting up adjacency matrix" << std::endl;
    Graph g(Y_in.cols());
    for (std::size_t j = 0; j < Y_in.cols(); ++j) {
        for (std::size_t i = 0; i < Y_in.cols() && i <= j ; ++i) {
            Edge e; bool inserted;
            tie(e, inserted) = add_edge(i, j, g);
        }
    }

    boost::property_map<Graph, boost::edge_weight_t>::type EdgeWeightMap = get(boost::edge_weight_t(), g);
    typedef boost::graph_traits<Graph>::edge_iterator edge_iter;

    MatrixXd B = MatrixXd::Zero(Y_in.cols(), Y_in.cols());

    std::vector < Edge > old_spanning_tree;

    std::vector<double> objective_vals;

    MatrixXd distsqMU;
    MatrixXd L;
    MatrixXd distZY;
    distZY.resize(X_in.cols(), num_clusters);

    MatrixXd min_dist;
    min_dist.resize(X_in.cols(), num_clusters);

    MatrixXd tmp_distZY;
    tmp_distZY.resize(X_in.cols(), num_clusters);

    //SpMat tmp_R(X_in.cols(), num_clusters);
    MatrixXd tmp_R;
    tmp_R.resize(X_in.cols(), num_clusters);

    //SpMat R(X_in.cols(), num_clusters);
    MatrixXd R;
    R.resize(tmp_R.cols(), num_clusters);

    //SpMat Gamma(R.cols(), R.cols());
    MatrixXd Gamma = MatrixXd::Zero(R.cols(), R.cols());

    SpMat tmp(Gamma.rows(), Gamma.cols());

    MatrixXd tmp_dense;
    tmp_dense.resize(Gamma.rows(), Gamma.cols());

    MatrixXd Q;
    Q.resize(X_in.cols(), X_in.cols());

    MatrixXd C;
    C.resize(X_in.rows(), Q.cols());

    MatrixXd tmp1;
    tmp1.resize(C.rows(), X_in.rows());

    Environment stats("package:DDRTree");
    Function pca_projection_R = stats["pca_projection_R"];

    Function get_major_eigenvalue = stats["get_major_eigenvalue"];

    for (int iter = 0; iter < maxIter; ++iter){
        if (verbose)
            Rcpp::Rcout << "************************************** " << std::endl;
        if (verbose)
            Rcpp::Rcout << "Iteration: " << iter << std::endl;

        sq_dist_cpp(Y_out, Y_out, distsqMU);
        //Rcpp::Rcout << "distsqMU: " << distsqMU<< std::endl;
        std::pair<edge_iter, edge_iter> edgePair;
        if (verbose)
            Rcpp::Rcout << "updating weights in graph" << std::endl;
        for(edgePair = edges(g); edgePair.first != edgePair.second; ++edgePair.first)
        {
            //Rcpp::Rcout << "edge: " << source(*edgePair.first,g) << " " << target(*edgePair.first,g) << " : " << distsqMU(source(*edgePair.first,g), target(*edgePair.first,g)) << std::endl;
            EdgeWeightMap[*edgePair.first] = distsqMU(source(*edgePair.first,g), target(*edgePair.first,g));
        }

        std::vector < Edge > spanning_tree;

        if (verbose)
            Rcpp::Rcout << "Finding MST" << std::endl;
        kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

        if (verbose)
            Rcpp::Rcout << "Refreshing B matrix" << std::endl;
        // update the adjacency matrix. First, erase the old edges
        for (std::vector < Edge >::iterator ei = old_spanning_tree.begin();
             ei != old_spanning_tree.end(); ++ei)
        {
            B(source(*ei, g), target(*ei, g)) = 0;
            B(target(*ei, g), source(*ei, g)) = 0;
        }

        // now add the new edges
        for (std::vector < Edge >::iterator ei = spanning_tree.begin();
             ei != spanning_tree.end(); ++ei)
        {
            B(source(*ei, g), target(*ei, g)) = 1;
            B(target(*ei, g), source(*ei, g)) = 1;
        }
        //Rcpp::Rcout << "B: " << std::endl << B << std::endl;

        old_spanning_tree = spanning_tree;

        L = B.colwise().sum().asDiagonal();
        L = L - B;

        if (verbose)
            Rcpp::Rcout << "Computing distZY" << std::endl;

        sq_dist_cpp(Z_out, Y_out, distZY);
        //Rcpp::Rcout << distZY << std::endl;

        if (verbose)
            Rcpp::Rcout << "Computing min_dist" << std::endl;
        //min_dist <- matrix(rep(apply(distZY, 1, min), times = K), ncol = K, byrow = F)

        VectorXd distZY_minCoeff = distZY.rowwise().minCoeff();
        for (int i=0; i < min_dist.cols(); i++)
        {
            min_dist.col(i) = distZY_minCoeff;
        }
        //Rcpp::Rcout << min_dist << std::endl;

        //tmp_distZY <- distZY - min_dist
        tmp_distZY = distZY - min_dist;
        //Rcpp::Rcout << tmp_distZY << std::endl;

        if (verbose)
            Rcpp::Rcout << "Computing tmpR" << std::endl;
        //tmp_R <- exp(-tmp_distZY / params$sigma)
        tmp_R = tmp_distZY.array() / (-1.0 * sigma);
        //Rcpp::Rcout << tmp_R << std::endl;

        tmp_R = tmp_R.array().exp().matrix();

        if (verbose)
            Rcpp::Rcout << "Computing R" << std::endl;
        //R <- tmp_R / matrix(rep(rowSums(tmp_R), times = K), byrow = F, ncol = K)

        VectorXd tmp_R_rowsums =  tmp_R.rowwise().sum();
        for (int i=0; i < R.cols(); i++)
        {
            R.col(i) = tmp_R_rowsums;
        }
        //Rcpp::Rcout << R << std::endl;
        //Rcpp::Rcout << "&&&&&" << std::endl;
        R = (tmp_R.array() / R.array()).matrix();
        //Rcpp::Rcout << R << std::endl;

        if (verbose)
            Rcpp::Rcout << "Computing Gamma" << std::endl;
        //Gamma <- matrix(rep(0, ncol(R) ^ 2), nrow = ncol(R))
        Gamma = MatrixXd::Zero(R.cols(), R.cols());
        //diag(Gamma) <- colSums(R)
        Gamma.diagonal() = R.colwise().sum();
        //Rcpp::Rcout << Gamma << std::endl;
        if (verbose)
            Rcpp::Rcout << "Computing obj1" << std::endl;
        //termination condition
        //obj1 <- - params$sigma * sum(log(rowSums(exp(-tmp_distZY / params$sigma))) - min_dist[, 1] / params$sigma)
        VectorXd x1 = (tmp_distZY.array() / -sigma).exp().rowwise().sum().log();
        //Rcpp::Rcout << "Computing x1 " << x1.transpose() << std::endl;
        double obj1 = -sigma * (x1 - min_dist.col(0) / sigma).sum();
        //Rcpp::Rcout << obj1 << std::endl;
        if (verbose)
            Rcpp::Rcout << "Computing obj2" << std::endl;
        //obj2 <- (norm(X - W %*% Z, '2'))^2 + params$lambda * sum(diag(Y %*% L %*% t(Y))) + params$gamma * obj1 #sum(diag(A))
        //Rcpp:Rcout << X_in - W_out * Z_out << std::endl;

        // TODO: replace this with IRLBA
        //JacobiSVD<MatrixXd> svd(X_in - W_out * Z_out, ComputeThinU | ComputeThinV);
        if (verbose){
            Rcpp::Rcout << "   X : (" << X_in.rows() << " x " << X_in.cols() << ")" << std::endl;
            Rcpp::Rcout << "   W : (" << W_out.rows() << " x " << W_out.cols() << ")" << std::endl;
            Rcpp::Rcout << "   Z : (" << Z_out.rows() << " x " << Z_out.cols() << ")" << std::endl;
        }
        double major_eigen_value = as<double>(get_major_eigenvalue(X_in - W_out * Z_out,dimensions));
        //const int X_n = W_R.nrow(), X_p = W_R.ncol();
        //Map<MatrixXd> W(W_R.begin(), X_n, X_p);


        //MatrixXf::Index maxRow, maxCol;
        //double obj2 = svd.singularValues().array().abs().maxCoeff();
        double obj2 = major_eigen_value;
        //Rcpp::Rcout << "norm = " << obj2 << std::endl;
        obj2 = obj2 * obj2;

        if (verbose){
            Rcpp::Rcout << "   L : (" << L.rows() << " x " << L.cols() << ")" << std::endl;
        }

        obj2 = obj2 + lambda * (Y_out * L * Y_out.transpose()).diagonal().sum() + gamma * obj1;
        //Rcpp::Rcout << obj2 << std::endl;
        //Rcpp::Rcout << "obj2 = " << obj2 << std::endl;
        objective_vals.push_back(obj2);

        if (verbose)
            Rcpp::Rcout << "Checking termination criterion" << std::endl;
        if(iter > 1) {
            double delta_obj = std::abs(objective_vals[iter] - objective_vals[iter - 1]);
            delta_obj /=  std::abs(objective_vals[iter - 1]);
            if(delta_obj < eps) {
                break;
            }
        }

        //Rcpp::Rcout << "L" << std::endl;
        //Rcpp::Rcout << L << std::endl;
        if (verbose)
            Rcpp::Rcout << "Computing tmp" << std::endl;
        //tmp <- t(solve( ( ( (params$gamma + 1) / params$gamma) * ((params$lambda / params$gamma) * L + Gamma) - t(R) %*% R), t(R)))

        //tmp = ((Gamma + (L * (lambda / gamma))) * ((gamma + 1.0) / gamma)).sparseView();
        //Rcpp::Rcout << tmp << std::endl;
        //tmp = tmp - (R.sparseView().transpose() * R.sparseView());
        //MatrixXd tmp = (L * (lambda / gamma) + Gamma) * ((gamma + 1.0) / gamma) - R.transpose() * R;
        //Rcpp::Rcout << tmp << std::endl;
        //tmp = tmp.llt().solve(R.transpose()).transpose();

        if (verbose)
            Rcpp::Rcout << "... stage 1" << std::endl;
        tmp = ((Gamma + (L * (lambda / gamma))) * ((gamma + 1.0) / gamma)).sparseView();
        //Rcpp::Rcout << tmp << std::endl;
        if (verbose)
            Rcpp::Rcout << "... stage 2" << std::endl;
        tmp = tmp - (R.transpose().sparseView() * R.sparseView());
        //tmp = tmp_dense.sparseView();

        if (verbose){
            Rcpp::Rcout << "Pre-computing LLT analysis" << std::endl;

            Rcpp::Rcout << "tmp is (" << tmp.rows() << "x" << tmp.cols() <<"), " << tmp.nonZeros() << " non-zero values" << std::endl;

        }

        SimplicialLLT <SparseMatrix<double> > solver;
        solver.compute(tmp);
        if(solver.info()!=Success) {
            // decomposition failed
            Rcpp::Rcout << "Error!" << std::endl;
        }
        if (verbose)
            Rcpp::Rcout << "Computing LLT" << std::endl;
        tmp_dense = solver.solve(R.transpose()).transpose();
        if(solver.info()!=Success) {
            // solving failed
            Rcpp::Rcout << "Error!" << std::endl;
        }

        //tmp_dense = tmp_dense.llt().solve(R.transpose()).transpose();

        if (verbose)
            Rcpp::Rcout << "Computing Q" << std::endl;
        //Q <- 1 / (params$gamma + 1) * (diag(1, N) + tmp %*% t(R))
        Q = ((MatrixXd::Identity(X_in.cols(), X_in.cols()) + tmp_dense * R.transpose() ).array() / (gamma + 1.0));



        if (verbose)
            Rcpp::Rcout << "Computing C" << std::endl;
        // C <- X %*% Q
        C = X_in * Q;
        //Rcpp::Rcout << C << std::endl;

        //Rcpp::Rcout << "Computing tmp1" << std::endl;
        //tmp1 <- C %*% t(X)
        tmp1 =  C * X_in.transpose();
        //Rcpp::Rcout << tmp1 << std::endl;

        if (verbose){
            Rcpp::Rcout << "Computing W" << std::endl;
            //Rcpp::Rcout << (tmp1 + tmp1.transpose()) / 2 << std::endl;
        }

        //W <- pca_projection_R((tmp1 + t(tmp1)) / 2, params$dim)

        NumericMatrix W_R = pca_projection_R((tmp1 + tmp1.transpose()) / 2,dimensions);
        const int X_n = W_R.nrow(), X_p = W_R.ncol();
        Map<MatrixXd> W(W_R.begin(), X_n, X_p);

        W_out = W;
        //pca_projection_cpp((tmp1 + tmp1.transpose()) / 2, dimensions, W_out);
        //Rcpp::Rcout << W_out << std::endl;

        if (verbose)
            Rcpp::Rcout << "Computing Z" << std::endl;
        //Z <- t(W) %*% C
        Z_out = W_out.transpose() * C;
        //Rcpp::Rcout << Z_out << std::endl;

        if (verbose)
            Rcpp::Rcout << "Computing Y" << std::endl;
        //Y <- t(solve((params$lambda / params$gamma * L + Gamma), t(Z %*% R)))
        Y_out = L * (lambda / gamma) + Gamma;
        Y_out = Y_out.llt().solve((Z_out * R).transpose()).transpose();

        //Rcpp::Rcout << Y_out << std::endl;
    }

    if (verbose)
        Rcpp::Rcout << "Clearing MST sparse matrix" << std::endl;
    stree.setZero();

    if (verbose){
        Rcpp::Rcout << "Setting up MST sparse matrix with " << old_spanning_tree.size() << std::endl;
    }
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(old_spanning_tree.size());
    // Send back the weighted MST as a sparse matrix
    for (std::vector < Edge >::iterator ei = old_spanning_tree.begin();
         ei != old_spanning_tree.end(); ++ei)
    {
        //stree.insert(source(*ei, g), target(*ei, g)) = 1;//distsqMU(source(*ei, g), target(*ei, g));
        tripletList.push_back(T( source(*ei, g), target(*ei, g), distsqMU(source(*ei, g), target(*ei, g))));
        tripletList.push_back(T( target(*ei, g), source(*ei, g), distsqMU(source(*ei, g), target(*ei, g))));
    }
    stree = SpMat(N_cells, N_cells);
    stree.setFromTriplets(tripletList.begin(), tripletList.end());
}


// [[Rcpp::export]]
Rcpp::List DDRTree_reduce_dim(SEXP R_X,
                              SEXP R_Z,
                              SEXP R_Y,
                              SEXP R_W,
                              SEXP R_dimensions,
                              SEXP R_maxiter,
                              SEXP R_num_clusters,
                              SEXP R_sigma,
                              SEXP R_lambda,
                              SEXP R_gamma,
                              SEXP R_eps,
                              SEXP R_verbose){

    //Rcpp::Rcout << "Mapping verbose" << std::endl;

    bool verbose = as<bool>(R_verbose);

    if (verbose)
        Rcpp::Rcout << "Mapping X" << std::endl;

    NumericMatrix Rcpp_X(R_X);
    const int X_n = Rcpp_X.nrow(), X_p = Rcpp_X.ncol();
    Map<MatrixXd> X(Rcpp_X.begin(), X_n, X_p);

    if (verbose)
        Rcpp::Rcout << "Mapping Z" << std::endl;

    NumericMatrix Rcpp_Z(R_Z);
    const int Z_n = Rcpp_Z.nrow(), Z_p = Rcpp_Z.ncol();
    Map<MatrixXd> Z(Rcpp_Z.begin(), Z_n, Z_p);

    if (verbose)
        Rcpp::Rcout << "Mapping Y" << std::endl;

    NumericMatrix Rcpp_Y(R_Y);
    const int Y_n = Rcpp_Y.nrow(), Y_p = Rcpp_Y.ncol();
    Map<MatrixXd> Y(Rcpp_Y.begin(), Y_n, Y_p);

    if (verbose)
        Rcpp::Rcout << "Mapping W" << std::endl;

    NumericMatrix Rcpp_W(R_W);
    const int W_n = Rcpp_W.nrow(), W_p = Rcpp_W.ncol();
    Map<MatrixXd> W(Rcpp_W.begin(), W_n, W_p);

    if (verbose)
        Rcpp::Rcout << "Mapping dimensions" << std::endl;

    int dimensions = as<int>(R_dimensions);

    if (verbose)
        Rcpp::Rcout << "Mapping maxIter" << std::endl;

    int maxiter = as<int>(R_maxiter);

    if (verbose)
        Rcpp::Rcout << "Mapping num_clusters" << std::endl;

    int num_clusters = as<int>(R_num_clusters);

    if (verbose)
        Rcpp::Rcout << "Mapping sigma" << std::endl;

    double sigma = as<double>(R_sigma);

    if (verbose)
        Rcpp::Rcout << "Mapping lambda" << std::endl;

    double lambda = as<double>(R_lambda);

    if (verbose)
        Rcpp::Rcout << "Mapping gamma" << std::endl;

    double gamma = as<double>(R_gamma);

    if (verbose)
        Rcpp::Rcout << "Mapping eps" << std::endl;

    double eps = as<double>(R_eps);



    MatrixXd Y_res;
    SpMat stree_res;
    MatrixXd Z_res;

    DDRTree_reduce_dim_cpp(X, Z, Y, W, dimensions, maxiter, num_clusters, sigma, lambda, gamma, eps, verbose, Y_res, stree_res, Z_res);

    NumericMatrix X_res;
    NumericMatrix stree;

    return Rcpp::List::create(Rcpp::Named("W") = X,
                              Rcpp::Named("Z") = Z_res,
                              Rcpp::Named("stree") = wrap(stree_res),
                              Rcpp::Named("Y") = wrap(Y_res));
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

