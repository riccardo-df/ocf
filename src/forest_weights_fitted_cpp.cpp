#include <Rcpp.h>
using namespace Rcpp;

//' Forest In-Sample Honest Weights
//'
//' Computes forest in-sample honest weights for a \code{ocf.forest} object relative to the m-th class.
//'
//' @param leaf_IDs_train_list List of size \code{n.trees}, storing the leaf of each tree where training units fall into.
//' @param leaf_IDs_honest_list List of size \code{n.trees}, storing the leaf of each tree where honest units fall into.
//' @param leaf_size_honest_list List of size \code{n.trees}, storing the size of the leaves of each tree computed with honest units.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix forest_weights_fitted_cpp(List leaf_IDs_train_list, List leaf_IDs_honest_list, List leaf_size_honest_list) { // Taken from https://github.com/okasag/orf/blob/master/orf/src/get_weights_rcpp.cpp
  // Declaring variables.
  int nlist = leaf_IDs_train_list.size(); // Number of trees.
  
  NumericVector f_rows = as<NumericVector>(leaf_IDs_train_list[1]);
  NumericVector f_cols = as<NumericVector>(leaf_IDs_honest_list[1]);
  
  int nf_rows = f_rows.size(); // Number of rows of training data.
  int nf_cols = f_cols.size(); // Number of rows of honest data.
  
  NumericMatrix forest_out_train(nf_rows, nf_cols); // Matrix where rows are rows of training and columns are rows of honest.
  NumericMatrix forest_out_honest(nf_cols, nf_cols); // Matrix where rows are rows of honest and columns are rows of honest.
  NumericMatrix forest_out_all(nf_cols+nf_rows, nf_cols); // Matrix where rows are rows of honest and train and columns are rows of honest.
  
  // Looping over trees.
  for(int l = 0; l < nlist; ++l) {
    NumericVector leaf_IDs_train = as<NumericVector>(leaf_IDs_train_list[l]); // Leaves ids of l-th tree with training units.
    NumericVector leaf_IDs_honest = as<NumericVector>(leaf_IDs_honest_list[l]); // Leaves ids of l-th tree with honest units.
    NumericVector leaf_size_honest = as<NumericVector>(leaf_size_honest_list[l]); // Size of leaves of l-th tree computed with honest units.
    
    int n_leaf_IDs_train = leaf_IDs_train.size(); // How many different leaves in this tree using training units.
    int n_leaf_IDs_honest = leaf_IDs_honest.size(); // How many different leaves in this tree using honest units.
    
    NumericMatrix tree_out_train(n_leaf_IDs_train, n_leaf_IDs_honest); // Matrix where rows are leaves in training and columns are leaves in honest.
    NumericMatrix tree_out_honest(n_leaf_IDs_honest, n_leaf_IDs_honest); // Matrix where rows are leaves in honest and columns are leaves in honest.
    
    // Loop to go element by element and check the equality of leaves in each tree for training.
    // This is to predict in training sample.
    for(int i = 0; i < n_leaf_IDs_train; ++i) { // For each unit in training sample.
      for(int j = 0; j < n_leaf_IDs_honest; ++j) {  // For each unit in honest sample.
        tree_out_train(i,j) = leaf_IDs_train[i] == leaf_IDs_honest[j]; // TRUE if this i-th training unit falls in same leaf as j-th honest unit.
        tree_out_train(i,j) = tree_out_train(i,j) / leaf_size_honest[j]; // Normalize by leaf size using honest.
      }
    }
    
    // Loop to add each tree weight to overall forest weight.
    for(int i = 0; i < n_leaf_IDs_train; ++i) {
      for(int j = 0; j < n_leaf_IDs_honest; ++j) {
        forest_out_train(i,j) = forest_out_train(i,j) + tree_out_train(i,j); // Add 0 if do not share leaf.
      }
    }
    
    // Now do the same weight computation for honest sample.
    // This is to predict in honest sample.
    for(int i = 0; i < n_leaf_IDs_honest; ++i) {
      for(int j = 0; j < n_leaf_IDs_honest; ++j) {
        tree_out_honest(i,j) = leaf_IDs_honest[i] == leaf_IDs_honest[j];
        tree_out_honest(i,j) = tree_out_honest(i,j) / leaf_size_honest[j]; 
      }
    }
    
    for(int i = 0; i < n_leaf_IDs_honest; ++i) {
      for(int j = 0; j < n_leaf_IDs_honest; ++j) {
        forest_out_honest(i,j) = forest_out_honest(i,j) + tree_out_honest(i,j);
      }
    }
    
    // check for user interruptions
    Rcpp::checkUserInterrupt();
  }
  
  // Loop to divide each element of a matrix by number of trees to get mean weights.
  for(int i = 0; i < nf_rows; ++i) {
    for(int j = 0; j < nf_cols; ++j) {
      forest_out_train(i,j) = forest_out_train(i,j) / nlist ;
    }
  }
  
  for(int i = 0; i < nf_cols; ++i) {
    for(int j = 0; j < nf_cols; ++j) {
      forest_out_honest(i,j) = forest_out_honest(i,j) / nlist ;
    }
  }
  
  // Binding matrices. We store first weights for honest sample and then for training sample.
  for (int row_idleaf_IDs_train = 0; row_idleaf_IDs_train < nf_cols + nf_rows; ++row_idleaf_IDs_train) {
    if (row_idleaf_IDs_train < nf_cols) {
      forest_out_all(row_idleaf_IDs_train,_) = forest_out_honest(row_idleaf_IDs_train,_);
    } else {
      forest_out_all(row_idleaf_IDs_train,_) = forest_out_train(row_idleaf_IDs_train-nf_cols,_);
    }
  }
  
  return forest_out_all;
}
