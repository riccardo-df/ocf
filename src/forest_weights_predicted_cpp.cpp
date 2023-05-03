#include <Rcpp.h>
using namespace Rcpp;

//' Forest Out-of-Sample Honest Weights
//'
//' Computes forest out-of-sample honest weights for a \code{ocf.forest} object relative to the m-th class.
//'
//' @param leaf_IDs_test_list List of size \code{n.trees}, storing the leaf of each tree where training units fall into.
//' @param leaf_IDs_honest_list List of size \code{n.trees}, storing the leaf of each tree where honest units fall into.
//' @param leaf_size_honest_list List of size \code{n.trees}, storing the size of the leaves of each tree computed with honest units.
//' @param w 1 if marginal effects are being computed, 0 otherwise for normal prediction.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector forest_weights_predicted_cpp(List leaf_IDs_test_list, List leaf_IDs_honest_list, List leaf_size_honest_list, int w) {
  // Declaring variables.
  int nlist = leaf_IDs_test_list.size();
  
  NumericVector f_rows = as<NumericVector>(leaf_IDs_test_list[1]);
  NumericVector f_cols = as<NumericVector>(leaf_IDs_honest_list[1]);
  
  int nf_rows = f_rows.size(); 
  int nf_cols = f_cols.size(); 
  
  NumericMatrix forest_out(nf_rows, nf_cols);
  
  // Looping over trees.
  for(int l = 0; l < nlist; ++l) {
    NumericVector leaf_IDs_pred = as<NumericVector>(leaf_IDs_test_list[l]);
    NumericVector leaf_IDs_honest = as<NumericVector>(leaf_IDs_honest_list[l]);
    NumericVector leaf_size = as<NumericVector>(leaf_size_honest_list[l]);
    
    int n_leaf_IDs_pred = leaf_IDs_pred.size(); 
    int n_leaf_IDs = leaf_IDs_honest.size(); 
    
    NumericMatrix tree_out(n_leaf_IDs_pred, n_leaf_IDs);
    
    for(int i = 0; i < n_leaf_IDs_pred; ++i) {
      for(int j = 0; j < n_leaf_IDs; ++j) {
        tree_out(i,j) = leaf_IDs_pred[i] == leaf_IDs_honest[j];
        tree_out(i,j) = tree_out(i,j) / leaf_size[j];
      }
    }
    
    for(int i = 0; i < n_leaf_IDs_pred; ++i) {
      for(int j = 0; j < n_leaf_IDs; ++j) {
        forest_out(i,j) = forest_out(i,j) + tree_out(i,j);
      }
    }
    
    Rcpp::checkUserInterrupt();
  }
  
  for(int i = 0; i < nf_rows; ++i) {
    for(int j = 0; j < nf_cols; ++j) {
      forest_out(i,j) = forest_out(i,j) / nlist ;
    }
  }
  
  if (w == 1) { // If we are computing marginal effects.
    NumericMatrix forest_out_mean(1, nf_cols); 
    
    for(int i = 0; i < nf_cols; ++i) {
      forest_out_mean(0, i) = mean(forest_out(_, i));
    }
    
    return forest_out_mean;
  } else { // If predictions.
    return forest_out;
  }
}
