#include <Rcpp.h>
using namespace Rcpp;

//' Honest In-Sample Predictions 
//'
//' Computes honest in-sample predictions for a ocf.forest object relative to the desired class.
//'
//' @param unique_leaves_honest List of size \code{n.trees}, storing the unique leaf ids of each tree relative to the honest sample.
//' @param y_m Indicator variable, equal to 1 if the \code{y} is lower or equal than the m-th class and zero otherwise.
//' @param y_m_1 Indicator variable, equal to 1 if the \code{y} is lower or equal than the (m-1)-th class and zero otherwise.
//' @param honest_leaves Matrix of size (\code{n.samples} x \code{n.trees}). The i-th row stores the id of the leaf where the i-th honest observation falls in each tree.
//' @param train_leaves Matrix of size (\code{n.samples} x \code{n.trees}). The i-th row stores the id of the leaf where the i-th training observation falls in each tree.
//' 
//' @keywords internal
// [[Rcpp::export]]
NumericVector honest_fitted_cpp(List unique_leaves_honest, NumericVector y_m, NumericVector y_m_1, 
                                     NumericMatrix honest_leaves, NumericMatrix train_leaves) { // Taken from https://github.com/okasag/orf/blob/master/orf/src/get_honest_rcpp.cpp
  // Declaring variables.
  int leaf_ID = 0;
  int n_unique_leaves_honest = unique_leaves_honest.size();
  
  int n_rows_train = train_leaves.nrow();
  int n_rows_honest = honest_leaves.nrow();

  NumericVector obs_same_all_train(n_rows_train);
  NumericVector obs_same_all_honest(n_rows_honest);
  
  NumericVector y_m_same(n_rows_honest);
  NumericVector y_m_1_same(n_rows_honest);
  
  double y_m_mean = 0;
  double y_m_1_mean = 0;
  
  NumericMatrix train_pred(n_rows_train, n_unique_leaves_honest);
  NumericMatrix honest_pred(n_rows_honest, n_unique_leaves_honest);
  
  NumericMatrix all_pred(n_rows_honest+n_rows_train, n_unique_leaves_honest);
  NumericVector all_pred_final(n_rows_honest+n_rows_train);
  
  // Looping over trees.
  for(int i = 0; i < n_unique_leaves_honest; ++ i) {
    // Select leaves ids of i-th tree. This is the i-th element of the unique_leaves_honest list.
    NumericVector leaves = as<NumericVector>(unique_leaves_honest[i]); 
    int n_leaves = leaves.size();
    
    // Looping over the unique leaves of each tree.
    for(int leaf_idx = 0; leaf_idx < n_leaves; ++leaf_idx) { 
      // Select one leaf.
      leaf_ID = leaves[leaf_idx];
      
      // Find observations in this leaf in training sample.
      for(int row_idx = 0; row_idx < n_rows_train; ++row_idx) {
        obs_same_all_train[row_idx] = train_leaves(row_idx, i) == leaf_ID; // True if the row_idx-th training unit falls into this leaf in this tree.
      }                                                                                       
      
      // Find observations in this leaf in honest sample.
      for(int row_idx = 0; row_idx < n_rows_honest; ++row_idx) {
        obs_same_all_honest[row_idx] = honest_leaves(row_idx, i) == leaf_ID; // True if the row_idx-th honest unit falls into this leaf in this tree.
      }
      
      // If honest unit in this leaf, save associated outcome, otherwise write a zero.
      for(int row_idx = 0; row_idx < n_rows_honest; ++row_idx) {
        if (obs_same_all_honest[row_idx] == 1) { 
          y_m_same[row_idx] = y_m[row_idx];
          y_m_1_same[row_idx] = y_m_1[row_idx];
        } else {
          y_m_same[row_idx] = 0;
          y_m_1_same[row_idx] = 0;
        }
      }
      
      double y_m_count = 0;
      double y_m_sum = 0;
      
      double y_m_1_count = 0;
      double y_m_1_sum = 0;
      
      // Count and sum honest outcomes within this leaf.
      for(int row_idx = 0; row_idx < n_rows_honest; ++row_idx) {
        if (obs_same_all_honest[row_idx] == 1) { 
          y_m_count = y_m_count + obs_same_all_honest[row_idx]; // Due to above loop, here we sum 1 if the row_idx-th honest unit is in this leaf, and 0 otherwise.
          y_m_sum = y_m_sum + y_m_same[row_idx]; // Due to the above loop, here we sum outcomes if the row_idx-th honest unit is in this leaf, and 0 otherwise.
           
          y_m_1_count = y_m_1_count + obs_same_all_honest[row_idx];
          y_m_1_sum = y_m_1_sum + y_m_1_same[row_idx];
        }
      }
      
      y_m_mean = y_m_sum / y_m_count;
      y_m_1_mean = y_m_1_sum / y_m_1_count; // Difference is the honest prediction in this leaf.
      
      // If row_idx-th honest unit is in this leaf, assign this prediction.
      for(int row_idx = 0; row_idx < n_rows_honest; ++row_idx) {
        if (obs_same_all_honest[row_idx] == 1) {
          honest_pred(row_idx, i) = y_m_mean - y_m_1_mean;
        }
      }
      
      // If row_idx-th training unit is in this leaf, assign this prediction.
      for(int row_idx = 0; row_idx < n_rows_train; ++row_idx) {
        if (obs_same_all_train[row_idx] == 1) {
          train_pred(row_idx, i) = y_m_mean - y_m_1_mean;
        }
      }
    } // Closing loop across leaves of a tree.
    
    // Checking user interruptions.
    Rcpp::checkUserInterrupt();
  } // closing loop across trees.
  
  // Binding matrices. We store first predictions on honest sample and then predictions on training sample.
  for (int row_idx = 0; row_idx < n_rows_honest + n_rows_train; ++row_idx) {
    if (row_idx < n_rows_honest) {
      all_pred(row_idx,_) = honest_pred(row_idx,_);
    } else {
      all_pred(row_idx,_) = train_pred(row_idx-n_rows_honest,_);
    }
  }
  
  // Computing rowmeans for the all_pred matrix.
  for (int row_idx = 0; row_idx < n_rows_honest + n_rows_train; ++row_idx) {
    double total = 0;

    for (int col_idx = 0; col_idx < n_unique_leaves_honest; ++col_idx) {
      total += all_pred(row_idx, col_idx);
    }
    
    all_pred_final[row_idx] = total/n_unique_leaves_honest;
  }
  
  return all_pred_final;
}
