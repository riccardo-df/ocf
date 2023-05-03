#include <RcppEigen.h>
#include <vector>
#include <sstream>
#include <memory>
#include <utility>

#include "globals.h"
#include "Forest.h"
#include "OrderedCorrelationForest.h"
#include "Data.h"
#include "DataChar.h"
#include "DataRcpp.h"
#include "DataFloat.h"
#include "DataSparse.h"
#include "utility.h"

using namespace ocf;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Rcpp::List ocfCpp(unsigned int treetype, Rcpp::NumericMatrix& input_x, Rcpp::NumericMatrix& input_y,
    std::vector<std::string> variable_names, unsigned int mtry, unsigned int num_trees, bool verbose, unsigned int seed, 
    unsigned int num_threads, bool write_forest, unsigned int importance_mode_r, unsigned int min_node_size,
    std::vector<std::vector<double>>& split_select_weights, bool use_split_select_weights,
    std::vector<std::string>& always_split_variable_names, bool use_always_split_variable_names,
    bool prediction_mode, Rcpp::List loaded_forest, Rcpp::RawMatrix snp_data,
    bool sample_with_replacement, bool probability, std::vector<std::string>& unordered_variable_names,
    bool use_unordered_variable_names, bool save_memory, unsigned int splitrule_r, std::vector<double>& case_weights,
    bool use_case_weights, std::vector<double>& class_weights, bool predict_all, bool keep_inbag,
    std::vector<double>& sample_fraction, double alpha, double minprop, bool holdout, unsigned int prediction_type_r,
    unsigned int num_random_splits, Eigen::SparseMatrix<double>& sparse_x, 
    bool use_sparse_data, bool order_snps, bool oob_error, unsigned int max_depth, 
    std::vector<std::vector<size_t>>& inbag, bool use_inbag,
    std::vector<double>& regularization_factor, bool use_regularization_factor, bool regularization_usedepth,
    std::vector<double>& alpha_balance) {
  
  Rcpp::List result;

  try {
    std::unique_ptr<Forest> forest { };
    std::unique_ptr<Data> data { };

    // Handling inputs.
    if (!use_split_select_weights) {
      split_select_weights.clear();
    }
    
    if (!use_always_split_variable_names) {
      always_split_variable_names.clear();
    }
    
    if (!use_unordered_variable_names) {
      unordered_variable_names.clear();
    }
    
    if (!use_case_weights) {
      
      case_weights.clear();
    }
    
    if (!use_inbag) {
      inbag.clear();
    }
    
    if (!use_regularization_factor) {
      regularization_factor.clear();
    }

    std::ostream* verbose_out;
    
    if (verbose) {
      verbose_out = &Rcpp::Rcout;
    } else {
      verbose_out = new std::stringstream;
    }

    size_t num_rows;
    size_t num_cols;
    
    if (use_sparse_data) {
      num_rows = sparse_x.rows();
      num_cols = sparse_x.cols();
    } else {
      num_rows = input_x.nrow();
      num_cols = input_x.ncol();
    }
    
    // Initializing data. 
    if (use_sparse_data) {
      data = std::make_unique<DataSparse>(sparse_x, input_y, variable_names, num_rows, num_cols);
    } else {
      data = std::make_unique<DataRcpp>(input_x, input_y, variable_names, num_rows, num_cols);
    }
    
    // Using ordered forest.
    forest = std::make_unique<ForestOrdered>();

    ImportanceMode importance_mode = (ImportanceMode) importance_mode_r;
    SplitRule splitrule = (SplitRule) splitrule_r;
    PredictionType prediction_type = (PredictionType) prediction_type_r;

    // Initializing forest.
    forest->initR(std::move(data), mtry, num_trees, verbose_out, seed, num_threads,
        importance_mode, min_node_size, split_select_weights, always_split_variable_names,
        prediction_mode, sample_with_replacement, unordered_variable_names, save_memory, splitrule, case_weights,
        inbag, predict_all, keep_inbag, sample_fraction, alpha, minprop, holdout, prediction_type, num_random_splits, 
        order_snps, max_depth, regularization_factor, regularization_usedepth, alpha_balance);

    // Load forest object if in prediction mode.
    if (prediction_mode) {
      std::vector<std::vector<std::vector<size_t>> > child_nodeIDs = loaded_forest["child.nodeIDs"];
      std::vector<std::vector<size_t>> split_varIDs = loaded_forest["split.varIDs"];
      std::vector<std::vector<double>> split_values = loaded_forest["split.values"];
      std::vector<bool> is_ordered = loaded_forest["is.ordered"];

      auto& temp = dynamic_cast<ForestOrdered&>(*forest);
      temp.loadForest(num_trees, child_nodeIDs, split_varIDs, split_values, is_ordered);
    }

    // Running ocf.
    forest->run(false, oob_error);

    if (use_split_select_weights && importance_mode != IMP_NONE) {
      if (verbose_out) {
        *verbose_out
            << "Warning: Split select weights used. Variable importance measures are only comparable for variables with equal weights."
            << std::endl;
      }
    }

    // Using first non-empty dimension of predictions.
    const std::vector<std::vector<std::vector<double>>>& predictions = forest->getPredictions();
    if (predictions.size() == 1) {
      if (predictions[0].size() == 1) {
        result.push_back(forest->getPredictions()[0][0], "predictions");
      } else {
        result.push_back(forest->getPredictions()[0], "predictions");
      }
    } else {
      result.push_back(forest->getPredictions(), "predictions");
    }
    
    // Returning output.
    result.push_back(forest->getNumTrees(), "num.trees");
    result.push_back(forest->getNumIndependentVariables(), "num.covariates");

    if (!prediction_mode) {
      result.push_back(forest->getMtry(), "mtry");
      result.push_back(forest->getMinNodeSize(), "min.node.size");
      if (importance_mode != IMP_NONE) {
        result.push_back(forest->getVariableImportance(), "variable.importance");
        if (importance_mode == IMP_PERM_CASEWISE) {
          result.push_back(forest->getVariableImportanceCasewise(), "variable.importance.local");
        }
      }
    }

    if (keep_inbag) {
      result.push_back(forest->getInbagCounts(), "inbag.counts");
    }

    // Save forest if needed.
    if (write_forest) {
      Rcpp::List forest_object;
      forest_object.push_back(forest->getNumTrees(), "num.trees");
      forest_object.push_back(forest->getChildNodeIDs(), "child.nodeIDs");
      forest_object.push_back(forest->getSplitVarIDs(), "split.varIDs");
      forest_object.push_back(forest->getSplitValues(), "split.values");
      forest_object.push_back(forest->getIsOrderedVariable(), "is.ordered");
      
      result.push_back(forest_object, "forest");
    }
    
    if (!verbose) {
      delete verbose_out;
    }
  } catch (std::exception& e) {
    if (strcmp(e.what(), "User interrupt.") != 0) {
      Rcpp::Rcerr << "Error: " << e.what() << " ocf will EXIT now." << std::endl;
    }
    return result;
  }

  return result;
}
