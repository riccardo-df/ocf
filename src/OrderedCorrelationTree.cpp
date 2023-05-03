#include <algorithm>
#include <iostream>
#include <iterator>

#include <ctime>

#include "utility.h"
#include "OrderedCorrelationTree.h"
#include "Data.h"

namespace ocf {

TreeOrdered::TreeOrdered(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
    std::vector<double>& split_values) :
    Tree(child_nodeIDs, split_varIDs, split_values), counter(0), sums_m(0), sums_m_1(0), prods(0) {
}

void TreeOrdered::allocateMemory() {
  // Initializing arrays if not in memory efficient mode.
  if (!memory_saving_splitting) {
    size_t max_num_splits = data->getMaxNumUniqueValues();
    
    counter.resize(max_num_splits);
    
    sums_m.resize(max_num_splits);
    sums_m_1.resize(max_num_splits);
    
    prods.resize(max_num_splits);
  }
}

double TreeOrdered::estimate(size_t nodeID) { // Predicts in each node. 
  size_t num_samples_in_node = end_pos[nodeID] - start_pos[nodeID];
  
  // Mean of class m-1.
  double sum_responses_in_node_m_1 = 0;
  for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
    size_t sampleID = sampleIDs[pos];
    sum_responses_in_node_m_1 += data->get_y(sampleID, 0); 
  }
  
  double mean_m_1 = sum_responses_in_node_m_1 / (double) num_samples_in_node;
  
  // Mean of class m.
  double sum_responses_in_node_m = 0;
  for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
    size_t sampleID = sampleIDs[pos];
    sum_responses_in_node_m += data->get_y(sampleID, 1); 
  }
  double mean_m = sum_responses_in_node_m / (double) num_samples_in_node;
  
  // Difference in means.
  return (mean_m - mean_m_1);
}

void TreeOrdered::appendToFileInternal(std::ofstream& file) { // #nocov start
  // Empty on purpose
} // #nocov end

bool TreeOrdered::splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {
  size_t num_samples_node = end_pos[nodeID] - start_pos[nodeID];

  // Stop if minimum node size or maximum depth is reached.
  if (num_samples_node <= min_node_size || (nodeID >= last_left_nodeID && max_depth > 0 && depth >= max_depth)) {
    split_values[nodeID] = estimate(nodeID);
    return true;
  }

  // Checking if node is pure and set split_value to estimate, and stopping if pure.
  bool pure = true;
  double pure_value = 0;
  
  for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
    size_t sampleID = sampleIDs[pos];
    double value = data->get_y(sampleID, 1) - data->get_y(sampleID, 0);
    
    if (pos != start_pos[nodeID] && value != pure_value) {
      pure = false;
      break;
    }
    
    pure_value = value;
  }
  
  if (pure) {
    split_values[nodeID] = pure_value;
    return true;
  }

  // Finding best split, stopping if no decrease in impurity.
  bool stop;
  stop = findBestSplit(nodeID, possible_split_varIDs);

  if (stop) {
    split_values[nodeID] = estimate(nodeID);
    return true;
  }

  return false;
}

void TreeOrdered::createEmptyNodeInternal() {
  // Empty on purpose
}

double TreeOrdered::computePredictionAccuracyInternal(std::vector<double>* prediction_error_casewise) {
  size_t num_predictions = prediction_terminal_nodeIDs.size();
  double sum_of_squares = 0;
  
  // This is to be corrected. Implement MSE for ordered choice models.
  for (size_t i = 0; i < num_predictions; ++i) {
    size_t terminal_nodeID = prediction_terminal_nodeIDs[i];
    double predicted_value = split_values[terminal_nodeID];
    double real_value = data->get_y(oob_sampleIDs[i], 1) - data->get_y(oob_sampleIDs[i], 0);
    
    if (predicted_value != real_value) {
      double diff = (predicted_value - real_value) * (predicted_value - real_value);
      if (prediction_error_casewise) {
        (*prediction_error_casewise)[i] = diff;
      }
      sum_of_squares += diff;
    }
  }
  return (1.0 - sum_of_squares / (double) num_predictions);
} 

bool TreeOrdered::findBestSplit(size_t nodeID, std::vector<size_t>& possible_split_varIDs) {
  size_t num_samples_node = end_pos[nodeID] - start_pos[nodeID];
  double best_decrease = -1;
  size_t best_varID = 0;
  double best_value = 0;
  
  // Computing sum of responses and product of indicators in the node.
  double sum_node_m = 0;
  double sum_node_m_1 = 0;
  double prod_node = 0;
  
  for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
    size_t sampleID = sampleIDs[pos];
    sum_node_m += data->get_y(sampleID, 1);
    sum_node_m_1 += data->get_y(sampleID, 0);
    prod_node += (data->get_y(sampleID, 1) * data->get_y(sampleID, 0));
  }
  
  for (auto& varID : possible_split_varIDs) { // For all possible split variables.

    //ing Find best split value. 
    if (data->isOrderedVariable(varID)) { // If ordered, consider all values as split values.
      if (memory_saving_splitting) { // Use memory saving method if option is set.
        findBestSplitValueSmallQ(nodeID, varID, sum_node_m, sum_node_m_1, prod_node, num_samples_node, best_value, best_varID, 
                                 best_decrease);
      } else { // If not in memory save mode.
        // Using the faster method.
        double q = (double) num_samples_node / (double) data->getNumUniqueDataValues(varID);
        
        if (q < Q_THRESHOLD) {
          findBestSplitValueSmallQ(nodeID, varID, sum_node_m, sum_node_m_1, prod_node, num_samples_node, best_value, best_varID, 
                                   best_decrease);
        } else {
          findBestSplitValueLargeQ(nodeID, varID, sum_node_m, sum_node_m_1, prod_node, num_samples_node, best_value, best_varID, 
                                   best_decrease);
        }
      }
    } else { // If not ordered, consider all 2-partitions.
      findBestSplitValueUnordered(nodeID, varID, sum_node_m, sum_node_m_1, prod_node, num_samples_node, best_value, best_varID, 
                                  best_decrease);
    }
  }

  // Stop if no good split found.
  if (best_decrease < 0) {
    return true;
  }

  // Saving best values.
  split_varIDs[nodeID] = best_varID;
  split_values[nodeID] = best_value;

  // Computing decrease of impurity for this node and add to variable importance if needed.
  if (importance_mode == IMP_GINI || importance_mode == IMP_GINI_CORRECTED) {
    addImpurityImportance(nodeID, best_varID, best_decrease);
  }

  // Regularization.
  saveSplitVarID(best_varID);

  return false;
}

void TreeOrdered::findBestSplitValueSmallQ(size_t nodeID, size_t varID, double sum_node_m, double sum_node_m_1,
                                              double prod_node, size_t num_samples_node, double& best_value, 
                                              size_t& best_varID, double& best_decrease) {
  // Creating possible split values.
  std::vector<double> possible_split_values;
  data->getAllValues(possible_split_values, sampleIDs, varID, start_pos[nodeID], end_pos[nodeID]);

  // Try next variable if all equal for this.
  if (possible_split_values.size() < 2) {
    return;
  }
  
  const size_t num_splits = possible_split_values.size();
  
  if (memory_saving_splitting) {
    std::vector<double> sums_right(num_splits);
    std::vector<size_t> n_right(num_splits);
    findBestSplitValueSmallQ(nodeID, varID, sum_node_m, sum_node_m_1, prod_node, num_samples_node, best_value, best_varID, 
                             best_decrease, possible_split_values, sums_m, sums_m_1, prods, counter);
  } else { 
    std::fill_n(sums_m.begin(), num_splits, 0); 
    std::fill_n(sums_m_1.begin(), num_splits, 0);
    std::fill_n(prods.begin(), num_splits, 0); 
    std::fill_n(counter.begin(), num_splits, 0);
    findBestSplitValueSmallQ(nodeID, varID, sum_node_m, sum_node_m_1, prod_node, num_samples_node, best_value, best_varID, 
                             best_decrease, possible_split_values, sums_m, sums_m_1, prods, counter);
  }
}

void TreeOrdered::findBestSplitValueSmallQ(size_t nodeID, size_t varID, double sum_node_m, double sum_node_m_1, 
                                              double prod_node, size_t num_samples_node, double& best_value, 
                                              size_t& best_varID, double& best_decrease, 
                                              std::vector<double> possible_split_values, std::vector<double>& sums_m, 
                                              std::vector<double>& sums_m_1, std::vector<double>& prods, 
                                              std::vector<size_t>& counter) {
  for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
    size_t sampleID = sampleIDs[pos]; 
    size_t idx = std::lower_bound(possible_split_values.begin(), possible_split_values.end(),
        data->get_x(sampleID, varID)) - possible_split_values.begin(); 
    
    sums_m[idx] += data->get_y(sampleID, 1);
    sums_m_1[idx] += data->get_y(sampleID, 0);
    prods[idx] += data->get_y(sampleID, 1) * data->get_y(sampleID, 0);
    ++counter[idx]; 
  }
  
  size_t n_left = 0;
  
  double sum_left_m = 0;
  double sum_left_m_1 = 0;
  double prod_left = 0;

  // Computing decrease of impurity for each split.
  for (size_t i = 0; i < possible_split_values.size() - 1; ++i) { 

    // Stop if nothing here.
    if (counter[i] == 0) { 
      continue;
    }

    n_left += counter[i]; // If we split at i, number of observations on the left.
    sum_left_m += sums_m[i]; // If we split at i, sums of indicators for m-th class on the left.
    sum_left_m_1 += sums_m_1[i];
    prod_left += prods[i]; // If we split at i, sums of products on the left.

    // Stop if right child is empty.
    size_t n_right = num_samples_node - n_left;
    
    if (n_right == 0) {
      break;
    }
    
    // Ignore this split value if alpha-regularity would be violated.
    double frac = alpha_balance[0];
    bool condition_left = n_left < frac * num_samples_node; // Conditions for alpha-regularity.
    bool condition_right = n_right < frac * num_samples_node;
    
    if (condition_left || condition_right) {
      continue;
    }

    double sum_right_m = sum_node_m - sum_left_m;
    double sum_right_m_1 = sum_node_m_1 - sum_left_m_1;
    double prod_right = prod_node - prod_left;
    
    double meanLeft_m = sum_left_m / n_left;
    double meanLeft_m_1 = sum_left_m_1 / n_left;
    double meanRight_m = sum_right_m / n_right;
    double meanRight_m_1 = sum_right_m_1 / n_right;
    
    // MSE for both classes.
    double mse_m = sum_left_m * sum_left_m / (double) n_left + sum_right_m * sum_right_m / (double) n_right;
    double mse_m_1 = sum_left_m_1 * sum_left_m_1 / (double) n_left + sum_right_m_1 * sum_right_m_1 / (double) n_right;
    
    // MCE.
    double mce = (prod_left / n_left) - (meanLeft_m * meanLeft_m_1) + (prod_right / n_right) - (meanRight_m * meanRight_m_1);
    
    // Total decrease in loss function.
    // double decrease = mse_m + mse_m_1 - 2 * mce;
    double decrease = mse_m + mse_m_1 + 2 * mce;

    // Regularization.
    regularize(decrease, varID);

    // If better than before, use this.
    if (decrease > best_decrease) {
      // Using mid-point split
      best_value = (possible_split_values[i] + possible_split_values[i + 1]) / 2;
      best_varID = varID;
      best_decrease = decrease;

      // Using smaller value if average is numerically the same as the larger value.
      if (best_value == possible_split_values[i + 1]) {
        best_value = possible_split_values[i];
      }
    }
  }
}

void TreeOrdered::findBestSplitValueLargeQ(size_t nodeID, size_t varID, double sum_node_m, double sum_node_m_1,
                                              double prod_node, size_t num_samples_node, double& best_value, 
                                              size_t& best_varID, double& best_decrease) {
  // Setting counters to 0.
  size_t num_unique = data->getNumUniqueDataValues(varID);
  std::fill_n(counter.begin(), num_unique, 0);
  std::fill_n(sums_m.begin(), num_unique, 0);
  std::fill_n(sums_m_1.begin(), num_unique, 0);
  std::fill_n(prods.begin(), num_unique, 0);

  // For each unit in the node, getting sums and products of class indicators.
  for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
    size_t sampleID = sampleIDs[pos];
    size_t index = data->getIndex(sampleID, varID);

    sums_m[index] += data->get_y(sampleID, 1); 
    sums_m_1[index] += data->get_y(sampleID, 0);
    prods[index] += data->get_y(sampleID, 1) * data->get_y(sampleID, 0);

    ++counter[index];
  }
  
  // Computing decrease of impurity for each splitting value.
  size_t n_left = 0;
  
  double sum_left_m = 0;
  double sum_left_m_1 = 0;
  double prod_left = 0;
  
  for (size_t i = 0; i < num_unique - 1; ++i) {

    // Stop if nothing here.
    if (counter[i] == 0) {
      continue;
    }

    n_left += counter[i]; // If we split at i, number of observations on the left.
    sum_left_m += sums_m[i]; // If we split at i, sums of indicators for m-th class on the left.
    sum_left_m_1 += sums_m_1[i];
    prod_left += prods[i]; // If we split at i, sums of products on the left.
    
    // Stop if right child is empty.
    size_t n_right = num_samples_node - n_left;
  
    if (n_right == 0) {
      break;
    }
    
    // Ignore this split value if alpha-regularity would be violated.
    double frac = alpha_balance[0];
    bool condition_left = n_left < frac * num_samples_node; // Conditions for alpha-regularity.
    bool condition_right = n_right < frac * num_samples_node;
    
    if (condition_left || condition_right) {
      continue;
    }

    double sum_right_m = sum_node_m - sum_left_m;
    double sum_right_m_1 = sum_node_m_1 - sum_left_m_1;
    double prod_right = prod_node - prod_left;
    
    double meanLeft_m = sum_left_m / n_left;
    double meanLeft_m_1 = sum_left_m_1 / n_left;
    double meanRight_m = sum_right_m / n_right;
    double meanRight_m_1 = sum_right_m_1 / n_right;
    
    // MSE for both classes.
    double mse_m = sum_left_m * sum_left_m / (double) n_left + sum_right_m * sum_right_m / (double) n_right;
    double mse_m_1 = sum_left_m_1 * sum_left_m_1 / (double) n_left + sum_right_m_1 * sum_right_m_1 / (double) n_right;
    
    // MCE.
    double mce = (prod_left / n_left) - (meanLeft_m * meanLeft_m_1) + (prod_right / n_right) - (meanRight_m * meanRight_m_1);
    
    // Total decrease in loss function.
    // double decrease = mse_m + mse_m_1 - 2 * mce;
    double decrease = mse_m + mse_m_1 + 2 * mce;
    
    // Regularization.
    regularize(decrease, varID);

    // If better than before, use this.
    if (decrease > best_decrease) {
      // Finding next value in this node.
      size_t j = i + 1;
      
      while (j < num_unique && counter[j] == 0) {
        ++j;
      }

      // Using mid-point split.
      best_value = (data->getUniqueDataValue(varID, i) + data->getUniqueDataValue(varID, j)) / 2;
      best_varID = varID;
      best_decrease = decrease;

      // Using smaller value if average is numerically the same as the larger value.
      if (best_value == data->getUniqueDataValue(varID, j)) {
        best_value = data->getUniqueDataValue(varID, i);
      }
    }
  }
}

void TreeOrdered::findBestSplitValueUnordered(size_t nodeID, size_t varID, double sum_node_m, double sum_node_m_1,
                                                 double prod_node, size_t num_samples_node, double& best_value, 
                                                 size_t& best_varID, double& best_decrease) {
  // Creating possible split values.
  std::vector<double> factor_levels;
  data->getAllValues(factor_levels, sampleIDs, varID, start_pos[nodeID], end_pos[nodeID]);

  // Try next variable if all equal for this.
  if (factor_levels.size() < 2) {
    return;
  }

  // Number of possible splits is 2^num_levels.
  size_t num_splits = (1ULL << factor_levels.size());

  /* Compute decrease of impurity for each possible split. Split where all left (0) or all right (1) are excluded.
     The second half of numbers is just left/right switched the first half -> Exclude second half. */
  for (size_t local_splitID = 1; local_splitID < num_splits / 2; ++local_splitID) {
    // Computing overall splitID by shifting local factorIDs to global positions.
    size_t splitID = 0;
    
    for (size_t j = 0; j < factor_levels.size(); ++j) {
      if ((local_splitID & (1ULL << j))) {
        double level = factor_levels[j];
        size_t factorID = floor(level) - 1;
        splitID = splitID | (1ULL << factorID);
      }
    }

    // Initializing.
    double sum_right_m = 0;
    double sum_right_m_1 = 0;
    double prod_right = 0;
    size_t n_right = 0;

    // Sum in right child.
    for (size_t pos = start_pos[nodeID]; pos < end_pos[nodeID]; ++pos) {
      size_t sampleID = sampleIDs[pos];
      double response_m = data->get_y(sampleID, 1);
      double response_m_1 = data->get_y(sampleID, 0);
      double prods = data->get_y(sampleID, 1) * data->get_y(sampleID, 0);
      double value = data->get_x(sampleID, varID);
      size_t factorID = floor(value) - 1;

      // If in right child, count. In right child, if bitwise splitID at position factorID is 1.
      if ((splitID & (1ULL << factorID))) {
        ++n_right;
        sum_right_m += response_m;
        sum_right_m_1 += response_m_1;
        prod_right += prods;
      }
    }
    
    size_t n_left = num_samples_node - n_right;
    
    // Ignore this split value if alpha-regularity would be violated.
    double frac = alpha_balance[0];
    bool condition_left = n_left < frac * num_samples_node; // Conditions for alpha-regularity.
    bool condition_right = n_right < frac * num_samples_node;
    
    if (condition_left || condition_right) {
      continue;
    }

    double sum_left_m = sum_node_m - sum_right_m;
    double sum_left_m_1 = sum_node_m_1 - sum_right_m_1;
    double prod_left = prod_node - prod_right;
    
    double meanLeft_m = sum_left_m / n_left;
    double meanLeft_m_1 = sum_left_m_1 / n_left;
    double meanRight_m = sum_right_m / n_right;
    double meanRight_m_1 = sum_right_m_1 / n_right;
    
    // MSE for both classes.
    double mse_m = sum_left_m * sum_left_m / (double) n_left + sum_right_m * sum_right_m / (double) n_right;
    double mse_m_1 = sum_left_m_1 * sum_left_m_1 / (double) n_left + sum_right_m_1 * sum_right_m_1 / (double) n_right;
    
    // MCE.
    double mce = (prod_left / n_left) - (meanLeft_m * meanLeft_m_1) + (prod_right / n_right) - (meanRight_m * meanRight_m_1);
    
    // Total decrease in loss function.
    // double decrease = mse_m + mse_m_1 - 2 * mce;
    double decrease = mse_m + mse_m_1 + 2 * mce;
    
    // Regularization.
    regularize(decrease, varID);

    // If better than before, use this.
    if (decrease > best_decrease) {
      best_value = splitID;
      best_varID = varID;
      best_decrease = decrease;
    }
  }
}

void TreeOrdered::addImpurityImportance(size_t nodeID, size_t varID, double decrease) {
  size_t num_samples_node = end_pos[nodeID] - start_pos[nodeID];
  double best_decrease = decrease;
  
  // No variable importance for no split variables.
  size_t tempvarID = data->getUnpermutedVarID(varID);

  // Subtracting if corrected importance and permuted variable, else adding.
  if (importance_mode == IMP_GINI_CORRECTED && varID >= data->getNumCols()) {
    (*variable_importance)[tempvarID] -= best_decrease;
  } else {
    (*variable_importance)[tempvarID] += best_decrease;
  }
} 

} // namespace ocf
