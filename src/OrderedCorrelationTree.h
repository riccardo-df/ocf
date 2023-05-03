#ifndef TREEORDERED_H_
#define TREEORDERED_H_

#include <vector>

#include "globals.h"
#include "Tree.h"

namespace ocf {

class TreeOrdered: public Tree {
public:
  TreeOrdered() = default;

  // Creating from loaded forest.
  TreeOrdered(std::vector<std::vector<size_t>>& child_nodeIDs, std::vector<size_t>& split_varIDs,
      std::vector<double>& split_values);

  TreeOrdered(const TreeOrdered&) = delete;
  TreeOrdered& operator=(const TreeOrdered&) = delete;

  virtual ~TreeOrdered() override = default;

  void allocateMemory() override;

  double estimate(size_t nodeID);
  void computePermutationImportanceInternal(std::vector<std::vector<size_t>>* permutations);
  void appendToFileInternal(std::ofstream& file) override;

  double getPrediction(size_t sampleID) const {
    size_t terminal_nodeID = prediction_terminal_nodeIDs[sampleID];
    return (split_values[terminal_nodeID]);
  }

  size_t getPredictionTerminalNodeID(size_t sampleID) const {
    return prediction_terminal_nodeIDs[sampleID];
  }

private:
  bool splitNodeInternal(size_t nodeID, std::vector<size_t>& possible_split_varIDs) override;
  void createEmptyNodeInternal() override;

  double computePredictionAccuracyInternal(std::vector<double>* prediction_error_casewise) override;
  
  // Called by splitNodeInternal(). Sets split_varIDs and split_values.
  bool findBestSplit(size_t nodeID, std::vector<size_t>& possible_split_varIDs);
  void findBestSplitValueSmallQ(size_t nodeID, size_t varID, double sum_node_m, double sum_node_m_1, double prod_node, 
                                size_t num_samples_node,
      double& best_value, size_t& best_varID, double& best_decrease);
  void findBestSplitValueSmallQ(size_t nodeID, size_t varID, double sum_node_m, double sum_node_m_1, double prod_node,
                                size_t num_samples_node,
      double& best_value, size_t& best_varID, double& best_decrease, std::vector<double> possible_split_values,
      std::vector<double>& sums_m, std::vector<double>& sums_m_1, std::vector<double>& prods, std::vector<size_t>& counter_m);
  void findBestSplitValueLargeQ(size_t nodeID, size_t varID, double sum_node_m, double sum_node_m_1, double prod_node,
                                size_t num_samples_node,
      double& best_value, size_t& best_varID, double& best_decrease);
  void findBestSplitValueUnordered(size_t nodeID, size_t varID, double sum_node_m, double sum_node_m_1, double prod_node,
                                   size_t num_samples_node,
      double& best_value, size_t& best_varID, double& best_decrease);

  bool findBestSplitMaxstat(size_t nodeID, std::vector<size_t>& possible_split_varIDs);

  void addImpurityImportance(size_t nodeID, size_t varID, double decrease);

  double computePredictionMSE();

  void cleanUpInternal() override {
    counter.clear();
    counter.shrink_to_fit();
    
    sums_m.clear();
    sums_m.shrink_to_fit();
    sums_m_1.clear();
    sums_m_1.shrink_to_fit();
    
    prods.clear();
    prods.shrink_to_fit();
  }
 
  /* These variables are used for the splitting criterion. When changing these objects, recall to modify accordingly 
     the splitting functions in TreeOrdered.cpp and the cleanUpInternal() function above in this file. */
  std::vector<size_t> counter;
  
  std::vector<double> sums_m;
  std::vector<double> sums_m_1;
  
  std::vector<double> prods;
};

} // namespace ocf

#endif /* TREEORDERED_H_ */
