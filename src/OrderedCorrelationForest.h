#ifndef FORESTORDERED_H_
#define FORESTORDERED_H_

#include <iostream>
#include <vector>

#include "globals.h"
#include "Forest.h"

namespace ocf {

class ForestOrdered: public Forest {
public:
  ForestOrdered() = default;

  ForestOrdered(const ForestOrdered&) = delete;
  ForestOrdered& operator=(const ForestOrdered&) = delete;

  virtual ~ForestOrdered() override = default;

  void loadForest(size_t num_trees, std::vector<std::vector<std::vector<size_t>> >& forest_child_nodeIDs,
      std::vector<std::vector<size_t>>& forest_split_varIDs, std::vector<std::vector<double>>& forest_split_values,
      std::vector<bool>& is_ordered_variable);

private:
  void initInternal() override;
  void growInternal() override;
  void allocatePredictMemory() override;
  void predictInternal(size_t sample_idx) override;
  void computePredictionErrorInternal() override;
  void writeOutputInternal() override;
  void writeConfusionFile() override;
  void writePredictionFile() override;
  void saveToFileInternal(std::ofstream& outfile) override;
  void loadFromFileInternal(std::ifstream& infile) override;

private:
  double getTreePrediction(size_t tree_idx, size_t sample_idx) const;
  size_t getTreePredictionTerminalNodeID(size_t tree_idx, size_t sample_idx) const;
};

} // namespace ocf

#endif /* FORESTORDERED_H_ */
