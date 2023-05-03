#ifndef FOREST_H_
#define FOREST_H_

#include <vector>
#include <iostream>
#include <random>
#include <ctime>
#include <memory>
#ifndef OLD_WIN_R_BUILD
#include <thread>
#include <chrono>
#include <mutex>
#include <condition_variable>
#endif

#include "globals.h"
#include "Tree.h"
#include "Data.h"

namespace ocf {

class Forest {
public:
  Forest();

  Forest(const Forest&) = delete;
  Forest& operator=(const Forest&) = delete;

  virtual ~Forest() = default;

  // Initializing from c++ main or Rcpp from R.
  void initCpp(std::string dependent_variable_name, MemoryMode memory_mode, std::string input_file, unsigned int mtry,
      std::string output_prefix, unsigned int num_trees, std::ostream* verbose_out, unsigned int seed, unsigned int num_threads,
      std::string load_forest_filename, ImportanceMode importance_mode, unsigned int min_node_size,
      std::string split_select_weights_file, const std::vector<std::string>& always_split_variable_names,
      std::string status_variable_name, bool sample_with_replacement,
      const std::vector<std::string>& unordered_variable_names, bool memory_saving_splitting, SplitRule splitrule,
      std::string case_weights_file, bool predict_all, double sample_fraction, double alpha, double minprop,
      bool holdout, PredictionType prediction_type, unsigned int num_random_splits, unsigned int max_depth,
      const std::vector<double>& regularization_factor, bool regularization_usedepth,
      std::vector<double>& alpha_balance);
  
  void initR(std::unique_ptr<Data> input_data, unsigned int mtry, unsigned int num_trees, std::ostream* verbose_out, unsigned int seed,
      unsigned int num_threads, ImportanceMode importance_mode, unsigned int min_node_size,
      std::vector<std::vector<double>>& split_select_weights,
      const std::vector<std::string>& always_split_variable_names, bool prediction_mode, bool sample_with_replacement,
      const std::vector<std::string>& unordered_variable_names, bool memory_saving_splitting, SplitRule splitrule,
      std::vector<double>& case_weights, std::vector<std::vector<size_t>>& manual_inbag, bool predict_all,
      bool keep_inbag, std::vector<double>& sample_fraction, double alpha, double minprop, bool holdout,
      PredictionType prediction_type, unsigned int num_random_splits, bool order_snps, unsigned int max_depth,
      const std::vector<double>& regularization_factor, bool regularization_usedepth,
      std::vector<double>& alpha_balance);
  
  void init(std::unique_ptr<Data> input_data, unsigned int mtry, std::string output_prefix,
      unsigned int num_trees, unsigned int seed, unsigned int num_threads, ImportanceMode importance_mode, unsigned int min_node_size,
      bool prediction_mode, bool sample_with_replacement, const std::vector<std::string>& unordered_variable_names,
      bool memory_saving_splitting, SplitRule splitrule, bool predict_all, std::vector<double>& sample_fraction,
      double alpha, double minprop, bool holdout, PredictionType prediction_type, unsigned int num_random_splits,
      bool order_snps, unsigned int max_depth, const std::vector<double>& regularization_factor, bool regularization_usedepth,
      std::vector<double>& alpha_balance);
  virtual void initInternal() = 0;

  // Growing or predicting.
  void run(bool verbose, bool compute_oob_error);

  // Writing results to output files.
  void writeOutput();
  virtual void writeOutputInternal() = 0;
  virtual void writeConfusionFile() = 0;
  virtual void writePredictionFile() = 0;
  void writeImportanceFile();

  // Saving forest to file.
  void saveToFile();
  virtual void saveToFileInternal(std::ofstream& outfile) = 0;

  std::vector<std::vector<std::vector<size_t>>> getChildNodeIDs() {
    std::vector<std::vector<std::vector<size_t>>> result;
    for (auto& tree : trees) {
      result.push_back(tree->getChildNodeIDs());
    }
    
    return result;
  }
  
  std::vector<std::vector<size_t>> getSplitVarIDs() {
    std::vector<std::vector<size_t>> result;
    for (auto& tree : trees) {
      result.push_back(tree->getSplitVarIDs());
    }
    return result;
  }
  
  std::vector<std::vector<double>> getSplitValues() {
    std::vector<std::vector<double>> result;
    for (auto& tree : trees) {
      result.push_back(tree->getSplitValues());
    }
    
    return result;
  }
  
  const std::vector<double>& getVariableImportance() const {
    return variable_importance;
  }
  
  const std::vector<double>& getVariableImportanceCasewise() const {
    return variable_importance_casewise;
  }
  
  double getOverallPredictionError() const {
    return overall_prediction_error;
  }
  
  const std::vector<std::vector<std::vector<double>>>& getPredictions() const {
    return predictions;
  }
  
  size_t getNumTrees() const {
    return num_trees;
  }
  
  unsigned int getMtry() const {
    return mtry;
  }
  
  unsigned int getMinNodeSize() const {
    return min_node_size;
  }
  
  size_t getNumIndependentVariables() const {
    return num_independent_variables;
  }

  const std::vector<bool>& getIsOrderedVariable() const {
    return data->getIsOrderedVariable();
  }

  std::vector<std::vector<size_t>> getInbagCounts() const {
    std::vector<std::vector<size_t>> result;
    for (auto& tree : trees) {
      result.push_back(tree->getInbagCounts());
    }
    
    return result;
  }

protected:
  void grow();
  virtual void growInternal() = 0;

  // Predicting using existing tree from file and data as prediction data.
  void predict();
  virtual void allocatePredictMemory() = 0;
  virtual void predictInternal(size_t sample_idx) = 0;

  void computePredictionError();
  virtual void computePredictionErrorInternal() = 0;

  void computePermutationImportance();

  // Multithreading methods for growing/prediction/importance, called by each thread.
  void growTreesInThread(unsigned int thread_idx, std::vector<double>* variable_importance);
  void predictTreesInThread(unsigned int thread_idx, const Data* prediction_data, bool oob_prediction);
  void predictInternalInThread(unsigned int thread_idx);
  void computeTreePermutationImportanceInThread(unsigned int thread_idx, std::vector<double>& importance,
      std::vector<double>& variance, std::vector<double>& importance_casewise);

  // Loading forest from file.
  void loadFromFile(std::string filename);
  virtual void loadFromFileInternal(std::ifstream& infile) = 0;
  void loadDependentVariableNamesFromFile(std::string filename);

  // Loading data from file.
  std::unique_ptr<Data> loadDataFromFile(const std::string& data_path);

  // Setting split select weights and variables to be always considered for splitting.
  void setSplitWeightVector(std::vector<std::vector<double>>& split_select_weights);
  void setAlwaysSplitVariables(const std::vector<std::string>& always_split_variable_names);

  // Showing progress every few seconds.
#ifdef OLD_WIN_R_BUILD
  void showProgress(std::string operation, clock_t start_time, clock_t& lap_time);
#else
  void showProgress(std::string operation, size_t max_progress);
#endif

  // Verbose output stream, cout if verbose==true, logfile if not.
  std::ostream* verbose_out;

  std::vector<std::string> dependent_variable_names; 
  size_t num_trees;
  unsigned int mtry;
  unsigned int min_node_size;
  size_t num_independent_variables;
  unsigned int seed;
  size_t num_samples;
  bool prediction_mode;
  MemoryMode memory_mode;
  bool sample_with_replacement;
  bool memory_saving_splitting;
  SplitRule splitrule;
  bool predict_all;
  bool keep_inbag;
  std::vector<double> sample_fraction;
  bool holdout;
  PredictionType prediction_type;
  unsigned int num_random_splits;
  unsigned int max_depth;

  // MAXSTAT splitrule.
  double alpha;
  double minprop;

  // Multithreading.
  unsigned int num_threads;
  std::vector<unsigned int> thread_ranges;
#ifndef OLD_WIN_R_BUILD
  std::mutex mutex;
  std::condition_variable condition_variable;
#endif

  std::vector<std::unique_ptr<Tree>> trees;
  std::unique_ptr<Data> data;

  std::vector<std::vector<std::vector<double>>> predictions;
  double overall_prediction_error;

  /* Weight vector for selecting possible split variables, one weight between 0 (never select) and 1 (always select) 
     for each variable. Deterministic variables are always selected. */
  std::vector<size_t> deterministic_varIDs;
  std::vector<std::vector<double>> split_select_weights;

  // Bootstrap weights.
  std::vector<double> case_weights;

  // Pre-selected bootstrap samples (per tree).
  std::vector<std::vector<size_t>> manual_inbag;

  // Random number generator.
  std::mt19937_64 random_number_generator;

  std::string output_prefix;
  ImportanceMode importance_mode;

  // Regularization.
  std::vector<double> regularization_factor;
  bool regularization_usedepth;
  std::vector<bool> split_varIDs_used;
  
  // Alpha-regularity.
  std::vector<double> alpha_balance;
  
  // Variable importance for all variables in forest.
  std::vector<double> variable_importance;

  // Casewise variable importance for all variables in forest.
  std::vector<double> variable_importance_casewise;

  // Computation progress (finished trees).
  size_t progress;
#ifdef R_BUILD
  size_t aborted_threads;
  bool aborted;
#endif
};

} // namespace ocf

#endif /* FOREST_H_ */
