#ifndef DATA_H_
#define DATA_H_

#include <vector>
#include <iostream>
#include <numeric>
#include <random>
#include <algorithm>

#include "globals.h"

namespace ocf {

class Data {
public:
  Data();

  Data(const Data&) = delete;
  Data& operator=(const Data&) = delete;

  virtual ~Data() = default;

  virtual double get_x(size_t row, size_t col) const = 0;
  virtual double get_y(size_t row, size_t col) const = 0;

  size_t getVariableID(const std::string& variable_name) const;

  virtual void reserveMemory(size_t y_cols) = 0;

  virtual void set_x(size_t col, size_t row, double value, bool& error) = 0;
  virtual void set_y(size_t col, size_t row, double value, bool& error) = 0;

  bool loadFromFile(std::string filename, std::vector<std::string>& dependent_variable_names);
  bool loadFromFileWhitespace(std::ifstream& input_file, std::string header_line,
      std::vector<std::string>& dependent_variable_names);
  bool loadFromFileOther(std::ifstream& input_file, std::string header_line,
      std::vector<std::string>& dependent_variable_names, char seperator);

  void getAllValues(std::vector<double>& all_values, std::vector<size_t>& sampleIDs, size_t varID, size_t start,
      size_t end) const;

  void getMinMaxValues(double& min, double&max, std::vector<size_t>& sampleIDs, size_t varID, size_t start,
      size_t end) const;

  size_t getIndex(size_t row, size_t col) const {
    // Using permuted data for corrected impurity importance.
    if (col >= num_cols) {
      col = getUnpermutedVarID(col);
      row = getPermutedSampleID(row);
    }

    return index_data[col * num_rows + row];
  }

  double getUniqueDataValue(size_t varID, size_t index) const {
    // Using permuted data for corrected impurity importance.
    if (varID >= num_cols) {
      varID = getUnpermutedVarID(varID);
    }

    return unique_data_values[varID][index];
  }

  size_t getNumUniqueDataValues(size_t varID) const {
    // Using permuted data for corrected impurity importance
    if (varID >= num_cols) {
      varID = getUnpermutedVarID(varID);
    }

    return unique_data_values[varID].size();
  }

  void sort();

  const std::vector<std::string>& getVariableNames() const {
    return variable_names;
  }
  size_t getNumCols() const {
    return num_cols;
  }
  size_t getNumRows() const {
    return num_rows;
  }

  size_t getMaxNumUniqueValues() const {
    return max_num_unique_values;
  }

  std::vector<bool>& getIsOrderedVariable() noexcept {
    return is_ordered_variable;
  }

  void setIsOrderedVariable(const std::vector<std::string>& unordered_variable_names) {
    is_ordered_variable.resize(num_cols, true);
    for (auto& variable_name : unordered_variable_names) {
      size_t varID = getVariableID(variable_name);
      is_ordered_variable[varID] = false;
    }
  }

  void setIsOrderedVariable(std::vector<bool>& is_ordered_variable) {
    this->is_ordered_variable = is_ordered_variable;
  }

  bool isOrderedVariable(size_t varID) const {
    // Using permuted data for corrected impurity importance.
    if (varID >= num_cols) {
      varID = getUnpermutedVarID(varID);
    }
    
    return is_ordered_variable[varID];
  }

  void permuteSampleIDs(std::mt19937_64 random_number_generator) {
    permuted_sampleIDs.resize(num_rows);
    std::iota(permuted_sampleIDs.begin(), permuted_sampleIDs.end(), 0);
    std::shuffle(permuted_sampleIDs.begin(), permuted_sampleIDs.end(), random_number_generator);
  }

  size_t getPermutedSampleID(size_t sampleID) const {
    return permuted_sampleIDs[sampleID];
  }

  size_t getUnpermutedVarID(size_t varID) const {
    if (varID >= num_cols) {
      varID -= num_cols;
    }
    
    return varID;
  }

protected:
  std::vector<std::string> variable_names;
  size_t num_rows;
  size_t num_rows_rounded;
  size_t num_cols;

  unsigned char* snp_data;
  size_t num_cols_no_snp;

  bool externalData;

  std::vector<size_t> index_data;
  std::vector<std::vector<double>> unique_data_values;
  size_t max_num_unique_values;

  // For each varID, true if ordered.
  std::vector<bool> is_ordered_variable;

  // Permuted samples for corrected impurity importance.
  std::vector<size_t> permuted_sampleIDs;

  // Order of 0/1/2 for ordered splitting.
  std::vector<std::vector<size_t>> snp_order;
  bool order_snps;
};

} // namespace ocf

#endif /* DATA_H_ */
