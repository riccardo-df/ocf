#ifndef DATARCPP_H_
#define DATARCPP_H_

#include <Rcpp.h>

#include "globals.h"
#include "utility.h"
#include "Data.h"

namespace ocf {

class DataRcpp: public Data {
public:
  DataRcpp() = default;
  DataRcpp(Rcpp::NumericMatrix& x, Rcpp::NumericMatrix& y, std::vector<std::string> variable_names, size_t num_rows, size_t num_cols) {
      this->x = x;
      this->y = y;
      this->variable_names = variable_names;
      this->num_rows = num_rows;
      this->num_cols = num_cols;
      this->num_cols_no_snp = num_cols;
    }
  
  DataRcpp(const DataRcpp&) = delete;
  DataRcpp& operator=(const DataRcpp&) = delete;
  
  virtual ~DataRcpp() override = default;
  
  double get_x(size_t row, size_t col) const override {
    // Using permuted data for corrected impurity importance.
    size_t col_permuted = col;
    if (col >= num_cols) {
      col = getUnpermutedVarID(col);
      row = getPermutedSampleID(row);
    }
    
    return x(row, col);
  }
  
  double get_y(size_t row, size_t col) const override {
    return y(row, col);
  }
  
  // #nocov start 
  void reserveMemory(size_t y_cols) override {
    // Not needed
  }
  
  void set_x(size_t col, size_t row, double value, bool& error) override {
    x(row, col) = value;
  }
  
  void set_y(size_t col, size_t row, double value, bool& error) override {
    y(row, col) = value;
  }
  // #nocov end 
  
private:
  Rcpp::NumericMatrix x;
  Rcpp::NumericMatrix y;
};

} // namespace ocf

#endif /* DATARCPP_H_ */
