#ifndef DATASPARSE_H_
#define DATASPARSE_H_

#include <RcppEigen.h>

#include "globals.h"
#include "utility.h"
#include "Data.h"

namespace ocf {

class DataSparse: public Data {
public:
  DataSparse() = default;
  
  DataSparse(Eigen::SparseMatrix<double>& x, Rcpp::NumericMatrix& y, std::vector<std::string> variable_names, size_t num_rows,
      size_t num_cols);

  DataSparse(const DataSparse&) = delete;
  DataSparse& operator=(const DataSparse&) = delete;

  virtual ~DataSparse() override = default;

  double get_x(size_t row, size_t col) const override {
    // Using permuted data for corrected impurity importance.
    if (col >= num_cols) {
      col = getUnpermutedVarID(col);
      row = getPermutedSampleID(row);
    }
    
    return x.coeff(row, col);
  }
  
  double get_y(size_t row, size_t col) const override {
    return y[col * num_rows + row];
  }

  // #nocov start 
  void reserveMemory(size_t y_cols) override {
    // Not needed
  }

  void set_x(size_t col, size_t row, double value, bool& error) override {
    x.coeffRef(row, col) = value;
  }
  
  void set_y(size_t col, size_t row, double value, bool& error) override {
    y[col * num_rows + row] = value;
  }
  // #nocov end 

private:
  Eigen::SparseMatrix<double> x;
  Rcpp::NumericMatrix y;
};

} // namespace ocf

#endif /* DATASPARSE_H_ */
