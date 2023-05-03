// Ignore in coverage report (not used in R package)
// #nocov start
#ifndef DATADOUBLE_H_
#define DATADOUBLE_H_

#include <vector>
#include <utility>

#include "globals.h"
#include "utility.h"
#include "Data.h"

namespace ocf {

class DataDouble: public Data {
public:
  DataDouble() = default;
  
  DataDouble(const DataDouble&) = delete;
  DataDouble& operator=(const DataDouble&) = delete;

  virtual ~DataDouble() override = default;

  double get_x(size_t row, size_t col) const override {
    // Using permuted data for corrected impurity importance.
    size_t col_permuted = col;
    
    if (col >= num_cols) {
      col = getUnpermutedVarID(col);
      row = getPermutedSampleID(row);
    }

    return x[col * num_rows + row];
  }

  double get_y(size_t row, size_t col) const override {
    return y[col * num_rows + row];
  }

  void reserveMemory(size_t y_cols) override {
    x.resize(num_cols * num_rows);
    y.resize(y_cols * num_rows);
  }

  void set_x(size_t col, size_t row, double value, bool& error) override {
    x[col * num_rows + row] = value;
  }

  void set_y(size_t col, size_t row, double value, bool& error) override {
    y[col * num_rows + row] = value;
  }

private:
  std::vector<double> x;
  std::vector<double> y;
};

} // namespace ocf

#endif /* DATADOUBLE_H_ */
// #nocov end
