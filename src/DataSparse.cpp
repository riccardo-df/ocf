#include "DataSparse.h"

namespace ocf {

DataSparse::DataSparse(Eigen::SparseMatrix<double>& x, Rcpp::NumericMatrix& y, std::vector<std::string> variable_names, 
                       size_t num_rows, size_t num_cols) : x { }{
  this->x.swap(x);
  this->y = y;
  this->variable_names = variable_names;
  this->num_rows = num_rows;
  this->num_cols = num_cols;
  this->num_cols_no_snp = num_cols;
}

} // namespace ocf
