// Author: Mingcheng Chen (linyufly@gmail.com)

#include "sparse_matrix.h"

#include <cstdio>

#include <algorithm>
#include <vector>

namespace {

void output_matrix(const SparseMatrix &matrix) {
  std::vector<double> column(matrix.get_num_cols());

  for (int r = 0; r < matrix.get_num_rows(); r++) {
    std::fill(column.begin(), column.end(), 0.0);

    for (int c = 0; c < matrix.get_num_cols(); c++) {
      column[c] = matrix.get_element(r, c);
    }

    for (int c = 0; c < matrix.get_num_cols(); c++) {
      printf(" %lf", column[c]);
    }

    printf("\n");
  }
}

}

void transpose_test() {
  printf("transpose_test {\n");

  SparseMatrix a(10, 10);

  a.set_element(9, 1, 3.4);
  a.set_element(2, 5, 3.5);
  a.set_element(2, 6, 7.6);
  a.set_element(3, 5, 88);
  a.set_element(9, 1, 7.7);

  output_matrix(a);

  printf("\n");

  output_matrix(a.transpose());

  printf("} transpose_test\n");
}

void multiply_column_test() {
  printf("multiply_column_test {\n");

  SparseMatrix a(3, 2);
  a.set_element(0, 1, 3.0);
  a.set_element(0, 0, 2.0);
  a.set_element(2, 1, -3.0);
  a.set_element(2, 0, -2.0);

  std::vector<double> column(2);
  column[0] = 2.0;
  column[1] = 3.0;

  std::vector<double> result(3);
  a.multiply_column(&column[0], &result[0]);

  for (int r = 0; r < 3; r++) {
    printf(" %lf", result[r]);
  }
  printf("\n");

  printf("} multiply_column_test\n");
}

int main() {
  transpose_test();
  multiply_column_test();

  return 0;
}

