// Author: Mingcheng Chen (linyufly@gmail.com)

#include "vtkMatrixMath.h"

#include "vtkSparseMatrix.h"

#include <cmath>
#include <cstdio>

#include <algorithm>

const double kEpsilon = 1e-8;

// Multiply A^tA with b.
void vtkMatrixMath::MultiplyATAb(
    const vtkSparseMatrix *a_trans, const vtkSparseMatrix &a,
    const double *b, double *c) {
  double *temp = new double[a.GetNumberOfRows()];

  a.MultiplyColumn(b, temp);
  a_trans->MultiplyColumn(temp, c);

  delete [] temp;
}

double vtkMatrixMath::InnerProduct(const double *a, const double *b, int n) {
  double result = 0.0;
  for (int c = 0; c < n; c++) {
    result += a[c] * b[c];
  }
  return result;
}

void vtkMatrixMath::Add(const double *a, const double *b, double beta, int n, double *c) {
  for (int i = 0; i < n; i++) {
    c[i] = a[i] + b[i] * beta;
  }
}


void vtkMatrixMath::ConjugateGradient(
    const vtkSparseMatrix &a, const double *b, int num_iterations, double *x) {
  vtkSparseMatrix *a_trans = a.Transpose();

  double *a_trans_b = new double[a_trans->GetNumberOfRows()];
  a_trans->MultiplyColumn(b, a_trans_b);

  // Solve a_trans * a * x = a_trans_b.
  double *r = new double[a_trans->GetNumberOfRows()];
  double *p = new double[a_trans->GetNumberOfRows()];

  double *temp = new double[a_trans->GetNumberOfRows()];

  // temp = A'A * x
  this->MultiplyATAb(a_trans, a, x, temp);

  // r = A'b - temp
  this->Add(a_trans_b, temp, -1.0, a_trans->GetNumberOfRows(), r);

  // p = r
  std::copy(r, r + a_trans->GetNumberOfRows(), p);

  // rs_old = r' * r
  double rs_old = InnerProduct(r, r, a_trans->GetNumberOfRows());

  if (sqrt(rs_old) < kEpsilon) {
    printf("The initial solution is good.\n");
    return;
  }

  for (int iteration = 0;
       iteration < num_iterations && iteration < a_trans->GetNumberOfRows();
       iteration++) {
    // temp = A'A * p
    this->MultiplyATAb(a_trans, a, p, temp);

    // alpha = rs_old / (p' * temp)
    double alpha = rs_old / InnerProduct(p, temp, a_trans->GetNumberOfRows());

    // x = x + alpha * p
    this->Add(x, p, alpha, a_trans->GetNumberOfRows(), x);

    // r = r - alpha * temp
    this->Add(r, temp, -alpha, a_trans->GetNumberOfRows(), r);

    // rs_new = r' * r
    double rs_new = InnerProduct(r, r, a_trans->GetNumberOfRows());

    // Traditionally, if norm(rs_new) is small enough, the iteration can stop.
    if (sqrt(rs_new) < kEpsilon) {
      /// DEBUG ///
      printf("rs_new = %.20lf\n", rs_new);

      break;
    }

    // p = r + (rs_new / rs_old) * p
    this->Add(r, p, rs_new / rs_old, a_trans->GetNumberOfRows(), p);

    // rs_old = rs_new
    rs_old = rs_new;

    /// DEBUG ///
    printf("  rs_old = %.20lf\n", rs_old);
  }

  /// DEBUG ///
  printf("rs_old = %.20lf\n", rs_old);

  delete [] r;
  delete [] p;
  delete [] temp;
  delete [] a_trans_b;
}

