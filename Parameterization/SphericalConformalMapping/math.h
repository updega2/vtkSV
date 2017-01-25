// Author: Mingcheng Chen (linyufly@gmail.com)

#ifndef MATH_H_
#define MATH_H_

class SparseMatrix;

class Math {
 public:
  static void conjugate_gradient(
      const SparseMatrix &a, const double *b, int num_iterations, double *x);
  static void steepest_descent(
      const SparseMatrix &a, const double *b, int num_iterations, double *x, const double alpha);
};

#endif  // MATH_H_
