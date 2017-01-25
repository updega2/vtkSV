// Author: Mingcheng Chen (linyufly@gmail.com)

#ifndef vtkMatrixMath_h
#define vtkMatrixMath_h

#include "vtkDataObject.h"

class vtkSparseMatrix;

class vtkMatrixMath : public vtkDataObject
{
 public:

  void MultiplyATAb(const vtkSparseMatrix *a_trans, const vtkSparseMatrix &a,
                    const double *b, double *c);

  double InnerProduct(const double *a, const double *b, int n);

  void Add(const double *a, const double *b, double beta, int n, double *c);

  void ConjugateGradient(
      const vtkSparseMatrix &a, const double *b, int num_iterations, double *x);
};

#endif
