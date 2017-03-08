// Author: Mingcheng Chen (linyufly@gmail.com)
/*=========================================================================
 *
 * Copyright (c) 2014 The Regents of the University of California.
 * All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *=========================================================================*/

#ifndef vtkSparseMatrix_h
#define vtkSparseMatrix_h

#include "vtkDataObject.h"

#include <vector>

class vtkSparseMatrix : public vtkDataObject
{
 public:
  static vtkSparseMatrix* New();
  vtkTypeMacro(vtkSparseMatrix, vtkDataObject);
  void PrintSelf(ostream &os, vtkIndent indent);

  // Description:
  // Return what type of dataset this is.
  int GetDataObjectType() {return VTK_TABLE;}

  void MultiplyColumn(const double *column, double *result) const;
  void SetElement(int row, int col, double value);

  double GetElement(int row, int col) const;

  vtkSparseMatrix* Transpose() const;

  int GetNumberOfElements() const;

  void SetNumberOfRows(int num_rows);
  int GetNumberOfRows() const { return num_rows_;}

  void SetNumberOfColumns(int num_cols);
  int GetNumberOfColumns() const { return num_cols_;}
protected:
  vtkSparseMatrix();
  ~vtkSparseMatrix();

 private:
  std::vector<std::vector<double> > data_;
  std::vector<std::vector<int> > col_;
  int num_rows_, num_cols_;
};

#endif
