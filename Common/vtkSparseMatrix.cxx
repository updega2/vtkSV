
// Author: Mingcheng Chen (linyufly@gmail.com)

#include "vtkSparseMatrix.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"

#include <algorithm>

vtkStandardNewMacro(vtkSparseMatrix);

vtkSparseMatrix::vtkSparseMatrix()
{
  num_rows_ = 0;
  num_cols_ = 0;
}

vtkSparseMatrix::~vtkSparseMatrix()
{
}

void vtkSparseMatrix::SetNumberOfRows(int num_rows)
{
  this->num_rows_ = num_rows;
  data_.resize(num_rows);
}

void vtkSparseMatrix::SetNumberOfColumns(int num_cols)
{
  this->num_cols_ = num_cols;
  col_.resize(num_cols);
}

void vtkSparseMatrix::PrintSelf(ostream &os, vtkIndent indent)
{
  vtkDataObject::PrintSelf(os, indent);
}

void vtkSparseMatrix::SetElement(int row, int col, double value) {
  if (value == 0.0) {
    for (int c = 0; c < col_[row].size(); c++) {
      if (col_[row][c] == col) {
        col_[row].erase(col_[row].begin() + c);
        data_[row].erase(data_[row].begin() + c);
        break;
      }
    }
    return;
  }

  for (int c = 0; c < col_[row].size(); c++) {
    if (col_[row][c] == col) {
      data_[row][c] = value;
      return;
    }
  }

  col_[row].push_back(col);
  data_[row].push_back(value);
}

double vtkSparseMatrix::GetElement(int row, int col) const {
  for (int c = 0; c < col_[row].size(); c++) {
    if (col_[row][c] == col) {
      return data_[row][c];
    }
  }

  return 0.0;
}

vtkSparseMatrix* vtkSparseMatrix::Transpose() const {
  vtkSparseMatrix* result = vtkSparseMatrix::New();
  result->SetNumberOfRows(num_cols_);
  result->SetNumberOfColumns(num_rows_);

  for (int r = 0; r < num_rows_; r++) {
    for (int c = 0; c < col_[r].size(); c++) {
      result->SetElement(col_[r][c], r, data_[r][c]);
    }
  }

  return result;
}

int vtkSparseMatrix::GetNumberOfElements() const {
  int result = 0;
  for (int r = 0; r < num_rows_; r++) {
    result += col_[r].size();
  }

  return result;
}

void vtkSparseMatrix::MultiplyColumn(
    const double *column, double *result) const {
  for (int r = 0; r < num_rows_; r++) {
    result[r] = 0.0;
    for (int c = 0; c < col_[r].size(); c++) {
      result[r] += data_[r][c] * column[col_[r][c]];
    }
  }
}

