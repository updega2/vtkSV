/* Copyright (c) Stanford University, The Regents of the University of
 *               California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
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
 */

/**
 *  \file TestSparseMatrix.cxx
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */

#include "vtkSVSparseMatrix.h"

#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"

#include <cstdio>

#include <algorithm>
#include <vector>

static int TestTranspose()
{
  vtkNew(vtkSVSparseMatrix, a);
  a->SetMatrixSize(10, 10);

  a->SetElement(9, 1, 3.4);
  a->SetElement(2, 5, 3.5);
  a->SetElement(2, 6, 7.6);
  a->SetElement(3, 5, 88);
  a->SetElement(9, 1, 7.7);

  vtkNew(vtkSVSparseMatrix, aTrans);
  a->Transpose(aTrans);

  for (int i=0; i< a->GetNumberOfRows(); i++)
  {
    for (int j=0; j<a->GetNumberOfColumns(); j++)
    {
      if (a->GetElement(i, j) != aTrans->GetElement(j, i))
      {
        fprintf(stdout,"Transpose elements not equal\n");
        return SV_ERROR;
      }
    }
  }

  return SV_OK;
}

static int TestMultiplyColumn()
{
  vtkNew(vtkSVSparseMatrix, a);
  a->SetMatrixSize(3, 2);
  a->SetElement(0, 1, 3.0);
  a->SetElement(0, 0, 2.0);
  a->SetElement(2, 1, -3.0);
  a->SetElement(2, 0, -2.0);

  std::vector<double> column(2);
  column[0] = 2.0;
  column[1] = 3.0;

  std::vector<double> result(3);
  a->MultiplyColumn(&column[0], &result[0]);

  std::vector<double> expected_result(3);
  expected_result[0] = 13;
  expected_result[1] = 0;
  expected_result[2] = -13;

  for (int i=0; i<3; i++)
  {
    if (result[i] != expected_result[i])
    {
      fprintf(stdout,"Expected %.4f, but result was %.4f\n", expected_result[i], result[i]);
      return SV_ERROR;
    }
  }
  return SV_OK;
}

int TestSparseMatrix(int argc, char *argv[])
{
  if (TestTranspose() != SV_OK)
  {
    fprintf(stdout,"Transpose of matrix failed\n");
    return EXIT_FAILURE;
  }
  if (TestMultiplyColumn() != SV_OK)
  {
    fprintf(stdout,"Multiply column test failed\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

