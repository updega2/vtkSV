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

/**
 *  \file ConugateGradient.cxx
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */

#include "vtkSVMathUtils.h"
#include "vtkSVSparseMatrix.h"

#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"

#include <algorithm>
#include <vector>

static int TestSolve()
{
  vtkNew(vtkSVSparseMatrix, a);
  a->SetMatrixSize(1, 2);
  a->SetElement(0, 0, 1.0);
  a->SetElement(0, 1, 2.0);

  std::vector<double> b(1);
  b[0] = 3.0;

  std::vector<double> x(2);
  x[0] = 7.0;
  x[1] = -3.0;

  vtkSVMathUtils::ConjugateGradient(a, &b[0], 1, &x[0], 1.0e-8);

  std::vector<double> test_b(4);
  a->MultiplyColumn(&x[0], &test_b[0]);

  if (b[0] != test_b[0])
  {
    fprintf(stdout,"Incorrect element\n");
    return SV_ERROR;
  }

  return SV_OK;
}

static int TestMinimize()
{
  vtkNew(vtkSVSparseMatrix, a);
  a->SetMatrixSize(4, 2);
  a->SetElement(0, 0, 1.0);
  a->SetElement(0, 1, 2.0);

  a->SetElement(1, 0, 2.0);
  a->SetElement(1, 1, -3.0);

  a->SetElement(2, 0, 4.0);
  a->SetElement(2, 1, -1.0);

  a->SetElement(3, 0, -5.0);
  a->SetElement(3, 1, 2.0);

  std::vector<double> b(4);
  b[0] = 3.0;  // 8.0
  b[1] = 2.0;  // -5.0
  b[2] = -4.0;
  b[3] = 5.0;

  std::vector<double> x(2);
  x[0] = x[1] = 0.0;

  double res = vtkSVMathUtils::ConjugateGradient(a, &b[0], 1000, &x[0], 1.0e-8);

  if (res > 1.0e-8)
  {
    fprintf(stdout,"System was not minimized to specified tolerance\n");
    return SV_ERROR;
  }
  std::vector<double> c(4);
  a->MultiplyColumn(&x[0], &c[0]);

  return SV_OK;
}

int TestConjugateGradient(int argc, char *argv[])
{
  if (TestSolve() != SV_OK)
  {
    fprintf(stdout,"Incorrect system solve\n");
    EXIT_FAILURE;
  }

  if (TestMinimize() != SV_OK)
  {
    fprintf(stdout,"Didn't minimize system\n");
    EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
