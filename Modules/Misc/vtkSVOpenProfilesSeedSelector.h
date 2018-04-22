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

#ifndef vtkSVOpenProfilesSeedSelector_h
#define vtkSVOpenProfilesSeedSelector_h

#include "vtkSVSeedSelector.h"
#include "vtkIdList.h"
#include "vtkSVRenderer.h"

//#include "vtkvmtkComputationalGeometryWin32Header.h"
#include "vtkvmtkWin32Header.h"
#include "vtkSVMiscModule.h" // For exports

#include "vtkUnstructuredGrid.h"

class vtkPolyData;
class vtkPoints;
class vtkIdList;
class vtkDataArray;

class VTKSVMISC_EXPORT vtkSVOpenProfilesSeedSelector : public vtkSVSeedSelector
{
  public:
  vtkTypeMacro(vtkSVOpenProfilesSeedSelector,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  vtkSetObjectMacro(SeedIds, vtkIdList);
  vtkGetObjectMacro(SeedIds, vtkIdList);

  static vtkSVOpenProfilesSeedSelector *New();

  protected:
  vtkSVOpenProfilesSeedSelector();
  ~vtkSVOpenProfilesSeedSelector();

  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  vtkIdList* SeedIds;

  vtkSVRenderer *SVRenderer;

  private:
  vtkSVOpenProfilesSeedSelector(const vtkSVOpenProfilesSeedSelector&);  // Not implemented.
  void operator=(const vtkSVOpenProfilesSeedSelector&);  // Not implemented.
};

#endif
