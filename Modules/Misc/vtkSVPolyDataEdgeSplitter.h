/*=========================================================================
 *
 * Copyright (c) 2014-2015 The Regents of the University of California.
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
 *  \class vtkSVPolyDataEdgeSplitter
 *  \brief Using a polydata centerlines, separate the polydata into regions
 *  based on the centerlines
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */

#ifndef vtkSVPolyDataEdgeSplitter_h
#define vtkSVPolyDataEdgeSplitter_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyData.h"
#include "vtkIdList.h"

#include "vtkSVGlobals.h"

#include "vtkSVMiscModule.h" // For export

class VTKSVMISC_EXPORT vtkSVPolyDataEdgeSplitter : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkSVPolyDataEdgeSplitter,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkSVPolyDataEdgeSplitter *New();

  //@{
  /// \brief Get/Set macro for merged centerlines
  vtkSetObjectMacro(SplitPointIds,vtkIdList);
  vtkGetObjectMacro(SplitPointIds,vtkIdList);
  //@}

  //@{
  /// \brief Get/Set macro for array name used by the filter. Must
  //  be present on the centerlines.
  vtkSetStringMacro(SplitPointsArrayName);
  vtkGetStringMacro(SplitPointsArrayName);
  //@}

protected:
  vtkSVPolyDataEdgeSplitter();
  ~vtkSVPolyDataEdgeSplitter();

  // Usual data generation method
  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *);

  int PrepFilter(); // Prep work.
  int RunFilter(); // Run filter operations.

  int SplitCellsAroundPoint(vtkPolyData *pd, int ptId);
  int SplitEdge(vtkPolyData *pd, int cellId, int ptId0, int ptId1);

  char *SplitPointsArrayName;

  vtkPolyData *WorkPd;

  vtkIdList *SplitPointIds;

  vtkCellArray *NewCells;

  vtkPoints *NewPoints;

  int SplitPointsArrayAdded;
  std::vector<int> CellBool;
  std::vector<std::vector<int> > SplitCellsInfo;

private:
  vtkSVPolyDataEdgeSplitter(const vtkSVPolyDataEdgeSplitter&);  // Not implemented.
  void operator=(const vtkSVPolyDataEdgeSplitter&);  // Not implemented.
};

#endif
