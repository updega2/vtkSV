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


/** @file vtkPullApartPolyData.h
 *  @brief This is a vtk filter to map a triangulated surface to a sphere.
 *  @details This filter uses the heat flow method to map a triangulated
 *  surface to a sphere. The first step is to compute the Tutte Energy, and
 *  the second step is to perform the conformal map. For more details, see
 *  Gu et al., Genus Zero Surface Conformal Mapping and Its
 *  Application to Brain Surface Mapping, 2004.
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#ifndef vtkPullApartPolyData_h
#define vtkPullApartPolyData_h

#include "vtkPolyDataAlgorithm.h"

#include "vtkEdgeTable.h"
#include "vtkDoubleArray.h"
#include "vtkPolyData.h"
#include "vtkIdList.h"

class vtkPullApartPolyData : public vtkPolyDataAlgorithm
{
public:
  static vtkPullApartPolyData* New();
  //vtkTypeRevisionMacro(vtkPullApartPolyData, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // String to separate the polydata at. 1 Indicates the points that are along
  // separation line, everything else should be 0
  vtkGetStringMacro(CutPointsArrayName);
  vtkSetStringMacro(CutPointsArrayName);

  // Description:
  // The list of points that are to be replaced
  vtkSetObjectMacro(SeamPointIds, vtkIntArray);
  vtkGetObjectMacro(SeamPointIds, vtkIntArray);

  // Description:
  // The list of points that are to be replaced
  vtkGetObjectMacro(ReplacePointList, vtkIdList);

  // Description:
  // The list that will be the same length as replacepoint list with the ids
  // corresponding to the new points
  vtkGetObjectMacro(NewPointList, vtkIdList);

  // Description:
  // Axis of the object to use on orientation with sphee map
  vtkSetVector3Macro(ObjectXAxis, double);
  vtkSetVector3Macro(ObjectZAxis, double);

  // Description:
  // If start point is provided, it helps the algorithm go quicker because
  // a start point does not need to be found
  vtkGetMacro(StartPtId, int);
  vtkSetMacro(StartPtId, int);

protected:
  vtkPullApartPolyData();
  ~vtkPullApartPolyData();

  // Usual data generation method
  int RequestData(vtkInformation *vtkNotUsed(request),
		  vtkInformationVector **inputVector,
		  vtkInformationVector *outputVector);

  int PrepFilter();
  int RunFilter();
  int FindEdgeCells();
  int PullApartCutEdges();
  int FindStartingEdge(int &p0, int &p1, int &p2, int &cellId);
  int FindNextEdge(int p0, int p1, int p2, int cellId, std::vector<int> &cellList, int first);
  int FixTheBadStartCell(vtkPolyData *pd, const int pointId, const int cellId);
  int CheckArrayExists(vtkPolyData *pd, int datatype, std::string arrayname);

private:
  vtkPullApartPolyData(const vtkPullApartPolyData&);  // Not implemented.
  void operator=(const vtkPullApartPolyData&);  // Not implemented.

  char *CutPointsArrayName;

  // TODO: Add start and end point ids

  int StartPtId;
  vtkPolyData  *WorkPd;
  vtkEdgeTable *EdgeTable;
  vtkIntArray  *SeamPointIds;
  vtkIdList    *ReplacePointList;
  vtkIdList    *NewPointList;

  double ObjectXAxis[3];
  double ObjectZAxis[3];

  std::vector<int> ReplacePointVector;
  std::vector<std::vector<int> > ReplaceCellVector;

};

#endif
