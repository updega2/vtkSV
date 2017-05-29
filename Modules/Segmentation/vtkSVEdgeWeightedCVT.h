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
 *  \class  vtkSVGeneralCVT
 *  \brief This is class to perform edge weighted CVT clustering of an input polydata and
 *  generators
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */

#ifndef vtkSVEdgeWeightedCVT_h
#define vtkSVEdgeWeightedCVT_h

#include "vtkSVGeneralCVT.h"
#include "vtkSVSegmentationModule.h" // For exports

#include "vtkPolyData.h"

#include <vector>

class VTKSVSEGMENTATION_EXPORT vtkSVEdgeWeightedCVT : public vtkSVGeneralCVT
{
public:
  static vtkSVEdgeWeightedCVT* New();
  vtkTypeMacro(vtkSVEdgeWeightedCVT,vtkSVGeneralCVT);
  void PrintSelf(ostream& os, vtkIndent indent);

  //@{
  /// \brief Number of elements to be considered part of the element neighborhood.
  vtkGetMacro(NumberOfRings, int);
  vtkSetMacro(NumberOfRings, int);
  //@}
  //@{
  /// \brief Set a threshold criteria. Default is 2 transferred groups.
  vtkGetMacro(EdgeWeight, double);
  vtkSetMacro(EdgeWeight, double);
  //@}

protected:
  vtkSVEdgeWeightedCVT();
  ~vtkSVEdgeWeightedCVT();

  // Derived functions
  int InitializeConnectivity();
  int InitializeGenerators();
  int GetClosestGenerator(const int evalId, int &newGenerator);
  int ComputeSurfaceMetric(double &evalMetric);
  int UpdateConnectivity();
  int UpdateGenerators();

  // Edge weights setup
  int GetPointCellValence();
  int GetCellRingNeighbors(int ringNumber=1);
  int GetCellDirectNeighbors();
  int GetCellGroupNeighbors();
  int AddCellGroupNeighbor(const int cellId, const int cellNeighborGroup, int &neighborLoc);

private:
  vtkSVEdgeWeightedCVT(const vtkSVEdgeWeightedCVT&);  // Not implemented.
  void operator=(const vtkSVEdgeWeightedCVT&);  // Not implemented.

  std::vector<int>          PointCellValenceNumber;
  std::vector<std::vector<int> > PointCellValence;

  std::vector<int>          NumberOfNeighbors;
  std::vector<std::vector<int> > Neighbors;

  std::vector<int>          NumberOfDirectNeighbors;
  std::vector<std::vector<int> > DirectNeighbors;

  std::vector<int>          NumberOfNeighborGroups;
  std::vector<std::vector<int> > NeighborGroupsNumberOfElements;
  std::vector<std::vector<int> > NeighborGroupsIds;

  int NumberOfRings;
  double EdgeWeight;

};

#endif
