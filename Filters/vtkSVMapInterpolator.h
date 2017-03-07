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


/** @file vtkSVMapInterpolator.h
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

#ifndef __vtkSVMapInterpolator_h
#define __vtkSVMapInterpolator_h

#include "vtkPolyDataAlgorithm.h"

#include "vtkEdgeTable.h"
#include "vtkFloatArray.h"
#include "vtkPolyData.h"

#include <complex>
#include <vector>

class vtkSVMapInterpolator : public vtkPolyDataAlgorithm
{
public:
  static vtkSVMapInterpolator* New();
  //vtkTypeRevisionMacro(vtkSVMapInterpolator, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Print statements used for debugging
  vtkGetMacro(Verbose, int);
  vtkSetMacro(Verbose, int);

  // Description:
  // Print statements used for dnumber of subdivisions
  vtkGetMacro(NumSourceSubdivisions, int);
  vtkSetMacro(NumSourceSubdivisions, int);

  // Functions to set up complex linear system for landmark constraint
  static int InterpolateMapOntoSource(vtkPolyData *mappedSourcePd,
                                      vtkPolyData *mappedTargetPd,
                                      vtkPolyData *originalTargetPd,
                                      vtkPolyData *sourceToTargetPd);
  static int GetTriangleUV(double f[3], double pt0[3], double pt1[3],
                           double pt2[3], double &a0, double &a1, double &a2);

  static int ComputeArea(double pt0[], double pt1[], double pt2[], double &area);
  static int PDCheckArrayName(vtkPolyData *pd, int datatype, std::string arrayname);

  // Setup and Check Functions
protected:
  vtkSVMapInterpolator();
  ~vtkSVMapInterpolator();

  // Usual data generation method
  int RequestData(vtkInformation *vtkNotUsed(request),
		  vtkInformationVector **inputVector,
		  vtkInformationVector *outputVector);

  // Main functions in filter
  int SubdivideAndInterpolate();

  int MatchBoundaries();
  int FindBoundary(vtkPolyData *pd, vtkIntArray *isBoundary);
  int MoveBoundaryPoints();
  int GetPointOnTargetBoundary(int targPtId, int srcCellId, double returnPt[]);
  int BoundaryPointsOnCell(vtkPolyData *pd, int srcCellId, vtkIdList *boundaryPts, vtkIntArray *isBoundary);
  int GetProjectedPoint(double pt0[], double pt1[], double projPt[], double returnPt[]);
  int GetClosestTwoPoints(vtkPolyData *pd, double projPt[], vtkIdList *boundaryPts, int &ptId0, int &ptId1);


private:
  vtkSVMapInterpolator(const vtkSVMapInterpolator&);  // Not implemented.
  void operator=(const vtkSVMapInterpolator&);  // Not implemented.

  int Verbose;
  int NumSourceSubdivisions;
  int HasBoundary;

  vtkPolyData *SourceS2Pd;
  vtkPolyData *TargetPd;
  vtkPolyData *TargetS2Pd;
  vtkPolyData *MappedPd;
  vtkPolyData *MappedS2Pd;

  vtkIntArray *SourceBoundary;
  vtkIntArray *TargetBoundary;
};

#endif


