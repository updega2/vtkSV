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

/** @file vtkPlanarMapper.h
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

#ifndef vtkPlanarMapper_h
#define vtkPlanarMapper_h

#include "vtkPolyDataAlgorithm.h"

#include "vtkEdgeTable.h"
#include "vtkFloatArray.h"
#include "vtkPolyData.h"

class vtkPlanarMapper : public vtkPolyDataAlgorithm
{
public:
  static vtkPlanarMapper* New();
  vtkTypeRevisionMacro(vtkPlanarMapper, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // CG Update method
  vtkGetMacro(BoundaryType, int);
  vtkSetMacro(BoundaryType, int);

  // Boundary Corners
  vtkGetVector4Macro(BoundaryCorners, int);
  vtkSetVector4Macro(BoundaryCorners, int);

  // Boundary Corners
  vtkGetVector4Macro(BoundaryLengths, double);
  vtkSetVector4Macro(BoundaryLengths, double);

  // Axis of the object to use on orientation with sphee map
  vtkSetVector3Macro(ObjectXAxis, double);
  vtkSetVector3Macro(ObjectZAxis, double);

  // Description:
  // Internal ids array name, generated by GenerateIdFilter
  vtkGetStringMacro(InternalIdsArrayName);
  vtkSetStringMacro(InternalIdsArrayName);

  //MAP
  enum WEIGHT_TYPE
  {
    HARMONIC = 0,
    MEAN_VALUE,
    TUTTE
  };

  //BOUNDARY_TYPE
  enum BOUNDARY_TYPE
  {
    SQUARE = 0,
    CIRCLE
  };

  static int CheckSurface(vtkPolyData *pd);

  // Edge functions
  static int CreateEdgeTable(vtkPolyData *pd, vtkEdgeTable *edgeTable,
                             vtkFloatArray *edgeWeights,
                             vtkIntArray *edgeNeighbors,
                             vtkIntArray *isBoundary);
  static int ComputeEdgeWeight(vtkPolyData *pd, vtkIdType cellId,
                               vtkIdType neighborCellId,
                               vtkIdType p0, vtkIdType p1, double &weight);
  static int GetEdgeCotangentAngle(double pt0[3], double pt1[3], double pt2[3], double &angle);
  static int ComputeArea(double pt0[3], double pt1[3], double pt2[3], double &area);

  // Helper functions
  static int GetPointNeighbors(vtkIdType p0, vtkPolyData *pd, vtkIdList *pointNeighbors);
  static int RunLoopFind(vtkPolyData *pd, vtkIdType startPt, vtkIdType nextCell,
                         vtkPolyData *loop);
  static int CheckArrayExists(vtkPolyData *pd, int datatype, std::string arrayname);

  // Matrix functions
  static int InvertSystem(std::vector<std::vector<double> > &mat,
                          std::vector<std::vector<double> > &invMat);
  static int MatrixVectorMultiply(std::vector<std::vector<double> > &mat,
                                  std::vector<double> &inVec,
                                  std::vector<double> &outVec);
  static int PrintMatrix(std::vector<std::vector<double> > &mat);

protected:
  vtkPlanarMapper();
  ~vtkPlanarMapper();

  // Usual data generation method
  int RequestData(vtkInformation *vtkNotUsed(request),
		  vtkInformationVector **inputVector,
		  vtkInformationVector *outputVector);

  // Main functions in filter
  int PrepFilter();
  int RunFilter();
  int SetBoundaries();
  int GetBoundaryLoop();
  int FindBoundaries();
  int CalculateSquareEdgeLengths();
  int SetSquareBoundary();
  int SetCircleBoundary();
  int SetInternalNodes();
  int SolveSystem();

  // Point and edge wise functions using discrete laplace-beltrami

private:
  vtkPlanarMapper(const vtkPlanarMapper&);  // Not implemented.
  void operator=(const vtkPlanarMapper&);  // Not implemented.

  int RemoveInternalIds;

  vtkPolyData   *InitialPd;
  vtkPolyData   *WorkPd;
  vtkPolyData   *PlanarPd;
  vtkEdgeTable  *EdgeTable;
  vtkFloatArray *EdgeWeights;
  vtkIntArray   *EdgeNeighbors;
  vtkIntArray   *IsBoundary;
  vtkPolyData   *Boundaries;
  vtkPolyData   *BoundaryLoop;

  std::vector<std::vector<double> > AHarm;
  std::vector<std::vector<double> > ATutte;
  std::vector<double> Xu;
  std::vector<double> Xv;
  std::vector<double> Bu;
  std::vector<double> Bv;

  char *InternalIdsArrayName;

  int BoundaryType;
  int BoundaryCorners[4];
  double BoundaryLengths[4];
  double ObjectXAxis[3];
  double ObjectZAxis[3];

  double Lambda;
  double Mu;
};

#endif
