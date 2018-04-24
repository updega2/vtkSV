
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
 *  \class vtkSVParameterizeVolumeOnPolycube
 *  \brief Using a polydata centerlines, separate the polydata into regions
 *  based on the centerlines
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */

#ifndef vtkSVParameterizeVolumeOnPolycube_h
#define vtkSVParameterizeVolumeOnPolycube_h

#include "vtkSVParameterizationModule.h" // For exports

#include "vtkIdList.h"
#include "vtkMatrix4x4.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyData.h"
#include "vtkStructuredGrid.h"
#include "vtkUnstructuredGrid.h"

#include "vtkSVGlobals.h"

class VTKSVPARAMETERIZATION_EXPORT vtkSVParameterizeVolumeOnPolycube : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkSVParameterizeVolumeOnPolycube,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  static vtkSVParameterizeVolumeOnPolycube *New();

  //@{
  /// \brief Get/Set macro for surface polycube
  vtkSetObjectMacro(SurfaceOnPolycubePd,vtkPolyData);
  vtkGetObjectMacro(SurfaceOnPolycubePd,vtkPolyData);
  //@}

  //@{
  /// \brief Get/Set macro for surface polycube
  vtkSetObjectMacro(FinalHexMesh,vtkUnstructuredGrid);
  vtkGetObjectMacro(FinalHexMesh,vtkUnstructuredGrid);
  //@}

  //@{
  /// \brief Get/Set macro for surface polycube
  vtkSetObjectMacro(PolycubeUg,vtkUnstructuredGrid);
  vtkGetObjectMacro(PolycubeUg,vtkUnstructuredGrid);
  //@}

  //@{
  /// \brief Get/Set macro for array name used by the filter. Must
  //  be present on the centerlines.
  vtkSetStringMacro(GroupIdsArrayName);
  vtkGetStringMacro(GroupIdsArrayName);
  //@}

  static int GetRegions(vtkPolyData *pd, std::string arrayName,
                        std::vector<Region> &allRegions);

  static int GetCCWPoint(vtkPolyData *pd, const int pointId, const int cellId);
  static int GetCWPoint(vtkPolyData *pd, const int pointId, const int cellId);

  static int CheckBoundaryEdge(vtkPolyData *pd, std::string arrayName, const int cellId, const int pointId0, const int pointId1);

  static int GetPointEdgeCells(vtkPolyData *pd, std::string arrayName,
                               const int cellId, const int pointId,
                               vtkIdList *sameCells);

  static int FindPointMatchingValues(vtkPointSet *ps, std::string arrayName, vtkIdList *matchingVals, int &returnPtId);

  static int RotateGroupToGlobalAxis(vtkPolyData *pd,
                                     const int thresholdId,
                                     std::string arrayName,
                                     vtkPolyData *rotPd,
                                     vtkMatrix4x4 *rotMatrix0,
                                     vtkMatrix4x4 *rotMatrix1);
  static int InterpolateMapOntoTarget(vtkPolyData *sourceBasePd,
                                      vtkPolyData *targetPd,
                                      vtkPolyData *targetBasePd,
                                      vtkPolyData *mappedPd,
                                      std::string dataMatchingArrayName);


protected:
  vtkSVParameterizeVolumeOnPolycube();
  ~vtkSVParameterizeVolumeOnPolycube();

  // Usual data generation method
  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *) override;

  int GetInteriorPointMaps(vtkPolyData *pdWithAllInterior,
                           vtkPolyData *pdWithCleanInterior,
                           vtkPolyData *pdWithoutInterior,
                           std::vector<int> &ptMap,
                           std::vector<std::vector<int> > &invPtMap);
  int GetVolumePointMaps(vtkUnstructuredGrid *ugAll,
                         vtkUnstructuredGrid *ugClean,
                         std::vector<int> &ptMap,
                         std::vector<std::vector<int> > &invPtMap);
  int MapInteriorBoundary(vtkStructuredGrid *paraHexVolume,
                          vtkPolyData *mappedSurface,
                          const std::vector<int> ptMap);
  int FixInteriorBoundary(vtkPolyData *mappedSurface,
                          const std::vector<std::vector<int> > invPtMap);
  int FixVolume(vtkUnstructuredGrid *mappedVolume,
                vtkUnstructuredGrid *cleanVolume,
                const std::vector<int> ptMap);
  int SetControlMeshBoundaries(vtkUnstructuredGrid *mappedVolume,
                               vtkUnstructuredGrid *cleanVolume,
                               const std::vector<int> ptMap,
                               const std::vector<std::vector<int> > invPtMap);
  int MapVolume(vtkStructuredGrid *paraHexVolume,
                vtkPolyData *mappedSurface,
                vtkStructuredGrid *mappedVolume);
  int ConvertUGToSG(vtkUnstructuredGrid *ug,
                    vtkStructuredGrid *sg,
                    const int w_div, const int h_div, const int l_div);
  int GetPointConnectivity(vtkUnstructuredGrid *hexMesh,
                           std::vector<std::vector<int> > &ptEdgeNeighbors);
  int SmoothStructuredGrid(vtkStructuredGrid *hexMesh, const int iters);
  int SmoothUnstructuredGrid(vtkUnstructuredGrid *hexMesh, const int iters,
                             std::string fixedPointsArrayName);
  int RemoveInteriorCells(vtkPolyData *quadMesh);

  vtkPolyData *WorkPd;
  vtkPolyData *SurfaceOnPolycubePd;

  vtkUnstructuredGrid *PolycubeUg;
  vtkUnstructuredGrid *FinalHexMesh;

  int PrepFilter(); // Prep work.
  int RunFilter(); // Run filter operations.

  char *GroupIdsArrayName;

private:
  vtkSVParameterizeVolumeOnPolycube(const vtkSVParameterizeVolumeOnPolycube&);  // Not implemented.
  void operator=(const vtkSVParameterizeVolumeOnPolycube&);  // Not implemented.
};

#endif
