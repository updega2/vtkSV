/*=========================================================================

Program:   VMTK
Module:    $RCSfile: vtkSVCenterlines.h,v $
Language:  C++
Date:      $Date: 2006/07/17 09:52:56 $
Version:   $Revision: 1.6 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
  // .NAME vtkSVCenterlines - Compute centerlines from surface.
  // .SECTION Description
  // ...

#ifndef __vtkSVCenterlines_h
#define __vtkSVCenterlines_h

#include "vtkPolyDataAlgorithm.h"
//#include "vtkvmtkComputationalGeometryWin32Header.h"
#include "vtkvmtkWin32Header.h"
#include "vtkSVSegmentationModule.h" // For exports

#include "vtkUnstructuredGrid.h"

class vtkPolyData;
class vtkPoints;
class vtkIdList;
class vtkDataArray;

class VTKSVSEGMENTATION_EXPORT vtkSVCenterlines : public vtkPolyDataAlgorithm
{
  public:
  vtkTypeMacro(vtkSVCenterlines,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkSVCenterlines *New();

  virtual void SetSourceSeedIds(vtkIdList*);
  vtkGetObjectMacro(SourceSeedIds,vtkIdList);

  virtual void SetTargetSeedIds(vtkIdList*);
  vtkGetObjectMacro(TargetSeedIds,vtkIdList);

  virtual void SetCapCenterIds(vtkIdList*);
  vtkGetObjectMacro(CapCenterIds,vtkIdList);

  vtkSetObjectMacro(DelaunayTessellation,vtkUnstructuredGrid);
  vtkGetObjectMacro(DelaunayTessellation,vtkUnstructuredGrid);

  vtkGetObjectMacro(VoronoiDiagram,vtkPolyData);

  vtkGetObjectMacro(PoleIds,vtkIdList);

  vtkSetStringMacro(RadiusArrayName);
  vtkGetStringMacro(RadiusArrayName);

  vtkSetStringMacro(CostFunction);
  vtkGetStringMacro(CostFunction);

  vtkSetStringMacro(EikonalSolutionArrayName);
  vtkGetStringMacro(EikonalSolutionArrayName);

  vtkSetStringMacro(EdgeArrayName);
  vtkGetStringMacro(EdgeArrayName);

  vtkSetStringMacro(EdgePCoordArrayName);
  vtkGetStringMacro(EdgePCoordArrayName);

  vtkSetStringMacro(CostFunctionArrayName);
  vtkGetStringMacro(CostFunctionArrayName);

  vtkSetMacro(FlipNormals,int);
  vtkGetMacro(FlipNormals,int);
  vtkBooleanMacro(FlipNormals,int);

  vtkSetMacro(SimplifyVoronoi,int);
  vtkGetMacro(SimplifyVoronoi,int);
  vtkBooleanMacro(SimplifyVoronoi,int);

  vtkSetMacro(CenterlineResampling,int);
  vtkGetMacro(CenterlineResampling,int);
  vtkBooleanMacro(CenterlineResampling,int);

  vtkSetMacro(ResamplingStepLength,double);
  vtkGetMacro(ResamplingStepLength,double);

  vtkSetMacro(AppendEndPointsToCenterlines,int);
  vtkGetMacro(AppendEndPointsToCenterlines,int);
  vtkBooleanMacro(AppendEndPointsToCenterlines,int);

  vtkSetMacro(GenerateDelaunayTessellation,int);
  vtkGetMacro(GenerateDelaunayTessellation,int);
  vtkBooleanMacro(GenerateDelaunayTessellation,int);

  vtkSetMacro(DelaunayTolerance,double);
  vtkGetMacro(DelaunayTolerance,double);


  protected:
  vtkSVCenterlines();
  ~vtkSVCenterlines();

  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  void FindVoronoiSeeds(vtkUnstructuredGrid *delaunay, vtkIdList *boundaryBaricenterIds, vtkDataArray *normals, vtkIdList *seedIds);
  void AppendEndPoints(vtkPoints* endPointPairs);
  void ResampleCenterlines();
  void ReverseCenterlines();

  vtkIdList* SourceSeedIds;
  vtkIdList* TargetSeedIds;

  vtkIdList* CapCenterIds;

  vtkUnstructuredGrid* DelaunayTessellation;

  vtkPolyData* VoronoiDiagram;

  vtkIdList* PoleIds;

  char* RadiusArrayName;
  char* CostFunction;
  char* EikonalSolutionArrayName;
  char* EdgeArrayName;
  char* EdgePCoordArrayName;
  char* CostFunctionArrayName;

  int FlipNormals;
  int SimplifyVoronoi;
  int AppendEndPointsToCenterlines;
  int CenterlineResampling;

  double ResamplingStepLength;

  int GenerateDelaunayTessellation;
  double DelaunayTolerance;

  private:
  vtkSVCenterlines(const vtkSVCenterlines&);  // Not implemented.
  void operator=(const vtkSVCenterlines&);  // Not implemented.
};

#endif
