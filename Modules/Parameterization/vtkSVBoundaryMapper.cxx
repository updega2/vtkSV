/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSVBoundaryMapper.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSVBoundaryMapper.h"

#include "vtkCellIterator.h"
#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkConnectivityFilter.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkEdgeTable.h"
#include "vtkFeatureEdges.h"
#include "vtkIdFilter.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkMergePoints.h"
#include "vtkPointData.h"
#include "vtkPointLocator.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"
#include "vtkUnsignedCharArray.h"
#include "vtkUnstructuredGrid.h"

#include <sstream>
#include <map>

vtkSVBoundaryMapper::vtkSVBoundaryMapper()
{
  this->RemoveInternalIds = 1;

  this->InitialPd     = vtkPolyData::New();
  this->BoundaryPd    = vtkPolyData::New();
  this->EdgeTable     = NULL;

  this->IsBoundary    = NULL;
  this->BoundaryIds   = NULL;
  this->Boundaries    = vtkPolyData::New();
  this->BoundaryLoop  = vtkPolyData::New();

  this->InternalIdsArrayName = NULL;

  this->SetObjectXAxis(1.0, 0.0, 0.0);
  this->SetObjectZAxis(0.0, 0.0, 1.0);
}

//---------------------------------------------------------------------------
vtkSVBoundaryMapper::~vtkSVBoundaryMapper()
{
  if (this->InitialPd != NULL)
  {
    InitialPd->Delete();
  }
  if (this->BoundaryPd != NULL)
  {
    BoundaryPd->Delete();
  }
  if (this->EdgeTable != NULL)
  {
    EdgeTable->Delete();
  }
  if (this->Boundaries != NULL)
  {
    this->Boundaries->Delete();
  }
  if (this->BoundaryLoop != NULL)
  {
    this->BoundaryLoop->Delete();
  }
  if (this->InternalIdsArrayName)
  {
    delete [] this->InternalIdsArrayName;
    this->InternalIdsArrayName = NULL;
  }
}

int vtkSVBoundaryMapper::RequestData(vtkInformation *vtkNotUsed(request),
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector)
{
  // get the input and output
  vtkPolyData *input = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  //Copy the input to operate on
  this->InitialPd->DeepCopy(input);

  if (this->PrepFilter() != SV_OK)
  {
    vtkErrorMacro("Error when mapping");
    output->DeepCopy(this->InitialPd);
    return SV_ERROR;
  }

  if (this->RunFilter() != SV_OK)
  {
    vtkErrorMacro("Error when mapping");
    output->DeepCopy(this->InitialPd);
    return SV_ERROR;
  }

  if (this->RemoveInternalIds)
  {
    this->BoundaryPd->GetPointData()->RemoveArray(this->InternalIdsArrayName);
    this->BoundaryPd->GetCellData()->RemoveArray(this->InternalIdsArrayName);
  }
  output->DeepCopy(this->BoundaryPd);

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVBoundaryMapper::PrepFilter()
{
  vtkIdType numPolys = this->InitialPd->GetNumberOfPolys();
  vtkIdType numPoints = this->InitialPd->GetNumberOfPoints();
  //Check the input to make sure it is there
  if (numPolys < 1)
  {
    vtkDebugMacro("No input!");
    return SV_ERROR;
  }

  //Check the input to make sure it is manifold and a triangulated surface
  if (vtkSVGeneralUtils::CheckSurface(this->InitialPd) != SV_OK)
  {
    vtkErrorMacro("Error when checking input surface");
    return SV_ERROR;
  }

  // Check if internal id array name is given
  if (!this->InternalIdsArrayName)
  {
    vtkDebugMacro("Internal Ids Array Name not given, setting to InternalIds");
    this->InternalIdsArrayName = new char[strlen("InternalIds") + 1];
    strcpy(this->InternalIdsArrayName, "InternalIds");
  }
  // Check if array internal ids is already on pd
  if (vtkSVGeneralUtils::CheckArrayExists(this->InitialPd, 0, this->InternalIdsArrayName))
  {
    this->RemoveInternalIds = 0;
  }
  else
  {
    vtkNew(vtkIdFilter, ider);
    ider->SetInputData(this->InitialPd);
    ider->SetIdsArrayName(this->InternalIdsArrayName);
    ider->Update();
    this->InitialPd->DeepCopy(ider->GetOutput());
  }

  //Create the edge table for the input surface
  this->InitialPd->BuildLinks();

  if (this->EdgeTable->GetNumberOfEdges() == 0)
  {
    vtkErrorMacro("No Edges! Use SetEdgeTable");
    return SV_ERROR;
  }

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVBoundaryMapper::RunFilter()
{
  if (this->FindBoundaries() != SV_OK)
  {
    vtkErrorMacro("Could not find boundaries");
    return SV_ERROR;
  }

  if (this->GetBoundaryLoop() != SV_OK)
  {
    vtkErrorMacro("Error orienting boundary loop");
    return SV_ERROR;
  }

  if (this->SetBoundaries() != SV_OK)
  {
    vtkErrorMacro("Error in mapping");
    return SV_ERROR;
  }
  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
//Determine type of intersection
int vtkSVBoundaryMapper::GetBoundaryLoop()
{
  vtkIdType nextCell;
  vtkNew(vtkIdList, cellIds);
  vtkDataArray *pointIds = this->Boundaries->GetPointData()->GetArray(this->InternalIdsArrayName);
  vtkDataArray *oPointIds = this->InitialPd->GetPointData()->GetArray(this->InternalIdsArrayName);
  int numInterPts = this->Boundaries->GetNumberOfPoints();
  int numInterLines = this->Boundaries->GetNumberOfLines();
  this->Boundaries->BuildLinks();

  int count = 0;
  vtkIdType startPt = pointIds->LookupValue(
    oPointIds->GetTuple1(this->BoundaryIds->GetValue(0)));
  fprintf(stdout,"Start Point is!: %d\n", this->BoundaryIds->GetValue(0));
  this->BoundaryLoop->SetPoints(this->Boundaries->GetPoints());
  this->BoundaryLoop->GetPointData()->PassData(this->Boundaries->GetPointData());
  this->BoundaryLoop->Allocate(this->Boundaries->GetNumberOfCells(), 1000);
  fprintf(stdout,"The value on this is!: %lld\n", startPt);
  this->Boundaries->GetPointCells(startPt,cellIds);

  nextCell = cellIds->GetId(0);

  vtkNew(vtkIdList, boundaryIds);
  boundaryIds->SetNumberOfIds(this->BoundaryIds->GetNumberOfTuples());
  for (int i=0; i<this->BoundaryIds->GetNumberOfTuples(); i++)
    boundaryIds->SetId(i, pointIds->LookupValue(this->BoundaryIds->GetTuple1(i)));
  if (vtkSVGeneralUtils::RunLoopFind(this->Boundaries, startPt, nextCell, this->BoundaryLoop, boundaryIds) != SV_OK)
  {
    fprintf(stdout,"Other direction!\n");
    nextCell = cellIds->GetId(1);
    this->BoundaryLoop->DeleteCells();
    if (vtkSVGeneralUtils::RunLoopFind(this->Boundaries, startPt, nextCell, this->BoundaryLoop, boundaryIds) != SV_OK)
    {
      fprintf(stdout,"Both directions didn't work!!\n");
      return SV_ERROR;
    }
  }

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVBoundaryMapper::FindBoundaries()
{
  vtkIndent indenter;
  vtkNew(vtkPointLocator, locator);
  vtkNew(vtkFeatureEdges, finder);
  finder->SetInputData(this->InitialPd);
  finder->FeatureEdgesOff();
  //finder->SetLocator(locator);
  finder->NonManifoldEdgesOff();
  finder->BoundaryEdgesOn();
  finder->Update();

  vtkNew(vtkConnectivityFilter, connector);
  connector->SetInputData(finder->GetOutput());
  connector->SetExtractionMode(VTK_EXTRACT_ALL_REGIONS);
  connector->ColorRegionsOn();
  connector->Update();

  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(connector->GetOutput());
  surfacer->Update();

  this->Boundaries->DeepCopy(surfacer->GetOutput());

  if (this->Boundaries->GetNumberOfCells() == 0)
  {
    vtkErrorMacro("No boundaries on polydata");
    return SV_ERROR;
  }

  return SV_OK;
}


void vtkSVBoundaryMapper::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  //os << indent << "Number of subdivisions: "
  //   << this->GetNumberOfSubdivisions() << endl;
}
