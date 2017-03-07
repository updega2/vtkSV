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
#include "vtkUnsignedCharArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLPolyDataWriter.h"

#define vtkNew(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

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

  if (this->PrepFilter() != 1)
  {
    vtkErrorMacro("Error when mapping");
    output->DeepCopy(this->InitialPd);
    return 0;
  }

  if (this->RunFilter() != 1)
  {
    vtkErrorMacro("Error when mapping");
    output->DeepCopy(this->InitialPd);
    return 0;
  }

  if (this->RemoveInternalIds)
  {
    this->BoundaryPd->GetPointData()->RemoveArray(this->InternalIdsArrayName);
    this->BoundaryPd->GetCellData()->RemoveArray(this->InternalIdsArrayName);
  }
  output->DeepCopy(this->BoundaryPd);

  return 1;
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
    return 0;
  }

  //Check the input to make sure it is manifold and a triangulated surface
  if (this->CheckSurface(this->InitialPd) != 1)
  {
    vtkErrorMacro("Error when checking input surface");
    return 0;
  }

  // Check if internal id array name is given
  if (!this->InternalIdsArrayName)
  {
    vtkDebugMacro("Internal Ids Array Name not given, setting to InternalIds");
    this->InternalIdsArrayName = new char[strlen("InternalIds") + 1];
    strcpy(this->InternalIdsArrayName, "InternalIds");
  }
  // Check if array internal ids is already on pd
  if (this->CheckArrayExists(this->InitialPd, 0, this->InternalIdsArrayName))
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
    return 0;
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVBoundaryMapper::RunFilter()
{
  if (this->FindBoundaries() != 1)
  {
    vtkErrorMacro("Could not find boundaries");
    return 0;
  }

  if (this->GetBoundaryLoop() != 1)
  {
    vtkErrorMacro("Error orienting boundary loop");
    return 0;
  }

  if (this->SetBoundaries() != 1)
  {
    vtkErrorMacro("Error in mapping");
    return 0;
  }
  return 1;
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
  vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
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
  vtkIdType npts, *pts;
  int testPt = -1;
  this->Boundaries->GetCellPoints(nextCell, npts, pts);
  if (pts[0] == startPt)
    testPt = pts[1];
  else
    testPt = pts[0];
  fprintf(stdout,"And the next ponit: %f\n", pointIds->GetTuple1(testPt));

  double pt0[3], pt1[3], vec0[3], vec1[3];
  this->Boundaries->GetPoint(startPt, pt0);
  this->Boundaries->GetPoint(testPt, pt1);
  vtkMath::Subtract(pt1, pt0, vec0);
  vtkMath::Normalize(vec0);
  vtkMath::Cross(this->ObjectZAxis, this->ObjectXAxis, vec1);
  vtkMath::Normalize(vec1);

  vtkIdType checknpts, *checkpts;
  this->Boundaries->GetCellPoints(cellIds->GetId(1), checknpts, checkpts);
  int doubleCheckPt = -1;
  if (checkpts[0] == startPt)
    doubleCheckPt = checkpts[1];
  else
    doubleCheckPt = checkpts[0];
  double pt2[3], vec2[3], vec3[3];
  this->Boundaries->GetPoint(doubleCheckPt, pt2);
  vtkMath::Subtract(pt2, pt0, vec2);
  vtkMath::Normalize(vec2);
  vtkMath::Cross(this->ObjectZAxis, this->ObjectXAxis, vec3);
  vtkMath::Normalize(vec3);
  fprintf(stdout,"And the other ponit: %f\n", pointIds->GetTuple1(doubleCheckPt));
  if (vtkMath::Dot(vec0, vec1) < vtkMath::Dot(vec2, vec3))
  {
    fprintf(stdout,"Fliippped\n");
    nextCell = cellIds->GetId(1);
  }

  vtkSVBoundaryMapper::RunLoopFind(this->Boundaries, startPt, nextCell, this->BoundaryLoop);

  return 1;
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
    return 0;
  }

  return 1;
}


//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVBoundaryMapper::RunLoopFind(vtkPolyData *pd,
                                   vtkIdType startPt,
                                   vtkIdType nextCell,
                                   vtkPolyData *loop)
{
  vtkIdType prevPt = startPt;
  vtkIdType nextPt = startPt;
  vtkSmartPointer<vtkIdList> pointIds = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();

  pd->GetCellPoints(nextCell,pointIds);

  if (pointIds->GetId(0) == nextPt)
    nextPt = pointIds->GetId(1);
  else
    nextPt = pointIds->GetId(0);
  vtkSmartPointer<vtkIdList> newline = vtkSmartPointer<vtkIdList>::New();
  newline->SetNumberOfIds(2);
  newline->SetId(0, prevPt);
  newline->SetId(1, nextPt);
  //newline.id = nextCell;
  loop->InsertNextCell(VTK_LINE, newline);

  while(nextPt != startPt)
  {
    pd->GetPointCells(nextPt,cellIds);
    if (cellIds->GetId(0) == nextCell)
      nextCell = cellIds->GetId(1);
    else
      nextCell = cellIds->GetId(0);

    pd->GetCellPoints(nextCell,pointIds);
    prevPt = nextPt;
    if (pointIds->GetId(0) == nextPt)
      nextPt = pointIds->GetId(1);
    else
      nextPt = pointIds->GetId(0);

    vtkSmartPointer<vtkIdList> newestline = vtkSmartPointer<vtkIdList>::New();
    newestline->SetNumberOfIds(2);
    newestline->InsertId(0, prevPt);
    newestline->InsertId(1, nextPt);
    //newestline.id = nextCell;
    loop->InsertNextCell(VTK_LINE, newestline);
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVBoundaryMapper::CheckSurface(vtkPolyData *pd)
{
  pd->BuildLinks();

  int numPts = pd->GetNumberOfPoints();
  int numPolys = pd->GetNumberOfCells();

  for (int i=0; i<numPolys; i++)
  {
    vtkIdType npts, *pts;
    pd->GetCellPoints(i, npts, pts);
    if (npts != 3)
    {
      //vtkErrorMacro("Surface contains elements that aren't triangles");
      return 0;
    }
    for (int j=0; j<npts; j++)
    {
      vtkIdType p0, p1;
      p0 = pts[j];
      p1 = pts[(j+1)%npts];

      vtkSmartPointer<vtkIdList> edgeNeighbor =
        vtkSmartPointer<vtkIdList>::New();
      pd->GetCellEdgeNeighbors(i, p0, p1, edgeNeighbor);

      if (edgeNeighbor->GetNumberOfIds() > 1)
      {
        //vtkErrorMacro("Surface contains triangles with multiple neighbors, not manifold");
        return 0;
      }
    }
  }
  return 1;
}

// ----------------------
// CheckArrayExists
// ----------------------
/**
 * @brief Function to check is array with name exists in cell or point data
 * @param object this is the object to check if the array exists
 * @param datatype this is point or cell. point =0,cell=1
 * @param arrayname this is the name of the array to check
 * @reutrn this returns 1 if the array exists and zero if it doesn't
 * or the function does not return properly.
 */

int vtkSVBoundaryMapper::CheckArrayExists(vtkPolyData *object,int datatype,std::string arrayname )
{
  vtkIdType i;
  int numArrays;
  int exists =0;

  if (datatype == 0)
  {
    numArrays = object->GetPointData()->GetNumberOfArrays();
    for (i=0;i<numArrays;i++)
    {
      if (!strcmp(object->GetPointData()->GetArrayName(i),arrayname.c_str()))
      {
        exists =1;
      }
    }
  }
  else
  {
    numArrays = object->GetCellData()->GetNumberOfArrays();
    for (i=0;i<numArrays;i++)
    {
      if (!strcmp(object->GetCellData()->GetArrayName(i),arrayname.c_str()))
      {
        exists =1;
      }
    }
  }

  return exists;
}

void vtkSVBoundaryMapper::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  //os << indent << "Number of subdivisions: "
  //   << this->GetNumberOfSubdivisions() << endl;
}