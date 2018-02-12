#include "vtkSVPolyDataSurfaceInspector.h"

#include "vtkConnectivityFilter.h"
#include "vtkDataObject.h"
#include "vtkDoubleArray.h"
#include "vtkEdgeTable.h"
#include "vtkFeatureEdges.h"
#include "vtkIdList.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkSVGlobals.h"

vtkStandardNewMacro(vtkSVPolyDataSurfaceInspector);

vtkSVPolyDataSurfaceInspector::vtkSVPolyDataSurfaceInspector()
{
  this->NumberOfElements  = 0;
  this->NumberOfPoints    = 0;
  this->NumberOfEdges     = 0;
  this->NumberOfOpenEdges = 0;
  this->NumberOfNonTriangularElements = 0;
  this->NumberOfNonManifoldEdges      = 0;
  this->SurfaceGenus                  = 0;
  this->NumberOfConnectedRegions      = 0;
  this->NumberOfHoles        = 0;

  this->CheckNumberOfConnectedRegions = 0;
  this->CheckNumberOfHoles            = 0;
}

int vtkSVPolyDataSurfaceInspector::RequestData(
                                          vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
  // get the input and output
  vtkPolyData *input = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);
  output->DeepCopy(input);

  input->BuildLinks();

  int numPts = input->GetNumberOfPoints();
  this->NumberOfPoints = input->GetNumberOfPoints();

  int numPolys = input->GetNumberOfCells();
  this->NumberOfElements = input->GetNumberOfCells();

  // Start edge insertion for edge table
  vtkNew(vtkEdgeTable, surfaceEdgeTable);
  surfaceEdgeTable->InitEdgeInsertion(numPts, 1);

  this->NumberOfOpenEdges = 0;
  this->NumberOfNonTriangularElements = 0;
  this->NumberOfNonManifoldEdges = 0;
  for (int i=0; i<numPolys; i++)
  {
    vtkIdType npts, *pts;
    input->GetCellPoints(i, npts, pts);
    if (npts != 3)
    {
      this->NumberOfNonTriangularElements++;
    }
    for (int j=0; j<npts; j++)
    {
      vtkIdType p0, p1;
      p0 = pts[j];
      p1 = pts[(j+1)%npts];

      vtkNew(vtkIdList, edgeNeighbor);
      input->GetCellEdgeNeighbors(i, p0, p1, edgeNeighbor);

      if (edgeNeighbor->GetNumberOfIds() == 0)
      {
        this->NumberOfOpenEdges++;
      }
      if (edgeNeighbor->GetNumberOfIds() > 1)
      {
        this->NumberOfNonManifoldEdges++;
      }

      // Check to see if edge has already been inserted
      vtkIdType checkEdge = surfaceEdgeTable->IsEdge(p0, p1);
      if (checkEdge == -1)
      {
        // Get new edge id and insert into table
        vtkIdType edgeId = surfaceEdgeTable->InsertEdge(p0, p1);
      }
    }
  }

  int ne = surfaceEdgeTable->GetNumberOfEdges();
  int nv = numPts;
  int nf = numPolys;

  this->NumberOfEdges = ne;

  if (this->NumberOfOpenEdges > 0)
  {
    if (this->NumberOfNonManifoldEdges == 0)
    {
      vtkNew(vtkFeatureEdges, featureEdges);
      featureEdges->SetInputData(input);
      featureEdges->BoundaryEdgesOn();
      featureEdges->FeatureEdgesOff();
      featureEdges->ManifoldEdgesOff();
      featureEdges->NonManifoldEdgesOff();
      featureEdges->Update();

      vtkNew(vtkConnectivityFilter, connector);
      connector->SetInputData(featureEdges->GetOutput());
      connector->SetExtractionModeToAllRegions();
      connector->Update();

      int numHoles = connector->GetNumberOfExtractedRegions();

      int numBoundaryLines = featureEdges->GetOutput()->GetNumberOfLines();
      if (numBoundaryLines != this->NumberOfOpenEdges)
      {
        vtkWarningMacro("Feature edges and manual processing detected different number of open edges");
      }

      nv = nv + numHoles;
      nf = nf + numBoundaryLines;
      ne = ne + numBoundaryLines;

      this->SurfaceGenus = ((ne - nv - nf)/2) + 1;
    }
    else
    {
      this->SurfaceGenus = -1;
    }
  }
  else
  {
    if (this->NumberOfNonManifoldEdges == 0)
    {
      this->SurfaceGenus = ((ne - nv - nf)/2) + 1;
    }
    else
    {
      this->SurfaceGenus = -1;
    }
  }

  if (this->CheckNumberOfConnectedRegions)
  {
    vtkNew(vtkConnectivityFilter, connector);
    connector->SetInputData(input);
    connector->SetExtractionModeToAllRegions();
    connector->Update();

    this->NumberOfConnectedRegions = connector->GetNumberOfExtractedRegions();
  }

  if (this->CheckNumberOfHoles)
  {
    vtkNew(vtkFeatureEdges, featureEdges);
    featureEdges->SetInputData(input);
    featureEdges->BoundaryEdgesOn();
    featureEdges->FeatureEdgesOff();
    featureEdges->ManifoldEdgesOff();
    featureEdges->NonManifoldEdgesOff();
    featureEdges->Update();

    vtkNew(vtkConnectivityFilter, connector);
    connector->SetInputData(featureEdges->GetOutput());
    connector->SetExtractionModeToAllRegions();
    connector->Update();

    this->NumberOfHoles = connector->GetNumberOfExtractedRegions();
  }

  return SV_OK;
}
