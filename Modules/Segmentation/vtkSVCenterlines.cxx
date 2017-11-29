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

#include "vtkSVCenterlines.h"
#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkConnectivityFilter.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkDelaunay3D.h"
#include "vtkEdgeTable.h"
#include "vtkPolyDataNormals.h"
#include "vtkvmtkInternalTetrahedraExtractor.h"
#include "vtkvmtkVoronoiDiagram3D.h"
#include "vtkvmtkSimplifyVoronoiDiagram.h"
#include "vtkArrayCalculator.h"
#include "vtkvmtkNonManifoldFastMarching.h"
#include "vtkvmtkSteepestDescentLineTracer.h"
#include "vtkMath.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkTetra.h"
#include "vtkTriangleFilter.h"
#include "vtkThreshold.h"
#include "vtkPointData.h"
#include "vtkPointLocator.h"
#include "vtkIdList.h"
#include "vtkIdFilter.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkVersion.h"
#include "vtkXMLPolyDataWriter.h"


vtkStandardNewMacro(vtkSVCenterlines);

vtkCxxSetObjectMacro(vtkSVCenterlines,SourceSeedIds,vtkIdList);
vtkCxxSetObjectMacro(vtkSVCenterlines,TargetSeedIds,vtkIdList);
vtkCxxSetObjectMacro(vtkSVCenterlines,CapCenterIds,vtkIdList);

vtkSVCenterlines::vtkSVCenterlines()
{
  this->SourceSeedIds = NULL;
  this->TargetSeedIds = NULL;
  this->CapCenterIds = NULL;

  this->RadiusArrayName = NULL;
  this->CostFunction = new char[16];
  strcpy(this->CostFunction,"1/R");

  this->CostFunctionArrayName = new char[256];
  strcpy(this->CostFunctionArrayName,"CostFunctionArray");

  this->EikonalSolutionArrayName = new char[256];
  strcpy(this->EikonalSolutionArrayName,"EikonalSolutionArray");

  this->EdgeArrayName = new char[256];
  strcpy(this->EdgeArrayName,"EdgeArray");

  this->EdgePCoordArrayName = new char[256];
  strcpy(this->EdgePCoordArrayName,"EdgePCoordArray");

  this->FlipNormals = 0;
  this->SimplifyVoronoi = 0;
  this->CenterlineResampling = 0;
  this->AppendEndPointsToCenterlines = 0;

  this->ResamplingStepLength = 1.0;

  this->GenerateDelaunayTessellation = 1;

  this->DelaunayTessellation = NULL;
  this->DelaunayTolerance = 1E-3;

  this->VoronoiDiagram = vtkPolyData::New();
  this->PoleIds = vtkIdList::New();
}

vtkSVCenterlines::~vtkSVCenterlines()
{
  if (this->SourceSeedIds)
    {
    this->SourceSeedIds->Delete();
    this->SourceSeedIds = NULL;
    }

  if (this->TargetSeedIds)
    {
    this->TargetSeedIds->Delete();
    this->TargetSeedIds = NULL;
    }

  if (this->CapCenterIds)
    {
    this->CapCenterIds->Delete();
    this->CapCenterIds = NULL;
    }

  if (this->CostFunction)
    {
    delete[] this->CostFunction;
    this->CostFunction = NULL;
    }

  if (this->CostFunctionArrayName)
    {
    delete[] this->CostFunctionArrayName;
    this->CostFunctionArrayName = NULL;
    }

  if (this->EikonalSolutionArrayName)
    {
    delete[] this->EikonalSolutionArrayName;
    this->EikonalSolutionArrayName = NULL;
    }

  if (this->EdgeArrayName)
    {
    delete[] this->EdgeArrayName;
    this->EdgeArrayName = NULL;
    }

  if (this->EdgePCoordArrayName)
    {
    delete[] this->EdgePCoordArrayName;
    this->EdgePCoordArrayName = NULL;
    }

  if (this->RadiusArrayName)
    {
    delete[] this->RadiusArrayName;
    this->RadiusArrayName = NULL;
    }

  if (this->DelaunayTessellation)
    {
    this->DelaunayTessellation->Delete();
    this->DelaunayTessellation = NULL;
    }

  this->VoronoiDiagram->Delete();
  this->VoronoiDiagram = NULL;

  this->PoleIds->Delete();
  this->PoleIds = NULL;
}

int vtkSVCenterlines::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  //if (!this->SourceSeedIds)
  //  {
  //  vtkErrorMacro(<< "No SourceSeedIds set.");
  //  return 1;
  //  }

  //if (!this->TargetSeedIds)
  //  {
  //  vtkErrorMacro(<< "No TargetSeedIds set.");
  //  return 1;
  //  }

  if (!this->RadiusArrayName)
    {
    vtkErrorMacro(<< "No RadiusArrayName set.");
    return 1;
    }

  if (!this->GenerateDelaunayTessellation && !this->DelaunayTessellation)
    {
    vtkErrorMacro(<< "GenerateDelaunayTessellation is off but a DelaunayTessellation has not been set.");
    return 1;
    }

  vtkPolyDataNormals* surfaceNormals = vtkPolyDataNormals::New();
#if (VTK_MAJOR_VERSION <= 5)
  surfaceNormals->SetInput(input);
#else
  surfaceNormals->SetInputData(input);
#endif
  surfaceNormals->SplittingOff();
  surfaceNormals->AutoOrientNormalsOn();
  surfaceNormals->SetFlipNormals(this->FlipNormals);
  surfaceNormals->ComputePointNormalsOn();
  surfaceNormals->ConsistencyOn();
  surfaceNormals->Update();

  fprintf(stdout,"GENERATING DELAUNAY TESSELATION...\n");
  if (this->GenerateDelaunayTessellation)
    {
    vtkDelaunay3D* delaunayTessellator = vtkDelaunay3D::New();
    delaunayTessellator->CreateDefaultLocator();
#if (VTK_MAJOR_VERSION <= 5)
    delaunayTessellator->SetInput(surfaceNormals->GetOutput());
#else
    delaunayTessellator->SetInputConnection(surfaceNormals->GetOutputPort());
#endif
    delaunayTessellator->SetTolerance(this->DelaunayTolerance);
    delaunayTessellator->Update();

    vtkUnstructuredGrid* delaunay = delaunayTessellator->GetOutput();
    delaunay->GetPointData()->AddArray(surfaceNormals->GetOutput()->GetPointData()->GetNormals());

    vtkvmtkInternalTetrahedraExtractor* internalTetrahedraExtractor = vtkvmtkInternalTetrahedraExtractor::New();
#if (VTK_MAJOR_VERSION <= 5)
    internalTetrahedraExtractor->SetInput(delaunayTessellator->GetOutput());
#else
    internalTetrahedraExtractor->SetInputConnection(delaunayTessellator->GetOutputPort());
#endif
    internalTetrahedraExtractor->SetOutwardNormalsArrayName(surfaceNormals->GetOutput()->GetPointData()->GetNormals()->GetName());
    if (this->CapCenterIds)
      {
      internalTetrahedraExtractor->UseCapsOn();
      internalTetrahedraExtractor->SetCapCenterIds(this->CapCenterIds);
      }
    internalTetrahedraExtractor->Update();

    this->DelaunayTessellation = internalTetrahedraExtractor->GetOutput();
    this->DelaunayTessellation->Register(this);

    delaunayTessellator->Delete();
    internalTetrahedraExtractor->Delete();
    }

  fprintf(stdout,"GENERATING VORONOI DIAGRAM...\n");
  vtkvmtkVoronoiDiagram3D* voronoiDiagramFilter = vtkvmtkVoronoiDiagram3D::New();
#if (VTK_MAJOR_VERSION <= 5)
  voronoiDiagramFilter->SetInput(this->DelaunayTessellation);
#else
  voronoiDiagramFilter->SetInputData(this->DelaunayTessellation);
#endif
  voronoiDiagramFilter->SetRadiusArrayName(this->RadiusArrayName);
  voronoiDiagramFilter->Update();

  this->PoleIds->DeepCopy(voronoiDiagramFilter->GetPoleIds());

  vtkPolyData* voronoiDiagram = voronoiDiagramFilter->GetOutput();

  if (this->SimplifyVoronoi)
    {
    vtkvmtkSimplifyVoronoiDiagram* voronoiDiagramSimplifier = vtkvmtkSimplifyVoronoiDiagram::New();
#if (VTK_MAJOR_VERSION <= 5)
    voronoiDiagramSimplifier->SetInput(voronoiDiagramFilter->GetOutput());
#else
    voronoiDiagramSimplifier->SetInputConnection(voronoiDiagramFilter->GetOutputPort());
#endif
    voronoiDiagramSimplifier->SetUnremovablePointIds(voronoiDiagramFilter->GetPoleIds());
    voronoiDiagramSimplifier->Update();
    voronoiDiagram = voronoiDiagramSimplifier->GetOutput();
    voronoiDiagram->Register(this);
    voronoiDiagramSimplifier->Delete();
    }

  fprintf(stdout,"PRUNING VORONOI DIAGRAM...\n");
  vtkTriangleFilter *triangulator = vtkTriangleFilter::New();
  triangulator->SetInputData(voronoiDiagram);
  triangulator->Update();

  //vtkCleanPolyData *cleaner0 = vtkCleanPolyData::New();
  //cleaner0->SetInputData(triangulator->GetOutput());
  //cleaner0->ToleranceIsAbsoluteOn();
  //cleaner0->SetAbsoluteTolerance(0.001);
  //cleaner0->Update();

  vtkPolyData *triPd = vtkPolyData::New();
  triPd->DeepCopy(triangulator->GetOutput());
  triPd->BuildLinks();
  triangulator->Delete();

  int numCells = triPd->GetNumberOfCells();
  int numPts = triPd->GetNumberOfPoints();

  vtkIntArray *tmpCellArray = vtkIntArray::New();
  tmpCellArray->SetNumberOfTuples(numCells);
  tmpCellArray->SetName("TmpInternalIds");
  for (int i=0; i<numCells; i++)
    tmpCellArray->SetTuple1(i, i);
  triPd->GetCellData()->AddArray(tmpCellArray);
  tmpCellArray->Delete();

  //Make edge table
  // ------------------------------------------------------------------------

  // Start edge insertion for edge table
  vtkEdgeTable *edgeTable = vtkEdgeTable::New();
  vtkIntArray *edgeNeighbors = vtkIntArray::New();
  edgeNeighbors->SetNumberOfComponents(2);
  edgeTable->InitEdgeInsertion(numPts, 1);

  // Loop through cells
  int totEdges = 0;
  for (int i=0; i<numCells; i++)
  {
    // Get cellpoints
    vtkIdType npts, *pts;
    triPd->GetCellPoints(i, npts, pts);
    for (int j=0; j<npts; j++)
    {
      // Get each edge of cell
      vtkIdType p0 = pts[j];
      vtkIdType p1 = pts[(j+1)%npts];

      vtkIdList *neighborCellIds = vtkIdList::New();
      triPd->GetCellEdgeNeighbors(i, p0, p1, neighborCellIds);
      vtkIdType neighborCellId = 0;

      // Check to see if it is a boundary edge
      if (neighborCellIds->GetNumberOfIds() > 0)
        neighborCellId = neighborCellIds->GetId(0);
      else
      {
        neighborCellId = -1;
      }

      // Check to see if edge has already been inserted
      vtkIdType checkEdge = edgeTable->IsEdge(p0, p1);
      if (checkEdge == -1)
      {
        totEdges++;
        // Get new edge id and insert into table
        vtkIdType edgeId = edgeTable->InsertEdge(p0, p1);

        // Insert edge weights and neighboring cells`
        edgeNeighbors->InsertComponent(edgeId, 0, i);
        edgeNeighbors->InsertComponent(edgeId, 1, neighborCellId);
      }
      neighborCellIds->Delete();
    }
  }

  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Now make polydata from edge table

  vtkCellArray *edgeCells = vtkCellArray::New();
  edgeTable->InitTraversal();
  int edgeId = 0;
  for (edgeId = 0; edgeId < totEdges;  edgeId++)
  {
    vtkIdType edgePtId0, edgePtId1;
    edgeTable->GetNextEdge(edgePtId0, edgePtId1);

    vtkIdList *newEdgeCell = vtkIdList::New();
    newEdgeCell->SetNumberOfIds(2);
    newEdgeCell->SetId(0, edgePtId0);
    newEdgeCell->SetId(1, edgePtId1);

    edgeCells->InsertNextCell(newEdgeCell);

    newEdgeCell->Delete();
  }
  // ------------------------------------------------------------------------

  vtkPolyData *edgePd = vtkPolyData::New();
  edgePd->SetLines(edgeCells);
  edgePd->SetPoints(triPd->GetPoints());
  edgePd->BuildLinks();

  vtkPolyData *newEdgePd = vtkPolyData::New();
  vtkPolyData *newTriPd = vtkPolyData::New();

  this->PruneVoronoiDiagram(triPd, edgePd, newTriPd, newEdgePd, "");
  fprintf(stdout,"DONE COMPUTING REMOVAL ITERATIONS\n");

  std::string fullfn = "/Users/adamupdegrove/Desktop/tmp/VORONOI_TRI_REMOVAL.vtp";
  vtkXMLPolyDataWriter *newWriter = vtkXMLPolyDataWriter::New();
  newWriter->SetInputData(newTriPd);
  newWriter->SetFileName(fullfn.c_str());
  newWriter->Write();
  newWriter->Delete();

  std::string edgefn = "/Users/adamupdegrove/Desktop/tmp/VORONOI_EDGE_REMOVAL.vtp";
  vtkXMLPolyDataWriter *myWriter = vtkXMLPolyDataWriter::New();
  myWriter->SetInputData(newEdgePd);
  myWriter->SetFileName(edgefn.c_str());
  myWriter->Write();
  myWriter->Delete();

  //vtkPolyData *thresholdCellPd = vtkPolyData::New();
  //vtkPolyData *thresholdEdgePd = vtkPolyData::New();
  //thresholdCellPd->DeepCopy(newTriPd);
  //thresholdEdgePd->DeepCopy(newEdgePd);

  edgeCells->Delete();
  edgeNeighbors->Delete();
  edgeTable->Delete();

  vtkIntArray *tmpEdgeArray = vtkIntArray::New();
  tmpEdgeArray->SetNumberOfTuples(newEdgePd->GetNumberOfCells());
  tmpEdgeArray->SetName("TmpInternalIds");
  for (int i=0; i<newEdgePd->GetNumberOfCells(); i++)
    tmpEdgeArray->SetTuple1(i, i);
  newEdgePd->GetCellData()->AddArray(tmpEdgeArray);
  tmpEdgeArray->Delete();

  int mAbsThr = 3;
  double mAbsRange[2];
  newEdgePd->GetCellData()->GetArray("MAbs")->GetRange(mAbsRange);
  vtkThreshold *mAbsThresholder = vtkThreshold::New();
  mAbsThresholder->SetInputData(newEdgePd);
  mAbsThresholder->SetInputArrayToProcess(0, 0, 0, 1, "MAbs");
  mAbsThresholder->ThresholdBetween(mAbsThr, mAbsRange[1]);
  mAbsThresholder->Update();
  fprintf(stdout,"Thresholded MAbs: %d\n", mAbsThresholder->GetOutput()->GetNumberOfCells());

  double mRelThr = 0.5;
  double mRelRange[2];
  newEdgePd->GetCellData()->GetArray("MRel")->GetRange(mRelRange);
  vtkThreshold *mRelThresholder = vtkThreshold::New();
  mRelThresholder->SetInputData(mAbsThresholder->GetOutput());
  mRelThresholder->SetInputArrayToProcess(0, 0, 0, 1, "MRel");
  mRelThresholder->ThresholdBetween(mRelThr, mRelRange[1]);
  mRelThresholder->Update();
  fprintf(stdout,"Thresholded MRel: %d\n", mRelThresholder->GetOutput()->GetNumberOfCells());

  vtkConnectivityFilter *connector = vtkConnectivityFilter::New();
  connector->SetInputData(mRelThresholder->GetOutput());
  connector->SetExtractionModeToAllRegions();
  connector->ColorRegionsOn();
  connector->Update();

  vtkDataSetSurfaceFilter *surfacer = vtkDataSetSurfaceFilter::New();
  surfacer->SetInputData(connector->GetOutput());
  surfacer->Update();

  vtkPolyData *leftOver = vtkPolyData::New();
  leftOver->DeepCopy(surfacer->GetOutput());
  surfacer->Delete();

  mAbsThresholder->Delete();
  mRelThresholder->Delete();

  vtkIntArray *leaveArray = vtkIntArray::New();
  leaveArray->SetNumberOfTuples(edgePd->GetNumberOfCells());
  leaveArray->SetName("MedialEdges");
  for (int i=0; i<edgePd->GetNumberOfCells(); i++)
    leaveArray->SetTuple1(i, 0);

  int connectThr = 5;
  for (int i=0; i<connector->GetNumberOfExtractedRegions(); i++)
  {
    vtkThreshold *regionThresholder = vtkThreshold::New();
    regionThresholder->SetInputData(leftOver);
    regionThresholder->SetInputArrayToProcess(0, 0, 0, 1, "RegionId");
    regionThresholder->ThresholdBetween(i, i);
    regionThresholder->Update();

    fprintf(stdout,"Thresholded Region %d: %d\n", i, regionThresholder->GetOutput()->GetNumberOfCells());
    if (regionThresholder->GetOutput()->GetNumberOfCells() > connectThr)
    {
      for (int j=0; j<regionThresholder->GetOutput()->GetNumberOfCells(); j++)
      {
        int origCellId = regionThresholder->GetOutput()->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(j);
        leaveArray->SetTuple1(origCellId, 1);
      }
    }
  }

  connector->Delete();
  edgePd->GetCellData()->AddArray(leaveArray);
  leaveArray->Delete();

  vtkPolyData *nextTriPd = vtkPolyData::New();
  vtkPolyData *nextEdgePd = vtkPolyData::New();

  // --------------------------------------------------------------
  // Now prune again

  fprintf(stdout,"THINNING AGAIN\n");
  this->PruneVoronoiDiagram(triPd, edgePd, nextTriPd, nextEdgePd, "MedialEdges");

  std::string finedgefn = "/Users/adamupdegrove/Desktop/tmp/VORONOI_FINAL_EDGE.vtp";
  vtkXMLPolyDataWriter *nextWriter = vtkXMLPolyDataWriter::New();
  nextWriter->SetInputData(nextEdgePd);
  nextWriter->SetFileName(finedgefn.c_str());
  nextWriter->Write();
  nextWriter->Delete();

  double finalRange[2];
  nextEdgePd->GetCellData()->GetArray("RemovalIteration")->GetRange(finalRange);

  vtkIntArray *keepCellArray = vtkIntArray::New();
  keepCellArray->SetNumberOfTuples(triPd->GetNumberOfCells());
  keepCellArray->SetName("KeepCellArray");
  for (int i=0; i<triPd->GetNumberOfCells(); i++)
    keepCellArray->SetTuple1(i, 0);
  for (int i=0; i<nextEdgePd->GetNumberOfCells(); i++)
  {
    int rVal = nextEdgePd->GetCellData()->GetArray("RemovalIteration")->GetTuple1(i);

    if (rVal == finalRange[1])
    {
      vtkIdType npts, *pts;
      nextEdgePd->GetCellPoints(i, npts, pts);

      for (int j=0; j<npts; j++)
      {
        vtkIdList *pointCellIds = vtkIdList::New();
        triPd->GetPointCells(pts[j], pointCellIds);

        for (int k=0; k<pointCellIds->GetNumberOfIds(); k++)
        {
          int cellId = pointCellIds->GetId(k);
          keepCellArray->SetTuple1(pointCellIds->GetId(k), 1);
        }

        pointCellIds->Delete();
      }

      //vtkIdList *pointIds = vtkIdList::New();
      //pointIds->SetNumberOfIds(2);
      //pointIds->SetId(0, pts[0]);
      //pointIds->SetId(1, pts[1]);

      //vtkIdList *attachedCellIds = vtkIdList::New();

      //triPd->GetCellNeighbors(-1, pointIds, attachedCellIds);

      //for (int j=0; j<attachedCellIds->GetNumberOfIds(); j++)
      //  keepCellArray->SetTuple1(attachedCellIds->GetId(j), 1);

      //pointIds->Delete();
      //attachedCellIds->Delete();
    }
  }

  triPd->GetCellData()->AddArray(keepCellArray);
  keepCellArray->Delete();

  vtkIntArray *tmpPtArray = vtkIntArray::New();
  tmpPtArray->SetNumberOfTuples(triPd->GetNumberOfPoints());
  tmpPtArray->SetName("TmpInternalIds");
  for (int i=0; i<triPd->GetNumberOfPoints(); i++)
    tmpPtArray->SetTuple1(i, i);
  triPd->GetPointData()->AddArray(tmpPtArray);
  tmpPtArray->Delete();

  vtkThreshold *keepThreshold = vtkThreshold::New();
  keepThreshold->SetInputData(triPd);
  keepThreshold->SetInputArrayToProcess(0, 0, 0, 1, "KeepCellArray");
  keepThreshold->ThresholdBetween(1, 1);
  keepThreshold->Update();

  vtkDataSetSurfaceFilter *surfacer2 = vtkDataSetSurfaceFilter::New();
  surfacer2->SetInputData(keepThreshold->GetOutput());
  surfacer2->Update();

  vtkCleanPolyData *cleaner = vtkCleanPolyData::New();
  cleaner->SetInputData(surfacer2->GetOutput());
  cleaner->Update();

  vtkPolyData *reducedVorPd = vtkPolyData::New();
  reducedVorPd->DeepCopy(cleaner->GetOutput());
  surfacer2->Delete();
  cleaner->Delete();

  ////
  //// --------------------------------------------------------------
  //// new seed ids
  //vtkIdList *newSourceSeeds = vtkIdList::New();
  //vtkIdList *newTargetSeeds = vtkIdList::New();

  //newTargetSeeds->SetNumberOfIds(voronoiTargetSeedIds->GetNumberOfIds());
  //newSourceSeeds->SetNumberOfIds(voronoiSourceSeedIds->GetNumberOfIds());

  //vtkPointLocator *locator = vtkPointLocator::New();
  //locator->SetDataSet(reducedVorPd);
  ////locator->SetDataSet(triPd);
  //locator->BuildLocator();

  //for (int i=0; i<voronoiTargetSeedIds->GetNumberOfIds(); i++)
  //{
  //  int vorPtId = voronoiTargetSeedIds->GetId(i);
  //  double pt[3];
  //  voronoiFastMarching->GetOutput()->GetPoint(vorPtId, pt);

  //  int ptId = locator->FindClosestPoint(pt);

  //  newTargetSeeds->SetId(i, ptId);
  //}

  //for (int i=0; i<voronoiSourceSeedIds->GetNumberOfIds(); i++)
  //{
  //  int vorPtId = voronoiSourceSeedIds->GetId(i);
  //  double pt[3];
  //  voronoiFastMarching->GetOutput()->GetPoint(vorPtId, pt);

  //  int ptId = locator->FindClosestPoint(pt);

  //  newSourceSeeds->SetId(i, ptId);
  //}
  //newSourceSeeds->Delete();
  //newTargetSeeds->Delete();
  // --------------------------------------------------------------
  //
  vtkPolyData *newVoronoiDiagram = vtkPolyData::New();
  newVoronoiDiagram->DeepCopy(voronoiDiagram);

  vtkDoubleArray *newCostFunctionArray = vtkDoubleArray::New();
  newCostFunctionArray->SetNumberOfTuples(newVoronoiDiagram->GetNumberOfPoints());
  newCostFunctionArray->SetName("NewCostArray");

  //for (int i=0; i<newVoronoiDiagram->GetNumberOfPoints(); i++)
  //{
  //  double radVal = newVoronoiDiagram->GetPointData()->GetArray(this->RadiusArrayName)->GetTuple1(i);
  //  newCostFunctionArray->SetTuple1(i, radVal);
  //}

  //reducedVorPd->BuildLinks();
  //for (int i=0; i<reducedVorPd->GetNumberOfPoints(); i++)
  //{
  //  vtkIdList *ptCellIds = vtkIdList::New();
  //  reducedVorPd->GetPointCells(i, ptCellIds);

  //  if (ptCellIds->GetNumberOfIds() > 2)
  //  {
  //    int origPtId = reducedVorPd->GetPointData()->GetArray("TmpInternalIds")->GetTuple1(i);
  //    newCostFunctionArray->SetTuple1(origPtId, 1000000000.0);
  //  }
  //  ptCellIds->Delete();
  //}

  for (int i=0; i<newVoronoiDiagram->GetNumberOfPoints(); i++)
  {
    //vtkIdList *pointCellIds = vtkIdList::New();
    //nextTriPd->GetPointCells(i, pointCellIds);

    //double mAbsVal = 0.0;
    //for (int j=0; j<pointCellIds->GetNumberOfIds(); j++)
    //{
    //  int cellId = pointCellIds->GetId(j);
    //  double rVal = nextTriPd->GetCellData()->GetArray("RemovalIteration")->GetTuple1(cellId);
    //  mAbsVal += rVal;
    //}

    //mAbsVal = mAbsVal / pointCellIds->GetNumberOfIds();
    //if (mAbsVal == 0.0)
    //  mAbsVal = 0.00000001;

    newCostFunctionArray->SetTuple1(i, 0.000000001);

    //pointCellIds->Delete();
  }

  for (int i=0; i<nextEdgePd->GetNumberOfCells(); i++)
  {
    int rVal = nextEdgePd->GetCellData()->GetArray("RemovalIteration")->GetTuple1(i);

    if (rVal == finalRange[1])
    {
      vtkIdType npts, *pts;
      nextEdgePd->GetCellPoints(i, npts, pts);

      int mVal = nextEdgePd->GetCellData()->GetArray("MAbs")->GetTuple1(i);
      if (mVal == 0)
        mVal = 0.000000001;

      for (int j=0; j<npts; j++)
        newCostFunctionArray->SetTuple1(pts[j], mVal);
    }
  }

  newVoronoiDiagram->GetPointData()->AddArray(newCostFunctionArray);
  newCostFunctionArray->Delete();

  nextEdgePd->GetPointData()->AddArray(triPd->GetPointData()->GetArray(this->RadiusArrayName));

  std::string newfn = "/Users/adamupdegrove/Desktop/tmp/EDGES_WITH_RADIUS.vtp";
  vtkXMLPolyDataWriter *heWriter = vtkXMLPolyDataWriter::New();
  heWriter->SetInputData(nextEdgePd);
  heWriter->SetFileName(newfn.c_str());
  heWriter->Write();
  heWriter->Delete();
  reducedVorPd->Delete();

  double radRange[2];
  voronoiDiagram->GetPointData()->GetArray(this->RadiusArrayName)->GetRange(radRange);
  fprintf(stdout,"RADIUS ARRAY RANGE: %.6f %.6f\n", radRange[0], radRange[1]);
  fprintf(stdout,"1/RADIUS ARRAY RANGE: %.6f %.6f\n", 1./radRange[0], 1./radRange[1]);

  // --------------------------------------------------------------
  edgePd->Delete();
  triPd->Delete();
  newTriPd->Delete();
  newEdgePd->Delete();
  nextTriPd->Delete();
  nextEdgePd->Delete();

//  vtkArrayCalculator* voronoiCostFunctionCalculator = vtkArrayCalculator::New();
////#if (VTK_MAJOR_VERSION <= 5)
////  voronoiCostFunctionCalculator->SetInput(voronoiDiagram);
////#else
////  voronoiCostFunctionCalculator->SetInputData(voronoiDiagram);
////#endif
//  voronoiCostFunctionCalculator->SetInputData(newVoronoiDiagram);
//  voronoiCostFunctionCalculator->SetAttributeModeToUsePointData();
//  //voronoiCostFunctionCalculator->AddScalarVariable("R",this->RadiusArrayName,0);
//  voronoiCostFunctionCalculator->AddScalarVariable("R","NewCostArray",0);
//  voronoiCostFunctionCalculator->SetFunction(this->CostFunction);
//  voronoiCostFunctionCalculator->SetResultArrayName(this->CostFunctionArrayName);
//  voronoiCostFunctionCalculator->Update();
//
//  vtkIdList* voronoiSourceSeedIds = vtkIdList::New();
//  vtkIdList* voronoiTargetSeedIds = vtkIdList::New();
//
//  vtkIdList* voronoiSeeds = vtkIdList::New();
//
//  int i;
//  if (this->CapCenterIds)
//    {
//    this->FindVoronoiSeeds(this->DelaunayTessellation,this->CapCenterIds,surfaceNormals->GetOutput()->GetPointData()->GetNormals(),voronoiSeeds);
//    for (i=0; i<this->SourceSeedIds->GetNumberOfIds(); i++)
//      {
//      voronoiSourceSeedIds->InsertNextId(voronoiSeeds->GetId(this->SourceSeedIds->GetId(i)));
//      }
//    for (i=0; i<this->TargetSeedIds->GetNumberOfIds(); i++)
//      {
//      voronoiTargetSeedIds->InsertNextId(voronoiSeeds->GetId(this->TargetSeedIds->GetId(i)));
//      }
//    }
//  else
//    {
//    for (i=0; i<this->SourceSeedIds->GetNumberOfIds(); i++)
//      {
//      voronoiSourceSeedIds->InsertNextId(this->PoleIds->GetId(this->SourceSeedIds->GetId(i)));
//      }
//    for (i=0; i<this->TargetSeedIds->GetNumberOfIds(); i++)
//      {
//      voronoiTargetSeedIds->InsertNextId(this->PoleIds->GetId(this->TargetSeedIds->GetId(i)));
//      }
//    }
//
//  vtkvmtkNonManifoldFastMarching* voronoiFastMarching = vtkvmtkNonManifoldFastMarching::New();
//#if (VTK_MAJOR_VERSION <= 5)
//  voronoiFastMarching->SetInput(vtkPolyData::SafeDownCast(voronoiCostFunctionCalculator->GetOutput()));
//#else
//  voronoiFastMarching->SetInputConnection(voronoiCostFunctionCalculator->GetOutputPort());
//#endif
//  voronoiFastMarching->SetCostFunctionArrayName(this->CostFunctionArrayName);
//  voronoiFastMarching->SetSolutionArrayName(this->EikonalSolutionArrayName);
//  voronoiFastMarching->SeedsBoundaryConditionsOn();
//  voronoiFastMarching->SetSeeds(voronoiSourceSeedIds);
//  voronoiFastMarching->Update();
//
//  this->VoronoiDiagram->ShallowCopy(voronoiFastMarching->GetOutput());
//#if (VTK_MAJOR_VERSION <= 5)
//  this->VoronoiDiagram->Update();
//#endif
//
//
//  vtkvmtkSteepestDescentLineTracer* centerlineBacktracing = vtkvmtkSteepestDescentLineTracer::New();
//#if (VTK_MAJOR_VERSION <= 5)
//  centerlineBacktracing->SetInput(voronoiFastMarching->GetOutput());
//#else
//  centerlineBacktracing->SetInputConnection(voronoiFastMarching->GetOutputPort());
//#endif
//  //centerlineBacktracing->SetInputData(reducedVorPd);
//  centerlineBacktracing->SetDataArrayName(this->RadiusArrayName);
//  centerlineBacktracing->SetDescentArrayName(this->EikonalSolutionArrayName);
//  //centerlineBacktracing->SetDescentArrayName("MAbs");
//  centerlineBacktracing->SetEdgeArrayName(this->EdgeArrayName);
//  centerlineBacktracing->SetEdgePCoordArrayName(this->EdgePCoordArrayName);
//  centerlineBacktracing->SetSeeds(voronoiTargetSeedIds);
//  //centerlineBacktracing->SetSeeds(newTargetSeeds);
//  centerlineBacktracing->MergePathsOff();
//  centerlineBacktracing->StopOnTargetsOn();
//  centerlineBacktracing->SetTargets(voronoiSourceSeedIds);
//  //centerlineBacktracing->SetTargets(newSourceSeeds);
//  centerlineBacktracing->Update();
//
//  output->ShallowCopy(centerlineBacktracing->GetOutput());
//
//  vtkIdList* hitTargets = centerlineBacktracing->GetHitTargets();
//
//  vtkPoints* endPointPairs = vtkPoints::New();
//
//  const vtkIdType numTargetSeedIds = this->TargetSeedIds->GetNumberOfIds();
//  const vtkIdType numHitTargets = hitTargets->GetNumberOfIds();
//  if(numHitTargets == numTargetSeedIds) {
//  if (this->AppendEndPointsToCenterlines)
//    {
//    for (i=0; i<numTargetSeedIds; i++)
//      {
//      if (this->CapCenterIds)
//        {
//        vtkIdType endPointId1 = this->CapCenterIds->GetId(this->TargetSeedIds->GetId(i));
//        vtkIdType hitTargetPointId = hitTargets->GetId(i);
//        vtkIdType targetId = voronoiSourceSeedIds->IsId(hitTargetPointId);
//        vtkIdType endPointId2 = this->CapCenterIds->GetId(this->SourceSeedIds->GetId(targetId));
//        endPointPairs->InsertNextPoint(input->GetPoint(endPointId1));
//        endPointPairs->InsertNextPoint(input->GetPoint(endPointId2));
//        }
//      else
//        {
//        vtkIdType endPointId1 = this->TargetSeedIds->GetId(i);
//        vtkIdType hitTargetPointId = hitTargets->GetId(i);
//        vtkIdType targetId = voronoiSourceSeedIds->IsId(hitTargetPointId);
//        vtkIdType endPointId2 = this->SourceSeedIds->GetId(targetId);
//        endPointPairs->InsertNextPoint(input->GetPoint(endPointId1));
//        endPointPairs->InsertNextPoint(input->GetPoint(endPointId2));
//        }
//      }
//
//    this->AppendEndPoints(endPointPairs);
//    }
//  }
//
//  std::string thefn = "/Users/adamupdegrove/Desktop/tmp/REG_CENTERLINES.vtp";
//  vtkXMLPolyDataWriter *duhWriter = vtkXMLPolyDataWriter::New();
//  duhWriter->SetInputData(output);
//  duhWriter->SetFileName(thefn.c_str());
//  duhWriter->Write();
//  duhWriter->Delete();
//
//  if (this->CenterlineResampling)
//    {
//    this->ResampleCenterlines();
//    }
//  this->ReverseCenterlines();
//
//  surfaceNormals->Delete();
//  voronoiDiagramFilter->Delete();
//  voronoiCostFunctionCalculator->Delete();
//  voronoiSeeds->Delete();
//  voronoiSourceSeedIds->Delete();
//  voronoiTargetSeedIds->Delete();
//  voronoiFastMarching->Delete();
//  centerlineBacktracing->Delete();
//  endPointPairs->Delete();

  return 1;
}

void vtkSVCenterlines::FindVoronoiSeeds(vtkUnstructuredGrid *delaunay, vtkIdList *boundaryBaricenterIds, vtkDataArray *normals, vtkIdList *seedIds)
{
  vtkIdType i, j;
  vtkIdList *pointCells;
  vtkIdType baricenterId;
  double baricenter[3], normal[3];
  double maxRadius, secondMaxRadius;
  vtkTetra* tetra;
  double p0[3], p1[3], p2[3], p3[3];
  double circumcenter[3], circumradius, tetraRadius;
  double referenceVector[3];
  double pole[3], poleVector[3], secondPole[3], secondPoleVector[3];
  pole[0] = pole[1] = pole[2] = 0.0;
  vtkIdType maxRadiusCellId, secondMaxRadiusCellId;

  pointCells = vtkIdList::New();

  for (i=0; i<boundaryBaricenterIds->GetNumberOfIds(); i++)
    {
    baricenterId = boundaryBaricenterIds->GetId(i);
    delaunay->GetPoint(baricenterId,baricenter);
    normals->GetTuple(baricenterId,normal);
    pointCells->Initialize();
    delaunay->GetPointCells(baricenterId,pointCells);
    maxRadius = 0.0;
    maxRadiusCellId = -1;
    secondMaxRadiusCellId = -1;

    for (j=0; j<pointCells->GetNumberOfIds(); j++)
      {
      tetra = vtkTetra::SafeDownCast(delaunay->GetCell(pointCells->GetId(j)));
      tetra->GetPoints()->GetPoint(0,p0);
      tetra->GetPoints()->GetPoint(1,p1);
      tetra->GetPoints()->GetPoint(2,p2);
      tetra->GetPoints()->GetPoint(3,p3);

      circumradius = vtkTetra::Circumsphere(p0,p1,p2,p3,circumcenter);
      tetraRadius = sqrt(circumradius);

      if (tetraRadius - maxRadius > VTK_VMTK_DOUBLE_TOL)
        {
        maxRadius = tetraRadius;
        maxRadiusCellId = pointCells->GetId(j);
        pole[0] = circumcenter[0];
        pole[1] = circumcenter[1];
        pole[2] = circumcenter[2];
        }
      }

    poleVector[0] = pole[0] - baricenter[0];
    poleVector[1] = pole[1] - baricenter[1];
    poleVector[2] = pole[2] - baricenter[2];

    secondMaxRadius = 0.0;

    for (j=0; j<pointCells->GetNumberOfIds(); j++)
      {
      tetra = vtkTetra::SafeDownCast(delaunay->GetCell(pointCells->GetId(j)));
      tetra->GetPoints()->GetPoint(0,p0);
      tetra->GetPoints()->GetPoint(1,p1);
      tetra->GetPoints()->GetPoint(2,p2);
      tetra->GetPoints()->GetPoint(3,p3);

      circumradius = vtkTetra::Circumsphere(p0,p1,p2,p3,circumcenter);
      tetraRadius = sqrt(circumradius);

      referenceVector[0] = circumcenter[0] - baricenter[0];
      referenceVector[1] = circumcenter[1] - baricenter[1];
      referenceVector[2] = circumcenter[2] - baricenter[2];

      if ((tetraRadius - secondMaxRadius > VTK_VMTK_DOUBLE_TOL) && (vtkMath::Dot(poleVector,referenceVector) < VTK_VMTK_DOUBLE_TOL))
        {
        secondMaxRadius = tetraRadius;
        secondMaxRadiusCellId = pointCells->GetId(j);
        secondPole[0] = circumcenter[0];
        secondPole[1] = circumcenter[1];
        secondPole[2] = circumcenter[2];
        }
      }

    secondPoleVector[0] = secondPole[0] - baricenter[0];
    secondPoleVector[1] = secondPole[1] - baricenter[1];
    secondPoleVector[2] = secondPole[2] - baricenter[2];

    if (vtkMath::Dot(poleVector,normal) < VTK_VMTK_DOUBLE_TOL)
      {
      seedIds->InsertNextId(maxRadiusCellId);
      }
    else
      {
      seedIds->InsertNextId(secondMaxRadiusCellId);
      }
    }

  pointCells->Delete();
}

void vtkSVCenterlines::AppendEndPoints(vtkPoints* endPointPairs)
{
  vtkIdType endPointId1, endPointId2;
  vtkPolyData* output = this->GetOutput();
  vtkPolyData* completeCenterlines = vtkPolyData::New();
  vtkPoints* completeCenterlinesPoints = vtkPoints::New();
  vtkCellArray* completeCenterlinesCellArray = vtkCellArray::New();
  vtkDoubleArray* completeCenterlinesRadiusArray = vtkDoubleArray::New();
  completeCenterlinesRadiusArray->SetName(this->RadiusArrayName);
  vtkIdList* completeCell = vtkIdList::New();

  vtkDoubleArray* centerlinesRadiusArray = vtkDoubleArray::SafeDownCast(output->GetPointData()->GetArray(this->RadiusArrayName));

  completeCenterlinesPoints->DeepCopy(output->GetPoints());
  completeCenterlinesRadiusArray->DeepCopy(centerlinesRadiusArray);

  for (int k=0; k<output->GetNumberOfCells(); k++)
    {
    vtkCell* cell = output->GetCell(k);

    endPointId1 = completeCenterlinesPoints->InsertNextPoint(endPointPairs->GetPoint(2*k));
    endPointId2 = completeCenterlinesPoints->InsertNextPoint(endPointPairs->GetPoint(2*k+1));

    completeCell->Initialize();
    completeCell->SetNumberOfIds(cell->GetNumberOfPoints()+2);

    completeCell->SetId(0,endPointId1);

    for (int i=0; i<cell->GetNumberOfPoints(); i++)
      {
      completeCell->SetId(i+1,cell->GetPointId(i));
      }
    completeCell->SetId(cell->GetNumberOfPoints()+1,endPointId2);

    completeCenterlinesCellArray->InsertNextCell(completeCell);

    completeCenterlinesRadiusArray->InsertNextValue(centerlinesRadiusArray->GetValue(cell->GetPointId(0)));
    completeCenterlinesRadiusArray->InsertNextValue(centerlinesRadiusArray->GetValue(cell->GetPointId(cell->GetNumberOfPoints()-1)));
    }

  completeCenterlines->SetPoints(completeCenterlinesPoints);
  completeCenterlines->SetLines(completeCenterlinesCellArray);
  completeCenterlines->GetPointData()->AddArray(completeCenterlinesRadiusArray);
#if (VTK_MAJOR_VERSION <= 5)
  completeCenterlines->Update();
#endif

  output->ShallowCopy(completeCenterlines);

  completeCell->Delete();
  completeCenterlines->Delete();
  completeCenterlinesPoints->Delete();
  completeCenterlinesCellArray->Delete();
  completeCenterlinesRadiusArray->Delete();
}

void vtkSVCenterlines::ResampleCenterlines()
{
  vtkPolyData* output = this->GetOutput();
  vtkPolyData* resampledCenterlines = vtkPolyData::New();
  vtkPoints* resampledCenterlinesPoints = vtkPoints::New();
  vtkCellArray* resampledCenterlinesCellArray = vtkCellArray::New();
  vtkDoubleArray* resampledCenterlinesRadiusArray = vtkDoubleArray::New();
  resampledCenterlinesRadiusArray->SetName(this->RadiusArrayName);
  vtkIdList* resampledCell = vtkIdList::New();

  vtkDoubleArray* centerlinesRadiusArray = vtkDoubleArray::SafeDownCast(output->GetPointData()->GetArray(this->RadiusArrayName));

  for (int k=0; k<output->GetNumberOfCells(); k++)
    {
    vtkCell* cell = output->GetCell(k);

    resampledCell->Initialize();

    vtkIdType id = resampledCenterlinesPoints->InsertNextPoint(cell->GetPoints()->GetPoint(0));
    resampledCell->InsertNextId(id);
    resampledCenterlinesRadiusArray->InsertNextValue(centerlinesRadiusArray->GetValue(cell->GetPointId(0)));

    double point0[3], point1[3], point[3];
    double abscissa, lineAbscissa, lineLength, stepAbscissa;

    abscissa = 0.0;
    lineAbscissa = 0.0;
    lineLength = 0.0;
    stepAbscissa = 0.0;

    for (int i=0; i<cell->GetNumberOfPoints()-1; i++)
      {
      cell->GetPoints()->GetPoint(i,point0);
      cell->GetPoints()->GetPoint(i+1,point1);

      double scalar0 = centerlinesRadiusArray->GetValue(cell->GetPointId(i));
      double scalar1 = centerlinesRadiusArray->GetValue(cell->GetPointId(i+1));

      double length = sqrt(vtkMath::Distance2BetweenPoints(point0,point1));

      if (length < this->ResamplingStepLength - stepAbscissa)
        {
        stepAbscissa = stepAbscissa + length;
        continue;
        }

      double pcoord = 0.0;
      double pcoordStep = this->ResamplingStepLength / length;
      while (pcoord < 1.0)
        {
        point[0] = point0[0] + (point1[0] - point0[0]) * pcoord;
        point[1] = point0[1] + (point1[1] - point0[1]) * pcoord;
        point[2] = point0[2] + (point1[2] - point0[2]) * pcoord;

        double scalar = scalar0 + (scalar1 - scalar0) * pcoord;

        vtkIdType id = resampledCenterlinesPoints->InsertNextPoint(point);
        resampledCell->InsertNextId(id);
        resampledCenterlinesRadiusArray->InsertNextValue(scalar);

        if (pcoord + pcoordStep > 1.0)
          {
          break;
          }
        pcoord = pcoord + pcoordStep;
        }
      stepAbscissa = (1.0 - pcoord) * length;
      }

    id = resampledCenterlinesPoints->InsertNextPoint(cell->GetPoints()->GetPoint(cell->GetNumberOfPoints()-1));
    resampledCell->InsertNextId(id);
    resampledCenterlinesRadiusArray->InsertNextValue(centerlinesRadiusArray->GetValue(cell->GetPointId(cell->GetNumberOfPoints()-1)));

    resampledCenterlinesCellArray->InsertNextCell(resampledCell);
    }
  resampledCenterlines->SetPoints(resampledCenterlinesPoints);
  resampledCenterlines->SetLines(resampledCenterlinesCellArray);
  resampledCenterlines->GetPointData()->AddArray(resampledCenterlinesRadiusArray);
#if (VTK_MAJOR_VERSION <= 5)
  resampledCenterlines->Update();
#endif

  output->ShallowCopy(resampledCenterlines);

  resampledCenterlines->Delete();
  resampledCenterlinesPoints->Delete();
  resampledCenterlinesCellArray->Delete();
  resampledCenterlinesRadiusArray->Delete();
  resampledCell->Delete();
}

void vtkSVCenterlines::ReverseCenterlines()
{
  vtkPolyData* output = this->GetOutput();

  vtkCellArray* reversedCenterlinesCellArray = vtkCellArray::New();
  vtkIdList* reversedCell = vtkIdList::New();

  for (int k=0; k<output->GetNumberOfCells(); k++)
    {
    vtkCell* cell = output->GetCell(k);

    reversedCell->Initialize();

    vtkIdType numberOfCellPoints = cell->GetNumberOfPoints();

    for (int i=0; i<numberOfCellPoints; i++)
      {
      vtkIdType id = cell->GetPointId(numberOfCellPoints-1-i);
      reversedCell->InsertNextId(id);
      }
    reversedCenterlinesCellArray->InsertNextCell(reversedCell);
    }

  output->SetLines(reversedCenterlinesCellArray);

  reversedCell->Delete();
  reversedCenterlinesCellArray->Delete();
}

void vtkSVCenterlines::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

int vtkSVCenterlines::PruneVoronoiDiagram(vtkPolyData *inTriPd,
                                                    vtkPolyData *inEdgePd,
                                                    vtkPolyData *outTriPd,
                                                    vtkPolyData *outEdgePd,
                                                    std::string medialEdgeArrayName)
{
  int dontTouch = 0;
  if (medialEdgeArrayName != "")
    dontTouch = 1;

  vtkPolyData *tmpTriPd = vtkPolyData::New(); tmpTriPd->DeepCopy(inTriPd);
  vtkPolyData *tmpEdgePd = vtkPolyData::New(); tmpEdgePd->DeepCopy(inEdgePd);
  outTriPd->DeepCopy(inTriPd);
  outEdgePd->DeepCopy(inEdgePd);

  int numCells = tmpTriPd->GetNumberOfCells();
  int numPts = tmpTriPd->GetNumberOfPoints();
  int numEdgeCells = tmpEdgePd->GetNumberOfCells();
  int numEdgePts = tmpEdgePd->GetNumberOfPoints();

  std::vector<int> deletedCell(numCells, 0);
  std::vector<int> deletedEdge(numEdgeCells, 0);

  vtkIntArray *tmpEdgeArray = vtkIntArray::New();
  tmpEdgeArray->SetNumberOfTuples(numEdgeCells);
  tmpEdgeArray->SetName("TmpInternalIds");
  for (int i=0; i<numEdgeCells; i++)
    tmpEdgeArray->SetTuple1(i, i);
  tmpEdgePd->GetCellData()->AddArray(tmpEdgeArray);
  tmpEdgeArray->Delete();

  vtkIntArray *edgeRemoveIterArray = vtkIntArray::New();
  edgeRemoveIterArray->SetNumberOfTuples(numEdgeCells);
  edgeRemoveIterArray->FillComponent(0, -1);
  edgeRemoveIterArray->SetName("RemovalIteration");

  vtkIntArray *edgeIsolatedIterArray = vtkIntArray::New();
  edgeIsolatedIterArray->SetNumberOfTuples(numEdgeCells);
  edgeIsolatedIterArray->FillComponent(0, -1);
  edgeIsolatedIterArray->SetName("IsolatedIteration");
  tmpEdgePd->GetCellData()->AddArray(edgeIsolatedIterArray);
  edgeIsolatedIterArray->Delete();

  vtkIntArray *endIsolatedIterArray = vtkIntArray::New();
  endIsolatedIterArray->SetNumberOfTuples(numEdgeCells);
  endIsolatedIterArray->FillComponent(0, -1);
  endIsolatedIterArray->SetName("IsolatedIteration");

  vtkIntArray *removeIterArray = vtkIntArray::New();
  removeIterArray->SetNumberOfTuples(numCells);
  removeIterArray->FillComponent(0, -1);
  removeIterArray->SetName("RemovalIteration");

  int iter = 0;
  int nDelTris = 0;
  int nDelEdges = 0;
  int nIsolated = 0;
  while ( nDelTris > 0 || nDelEdges > 0 || iter == 0 )
  {
    std::vector<int> tmpDeletedCells;
    std::vector<int> tmpDeletedEdges;
    // --------------------------------------------------------------
    // Do edges before
    nDelEdges = 0;
    for (int i=0; i<tmpEdgePd->GetNumberOfCells(); i++)
    {
      if (!deletedEdge[i])
      {
        vtkIdType npts, *pts;
        tmpEdgePd->GetCellPoints(i, npts, pts);

        if (npts == 2)
        {
          vtkIdList *pt0CellIds = vtkIdList::New();
          vtkIdList *pt1CellIds = vtkIdList::New();

          tmpEdgePd->GetPointCells(pts[0], pt0CellIds);
          tmpEdgePd->GetPointCells(pts[1], pt1CellIds);

          int numNotDeletedNeighbors0 = 0;
          int numNotDeletedNeighbors1 = 0;
          for (int j=0; j<pt0CellIds->GetNumberOfIds(); j++)
          {
            if (!deletedEdge[pt0CellIds->GetId(j)])
              numNotDeletedNeighbors0++;
          }
          for (int j=0; j<pt1CellIds->GetNumberOfIds(); j++)
          {
            if (!deletedEdge[pt1CellIds->GetId(j)])
              numNotDeletedNeighbors1++;
          }

          //if (pt0CellIds->GetNumberOfIds() == 1 ||
          //    pt1CellIds->GetNumberOfIds() == 1)
          if (numNotDeletedNeighbors0 == 1 ||
              numNotDeletedNeighbors1 == 1)
          {
            int delEdge = 1;
            if (dontTouch)
            {
              int isMedEdge = tmpEdgePd->GetCellData()->GetArray("MedialEdges")->GetTuple1(i);
              if (isMedEdge == 1)
              {
                delEdge = 0;
              }
            }

            if (delEdge)
            {
              nDelEdges++;
              tmpDeletedEdges.push_back(i);
              //tmpEdgePd->DeleteCell(i);
              int origEdgeId = tmpEdgePd->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(i);
              //fprintf(stdout,"DELETING EDGE: %d\n", origEdgeId);
              edgeRemoveIterArray->SetTuple1(origEdgeId, iter);
            }
          }

          pt0CellIds->Delete();
          pt1CellIds->Delete();
        }
      }
    }
    if (iter == 0)
      fprintf(stdout,"NUM DEL SHOULD BE 0 TO START: %d\n", nDelEdges);
    // --------------------------------------------------------------
    nDelTris = 0;
    for (int i=0; i<tmpTriPd->GetNumberOfCells(); i++)
    {
      if (!deletedCell[i])
      {
        vtkIdType npts, *pts;
        tmpTriPd->GetCellPoints(i, npts, pts);

        if (npts == 3)
        {
          vtkIdList *openEdges = vtkIdList::New();
          for (int j=0; j<npts; j++)
          {
            int ptId0 = pts[j];
            int ptId1 = pts[(j+1)%npts];

            vtkIdList *cellNeighborIds = vtkIdList::New();
            tmpTriPd->GetCellEdgeNeighbors(i, ptId0, ptId1, cellNeighborIds);

            if (cellNeighborIds->GetNumberOfIds() == 0)
              openEdges->InsertNextId(j);
            else
            {
              int numDeletedNeighbors = 0;
              for (int k=0; k<cellNeighborIds->GetNumberOfIds(); k++)
              {
                if (deletedCell[cellNeighborIds->GetId(k)])
                  numDeletedNeighbors++;
              }
              if (numDeletedNeighbors == cellNeighborIds->GetNumberOfIds())
                openEdges->InsertNextId(j);
            }
            cellNeighborIds->Delete();
          }
          if (openEdges->GetNumberOfIds() == 3)
          {
            nDelTris++;
            tmpDeletedCells.push_back(i);
            //tmpTriPd->DeleteCell(i);
            int origCellId = tmpTriPd->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(i);
            //fprintf(stdout,"DELETING CELL: %d\n", origCellId);
            removeIterArray->SetTuple1(origCellId, iter);

            // --------------------------------------------------------------
            // Remove on edge pd
            int ptId0 = pts[0];
            int ptId1 = pts[1];
            vtkIdList *pointIds = vtkIdList::New();
            pointIds->SetNumberOfIds(2);
            pointIds->SetId(0, ptId0);
            pointIds->SetId(1, ptId1);

            vtkIdList *edgeCell = vtkIdList::New();
            tmpEdgePd->GetCellNeighbors(-1, pointIds, edgeCell);

            if (edgeCell->GetNumberOfIds() != 1)
              fprintf(stderr,"WEE HAVE PROBLEMMM 0: %d!\n", edgeCell->GetNumberOfIds());
            else
            {
              int delEdge = 1;
              if (dontTouch)
              {
                int isMedEdge = tmpEdgePd->GetCellData()->GetArray("MedialEdges")->GetTuple1(edgeCell->GetId(0));
                if (isMedEdge == 1)
                {
                  delEdge = 0;
                }
              }

              if (delEdge)
              {
                nDelEdges++;
                tmpDeletedEdges.push_back(edgeCell->GetId(0));
                //tmpEdgePd->DeleteCell(edgeCell->GetId(0));
                int origEdgeId = tmpEdgePd->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(edgeCell->GetId(0));
                //fprintf(stdout,"DELETING EDGE: %d\n", origEdgeId);
                edgeRemoveIterArray->SetTuple1(origEdgeId, iter);
              }
            }

            pointIds->Delete();
            edgeCell->Delete();
            // --------------------------------------------------------------
          }
          else if (openEdges->GetNumberOfIds() == 2)
          {
            nDelTris++;
            int loc;
            for (int j=0; j<npts; j++)
            {
              if (j != openEdges->GetId(0) && j != openEdges->GetId(1))
                loc = j;
            }

            int ptId0 = pts[loc];
            int ptId1 = pts[(loc+1)%npts];
            int ptId2 = pts[(loc+2)%npts];

            tmpDeletedCells.push_back(i);
            //tmpTriPd->DeleteCell(i);
            int origCellId = tmpTriPd->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(i);
            //fprintf(stdout,"DELETING CELL: %d\n", origCellId);
            removeIterArray->SetTuple1(origCellId, iter);

            // --------------------------------------------------------------
            // Remove on edge pd
            vtkIdList *pointIds = vtkIdList::New();
            pointIds->SetNumberOfIds(2);
            pointIds->SetId(0, ptId0);
            pointIds->SetId(1, ptId2);


            vtkIdList *edgeCell = vtkIdList::New();
            tmpEdgePd->GetCellNeighbors(-1, pointIds, edgeCell);

            if (edgeCell->GetNumberOfIds() != 1)
              fprintf(stderr,"WEE HAVE PROBLEMMM 1: %d!\n", edgeCell->GetNumberOfIds());
            else
            {
              int delEdge = 1;
              if (dontTouch)
              {
                int isMedEdge = tmpEdgePd->GetCellData()->GetArray("MedialEdges")->GetTuple1(edgeCell->GetId(0));
                if (isMedEdge == 1)
                {
                  delEdge = 0;
                }
              }

              if (delEdge)
              {
                nDelEdges++;
                tmpDeletedEdges.push_back(edgeCell->GetId(0));
                //tmpEdgePd->DeleteCell(edgeCell->GetId(0));
                int origEdgeId = tmpEdgePd->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(edgeCell->GetId(0));
                //fprintf(stdout,"DELETING EDGE: %d\n", origEdgeId);
                edgeRemoveIterArray->SetTuple1(origEdgeId, iter);
              }
            }

            pointIds->Delete();
            edgeCell->Delete();

            // --------------------------------------------------------------
          }
          else if (openEdges->GetNumberOfIds() == 1)
          {
            nDelTris++;
            int loc = openEdges->GetId(0);

            int ptId0 = pts[loc];
            int ptId1 = pts[(loc+1)%npts];
            int ptId2 = pts[(loc+2)%npts];

            tmpDeletedCells.push_back(i);
            //tmpTriPd->DeleteCell(i);
            int origCellId = tmpTriPd->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(i);
            //fprintf(stdout,"DELETING CELL: %d\n", origCellId);
            removeIterArray->SetTuple1(origCellId, iter);

            // --------------------------------------------------------------
            // Remove on edge pd
            vtkIdList *pointIds = vtkIdList::New();
            pointIds->SetNumberOfIds(2);
            pointIds->SetId(0, ptId0);
            pointIds->SetId(1, ptId1);

            vtkIdList *edgeCell = vtkIdList::New();
            tmpEdgePd->GetCellNeighbors(-1, pointIds, edgeCell);

            if (edgeCell->GetNumberOfIds() != 1)
              fprintf(stderr,"WEE HAVE PROBLEMMM 2: %d!\n", edgeCell->GetNumberOfIds());
            else
            {
              int delEdge = 1;
              if (dontTouch)
              {
                int isMedEdge = tmpEdgePd->GetCellData()->GetArray("MedialEdges")->GetTuple1(edgeCell->GetId(0));
                if (isMedEdge == 1)
                {
                  delEdge = 0;
                }
              }

              if (delEdge)
              {
                nDelEdges++;
                tmpDeletedEdges.push_back(edgeCell->GetId(0));
                //tmpEdgePd->DeleteCell(edgeCell->GetId(0));
                int origEdgeId = tmpEdgePd->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(edgeCell->GetId(0));
                //fprintf(stdout,"DELETING EDGE: %d\n", origEdgeId);
                edgeRemoveIterArray->SetTuple1(origEdgeId, iter);
              }
            }

            pointIds->Delete();
            edgeCell->Delete();
            // --------------------------------------------------------------
          }
          openEdges->Delete();
        }
      }
    }
    fprintf(stdout,"ITER %d, TRIS REMOVED: %d, EDGES REMOVED: %d\n", iter, nDelTris, nDelEdges);

    //tmpTriPd->RemoveDeletedCells();
    //tmpTriPd->BuildLinks();

    //tmpEdgePd->RemoveDeletedCells();
    //tmpEdgePd->BuildLinks();

    for (int i=0; i<tmpDeletedCells.size(); i++)
      deletedCell[tmpDeletedCells[i]] = 1;
    for (int i=0; i<tmpDeletedEdges.size(); i++)
      deletedEdge[tmpDeletedEdges[i]] = 1;

    //std::string fn = "/Users/adamupdegrove/Desktop/tmp/VORONOI_PRUNE" + std::to_string(iter) + ".vtp";
    //vtkXMLPolyDataWriter *inWriter = vtkXMLPolyDataWriter::New();
    //inWriter->SetInputData(tmpTriPd);
    //inWriter->SetFileName(fn.c_str());
    //inWriter->Write();
    //inWriter->Delete();

    //std::string inefn = "/Users/adamupdegrove/Desktop/tmp/VORONOI_EDGE_PRUNE" + std::to_string(iter) + ".vtp";
    //vtkXMLPolyDataWriter *inEdgeWriter = vtkXMLPolyDataWriter::New();
    //inEdgeWriter->SetInputData(tmpEdgePd);
    //inEdgeWriter->SetFileName(inefn.c_str());
    //inEdgeWriter->Write();
    //inEdgeWriter->Delete();

    // --------------------------------------------------------------
    // Now add to edge isolated list
    if (nIsolated != tmpEdgePd->GetNumberOfCells())
    {
      for (int i=0; i<tmpEdgePd->GetNumberOfCells(); i++)
      {
        int currVal = tmpEdgePd->GetCellData()->GetArray("IsolatedIteration")->GetTuple1(i);

        if (currVal == -1)
        {
          vtkIdType npts, *pts;
          tmpEdgePd->GetCellPoints(i, npts, pts);

          if (npts == 2)
          {
            vtkIdList *pointIds = vtkIdList::New();
            pointIds->SetNumberOfIds(2);
            pointIds->SetId(0, pts[0]);
            pointIds->SetId(1, pts[1]);

            vtkIdList *edgeCellIds = vtkIdList::New();

            tmpTriPd->GetCellNeighbors(-1, pointIds, edgeCellIds);

            int numDeletedNeighbors = 0;
            for (int j=0; j<edgeCellIds->GetNumberOfIds(); j++)
            {
              if (deletedCell[edgeCellIds->GetId(j)])
                numDeletedNeighbors++;
            }

            //if (edgeCellIds->GetNumberOfIds() == 0)
            if (numDeletedNeighbors == edgeCellIds->GetNumberOfIds())
            {
              int origEdgeId = tmpEdgePd->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(i);
              endIsolatedIterArray->SetTuple1(origEdgeId, iter);
              tmpEdgePd->GetCellData()->GetArray("IsolatedIteration")->SetTuple1(i, iter);
              nIsolated++;
            }
            pointIds->Delete();
            edgeCellIds->Delete();
          }
        }
      }
    }
    // --------------------------------------------------------------

    iter++;
  }

  tmpTriPd->Delete();

  outTriPd->GetCellData()->AddArray(removeIterArray);
  removeIterArray->Delete();

  outEdgePd->GetCellData()->AddArray(edgeRemoveIterArray);
  outEdgePd->GetCellData()->AddArray(endIsolatedIterArray);

  vtkIntArray *mAbsArray = vtkIntArray::New();
  mAbsArray->SetNumberOfTuples(numEdgeCells);
  mAbsArray->FillComponent(0, -1);
  mAbsArray->SetName("MAbs");

  vtkDoubleArray *mRelArray = vtkDoubleArray::New();
  mRelArray->SetNumberOfTuples(numEdgeCells);
  mRelArray->FillComponent(0, -1.0);
  mRelArray->SetName("MRel");

  for (int i=0; i<numEdgeCells; i++)
  {
    double currIVal = endIsolatedIterArray->GetTuple1(i);
    double edgeRVal = edgeRemoveIterArray->GetTuple1(i);

    if (currIVal == -1)
      endIsolatedIterArray->SetTuple1(i, 0);
    if (edgeRVal == -1)
      edgeRemoveIterArray->SetTuple1(i, iter);
  }

  for (int i=0; i<numEdgeCells; i++)
  {
    double iVal = endIsolatedIterArray->GetTuple1(i);
    double rVal = edgeRemoveIterArray->GetTuple1(i);

    int mAbsVal = rVal - iVal;
    double mRelVal = 1.0 - ((iVal+1)/(rVal+1));

    mAbsArray->SetTuple1(i, mAbsVal);
    mRelArray->SetTuple1(i, mRelVal);
  }

  outEdgePd->GetCellData()->AddArray(mAbsArray);
  outEdgePd->GetCellData()->AddArray(mRelArray);

  mAbsArray->Delete();
  mRelArray->Delete();
  edgeRemoveIterArray->Delete();
  endIsolatedIterArray->Delete();

  return 1;
}
