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
#include "vtkPolyLine.h"
#include "vtkvmtkInternalTetrahedraExtractor.h"
#include "vtkvmtkVoronoiDiagram3D.h"
#include "vtkvmtkSimplifyVoronoiDiagram.h"
#include "vtkArrayCalculator.h"
#include "vtkAppendPolyData.h"
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
#include "vtkSmartPointer.h"
#include "vtkIdList.h"
#include "vtkIdFilter.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkVersion.h"

#include "vtkSVIOUtils.h"
#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"

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

  if (!this->SourceSeedIds)
  {
    vtkDebugMacro(<< "No SourceSeedIds set.");
  }
  if (this->SourceSeedIds)
  {
    if (this->SourceSeedIds->GetNumberOfIds() != 1)
    {
      fprintf(stderr,"Only one source seed can be provided with this method\n");
      return SV_ERROR;
    }
  }

  if (!this->TargetSeedIds)
  {
    vtkDebugMacro(<< "No TargetSeedIds set.");
  }

  if (!this->RadiusArrayName)
  {
    vtkErrorMacro(<< "No RadiusArrayName set.");
    return SV_ERROR;
  }

  if (!this->GenerateDelaunayTessellation && !this->DelaunayTessellation)
  {
    vtkErrorMacro(<< "GenerateDelaunayTessellation is off but a DelaunayTessellation has not been set.");
    return SV_ERROR;
  }

  vtkNew(vtkPolyDataNormals, surfaceNormals);
  surfaceNormals->SetInputData(input);
  surfaceNormals->SplittingOff();
  surfaceNormals->AutoOrientNormalsOn();
  surfaceNormals->SetFlipNormals(this->FlipNormals);
  surfaceNormals->ComputePointNormalsOn();
  surfaceNormals->ConsistencyOn();
  surfaceNormals->Update();

  fprintf(stdout,"GENERATING DELAUNAY TESSELATION...\n");
  if (this->GenerateDelaunayTessellation)
  {
    vtkNew(vtkDelaunay3D, delaunayTessellator);
    delaunayTessellator->CreateDefaultLocator();
    delaunayTessellator->SetInputConnection(surfaceNormals->GetOutputPort());
    delaunayTessellator->SetTolerance(this->DelaunayTolerance);
    delaunayTessellator->Update();

    vtkUnstructuredGrid* delaunay = delaunayTessellator->GetOutput();
    delaunay->GetPointData()->AddArray(surfaceNormals->GetOutput()->GetPointData()->GetNormals());

    vtkNew(vtkvmtkInternalTetrahedraExtractor, internalTetrahedraExtractor);
    internalTetrahedraExtractor->SetInputConnection(delaunayTessellator->GetOutputPort());
    internalTetrahedraExtractor->SetOutwardNormalsArrayName(surfaceNormals->GetOutput()->GetPointData()->GetNormals()->GetName());
    if (this->CapCenterIds)
    {
      internalTetrahedraExtractor->UseCapsOn();
      internalTetrahedraExtractor->SetCapCenterIds(this->CapCenterIds);
    }
    internalTetrahedraExtractor->Update();

    this->DelaunayTessellation = internalTetrahedraExtractor->GetOutput();
    this->DelaunayTessellation->Register(this);
  }

  fprintf(stdout,"GENERATING VORONOI DIAGRAM...\n");
  vtkNew(vtkvmtkVoronoiDiagram3D, voronoiDiagramFilter);
  voronoiDiagramFilter->SetInputData(this->DelaunayTessellation);
  voronoiDiagramFilter->SetRadiusArrayName(this->RadiusArrayName);
  voronoiDiagramFilter->Update();

  this->PoleIds->DeepCopy(voronoiDiagramFilter->GetPoleIds());

  vtkPolyData* voronoiDiagram = voronoiDiagramFilter->GetOutput();
  if (this->SimplifyVoronoi)
  {
    vtkNew(vtkvmtkSimplifyVoronoiDiagram, voronoiDiagramSimplifier);
    voronoiDiagramSimplifier->SetInputConnection(voronoiDiagramFilter->GetOutputPort());
    voronoiDiagramSimplifier->SetUnremovablePointIds(voronoiDiagramFilter->GetPoleIds());
    voronoiDiagramSimplifier->Update();
    voronoiDiagram = voronoiDiagramSimplifier->GetOutput();
    voronoiDiagram->Register(this);
  }

  // ------------------------------------------------------------------------
  // Calculating genus
  // Start edge insertion for edge table
  fprintf(stdout,"CHECKING INPUT...\n");
  input->BuildLinks();
  int numNonTriangleCells, numNonManifoldEdges, numOpenEdges, surfaceGenus;
  vtkSVGeneralUtils::CheckSurface(input, numNonTriangleCells, numNonManifoldEdges, numOpenEdges, surfaceGenus);

  fprintf(stdout,"NUM NON TRIANGLE CELLS: %d\n", numNonTriangleCells);
  fprintf(stdout,"NUM NON MANIFOLD EDGES: %d\n", numNonManifoldEdges);
  fprintf(stdout,"NUM OPEN EDGES: %d\n", numOpenEdges);
  fprintf(stdout,"SURFACE GENUS: %d\n", surfaceGenus);

  if (numNonTriangleCells > 0)
  {
    vtkErrorMacro("Surface contains non-triangular cells. Number of non-triangular cells: " << numNonTriangleCells);
    return SV_ERROR;
  }
  if (numNonManifoldEdges > 0)
  {
    vtkErrorMacro("Surface contains non-manifold edges. Number of non-manifold edges: " << numNonManifoldEdges);
    return SV_ERROR;
  }
  if (numOpenEdges > 0)
  {
    vtkErrorMacro("Surface contains open edges. Number of open edges: " << numOpenEdges);
    return SV_ERROR;
  }
  if (surfaceGenus > 0)
  {
    vtkErrorMacro("Surface genus is greater than 0. Surface genus is: " << surfaceGenus);
    return SV_ERROR;
  }

  // ------------------------------------------------------------------------
  // Set up for pruning
  vtkNew(vtkTriangleFilter, triangulator);
  triangulator->SetInputData(voronoiDiagram);
  triangulator->Update();

  vtkNew(vtkPolyData, triPd);
  triPd->DeepCopy(triangulator->GetOutput());
  triPd->BuildLinks();

  int numCells = triPd->GetNumberOfCells();
  int numPts = triPd->GetNumberOfPoints();

  vtkNew(vtkPolyData, edgePd);
  vtkSVGeneralUtils::GetEdgePolyData(triPd, edgePd);
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Pruning voronoi diagram
  fprintf(stdout,"PRUNING VORONOI DIAGRAM...\n");
  vtkNew(vtkPolyData, newEdgePd);
  vtkNew(vtkPolyData, newTriPd);
  this->PruneVoronoiDiagram(triPd, edgePd, newTriPd, newEdgePd, "");
  fprintf(stdout,"DONE COMPUTING REMOVAL ITERATIONS\n");
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Give temporary ids for thresholding
  vtkNew(vtkIntArray, tmpEdgeArray);
  tmpEdgeArray->SetNumberOfTuples(newEdgePd->GetNumberOfCells());
  tmpEdgeArray->SetName("TmpInternalIds");
  for (int i=0; i<newEdgePd->GetNumberOfCells(); i++)
    tmpEdgeArray->SetTuple1(i, i);
  newEdgePd->GetCellData()->AddArray(tmpEdgeArray);
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Threshold based on absolute retention
  int mAbsThr = 3; // TODO: FIGURE OUT GOOD VALUE FOR THIS!!!
  double mAbsRange[2];
  newEdgePd->GetCellData()->GetArray("MAbs")->GetRange(mAbsRange);
  vtkNew(vtkThreshold, mAbsThresholder);
  mAbsThresholder->SetInputData(newEdgePd);
  mAbsThresholder->SetInputArrayToProcess(0, 0, 0, 1, "MAbs");
  mAbsThresholder->ThresholdBetween(mAbsThr, mAbsRange[1]);
  mAbsThresholder->Update();
  //fprintf(stdout,"Thresholded MAbs: %d\n", mAbsThresholder->GetOutput()->GetNumberOfCells());
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Threshold based on relative retention
  double mRelThr = 0.5;
  double mRelRange[2];
  newEdgePd->GetCellData()->GetArray("MRel")->GetRange(mRelRange);
  vtkNew(vtkThreshold, mRelThresholder);
  mRelThresholder->SetInputData(mAbsThresholder->GetOutput());
  mRelThresholder->SetInputArrayToProcess(0, 0, 0, 1, "MRel");
  mRelThresholder->ThresholdBetween(mRelThr, mRelRange[1]);
  mRelThresholder->Update();
  //fprintf(stdout,"Thresholded MRel: %d\n", mRelThresholder->GetOutput()->GetNumberOfCells());
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Get the medial edges
  vtkNew(vtkConnectivityFilter, connector);
  connector->SetInputData(mRelThresholder->GetOutput());
  connector->SetExtractionModeToAllRegions();
  connector->ColorRegionsOn();
  connector->Update();

  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(connector->GetOutput());
  surfacer->Update();

  vtkNew(vtkPolyData, leftOver);
  leftOver->DeepCopy(surfacer->GetOutput());

  vtkNew(vtkIntArray, leaveArray);
  leaveArray->SetNumberOfTuples(edgePd->GetNumberOfCells());
  leaveArray->SetName("MedialEdges");
  for (int i=0; i<edgePd->GetNumberOfCells(); i++)
    leaveArray->SetTuple1(i, 0);

  int connectThr = 5;
  vtkNew(vtkThreshold, regionThresholder);
  regionThresholder->SetInputData(leftOver);
  regionThresholder->SetInputArrayToProcess(0, 0, 0, 1, "RegionId");

  for (int i=0; i<connector->GetNumberOfExtractedRegions(); i++)
  {
    regionThresholder->ThresholdBetween(i, i);
    regionThresholder->Update();

    fprintf(stdout,"THRESHOLDED REGION %d: %d\n", i, regionThresholder->GetOutput()->GetNumberOfCells());
    if (regionThresholder->GetOutput()->GetNumberOfCells() > connectThr)
    {
      for (int j=0; j<regionThresholder->GetOutput()->GetNumberOfCells(); j++)
      {
        int origCellId = regionThresholder->GetOutput()->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(j);
        leaveArray->SetTuple1(origCellId, 1);
      }
    }
  }

  edgePd->GetCellData()->AddArray(leaveArray);
  // ------------------------------------------------------------------------


  // ------------------------------------------------------------------------
  // Now prune again
  fprintf(stdout,"PRUNE VORONOI DIAGRAM WHILE RETAINING MEDIAL EDGES\n");
  vtkNew(vtkPolyData, nextTriPd);
  vtkNew(vtkPolyData, nextEdgePd);
  this->PruneVoronoiDiagram(triPd, edgePd, nextTriPd, nextEdgePd, "MedialEdges");
  fprintf(stdout,"DONE COMPUTING REMOVAL ITERATIONS\n");
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Threshold last removal iteration cells to get centerlines
  double finalRange[2];
  nextEdgePd->GetCellData()->GetArray("RemovalIteration")->GetRange(finalRange);

  vtkNew(vtkIntArray, keepCellArray);
  keepCellArray->SetNumberOfTuples(triPd->GetNumberOfCells());
  keepCellArray->SetName("KeepCellArray");
  for (int i=0; i<triPd->GetNumberOfCells(); i++)
    keepCellArray->SetTuple1(i, 0);

  vtkNew(vtkIdList, pointCellIds);
  for (int i=0; i<nextEdgePd->GetNumberOfCells(); i++)
  {
    int rVal = nextEdgePd->GetCellData()->GetArray("RemovalIteration")->GetTuple1(i);

    if (rVal == finalRange[1])
    {
      vtkIdType npts, *pts;
      nextEdgePd->GetCellPoints(i, npts, pts);

      for (int j=0; j<npts; j++)
      {
        triPd->GetPointCells(pts[j], pointCellIds);

        for (int k=0; k<pointCellIds->GetNumberOfIds(); k++)
        {
          int cellId = pointCellIds->GetId(k);
          keepCellArray->SetTuple1(pointCellIds->GetId(k), 1);
        }

      }
    }
  }

  triPd->GetCellData()->AddArray(keepCellArray);
  nextEdgePd->GetPointData()->AddArray(triPd->GetPointData()->GetArray(this->RadiusArrayName));
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Threshold last removal iteration cells to get centerlines
  double radRange[2];
  voronoiDiagram->GetPointData()->GetArray(this->RadiusArrayName)->GetRange(radRange);
  //fprintf(stdout,"RADIUS ARRAY RANGE: %.6f %.6f\n", radRange[0], radRange[1]);
  //fprintf(stdout,"1/RADIUS ARRAY RANGE: %.6f %.6f\n", 1./radRange[0], 1./radRange[1]);

  vtkNew(vtkIntArray, tmpPtArray);
  tmpPtArray->Reset();
  tmpPtArray->SetNumberOfTuples(nextEdgePd->GetNumberOfPoints());
  tmpPtArray->SetName("TmpInternalIds");
  for (int i=0; i<nextEdgePd->GetNumberOfPoints(); i++)
    tmpPtArray->SetTuple1(i, i);

  nextEdgePd->GetPointData()->AddArray(tmpPtArray);

  vtkNew(vtkThreshold, finalThreshold);
  finalThreshold->SetInputData(nextEdgePd);
  finalThreshold->SetInputArrayToProcess(0, 0, 0, 1, "RemovalIteration");
  finalThreshold->ThresholdBetween(finalRange[1], finalRange[1]);
  finalThreshold->Update();

  surfacer->SetInputData(finalThreshold->GetOutput());
  surfacer->Update();

  vtkNew(vtkCleanPolyData, cleaner);
  cleaner->SetInputData(surfacer->GetOutput());
  cleaner->Update();

  vtkNew(vtkPolyData, linesPd);
  linesPd->DeepCopy(cleaner->GetOutput());
  linesPd->BuildLinks();
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Remove cells that may have potentially been reduced to vertices
  for (int i=0; i<linesPd->GetNumberOfCells(); i++)
  {
    if (linesPd->GetCellType(i) != VTK_LINE)
    {
      linesPd->DeleteCell(i);
    }
  }
  linesPd->RemoveDeletedCells();
  linesPd->BuildLinks();
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Find the end points of the lines
	int firstVertex = -1;
  std::vector<std::vector<int> > connectedEdgePts(linesPd->GetNumberOfPoints());
  vtkNew(vtkIdList, linesEndPointIds);
  vtkNew(vtkPoints, linesEndPoints);

  this->GetLinesEndPoints(linesPd, linesEndPointIds, linesEndPoints, connectedEdgePts, firstVertex);
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Set the end point locator
  vtkNew(vtkPointLocator, linesEndPointLocator);
  vtkNew(vtkPolyData, linesEndPointsPd);  linesEndPointsPd->SetPoints(linesEndPoints);
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Check if an end point found
  if (firstVertex == -1)
  {
    vtkErrorMacro("No first vertex found, lines must form loop");
    return SV_ERROR;
  }
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Set up to get edge data structures from leftover polydata
  std::vector<int> pointUsed(linesPd->GetNumberOfPoints(), 0);

  pointUsed[firstVertex] = 1;
  int startVertex = connectedEdgePts[firstVertex][0];

  std::vector<std::vector<int> > allEdges;
  std::vector<int> thisEdge;
  thisEdge.push_back(firstVertex);

  this->RecursiveGetPolylines(linesPd, connectedEdgePts, startVertex, pointUsed, allEdges, thisEdge);
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Check the output lines and use graph cycles to remove duplicate lines
  vtkNew(vtkIdList, allEndIds);
  std::vector<int> nodeCount;
  std::vector<int> needToDelete(allEdges.size(), 0);
  for (int i=0; i<allEdges.size(); i++)
  {
    int edgeSize = allEdges[i].size();
    int edgeId0 = allEdges[i][0];
    int edgeIdN = allEdges[i][edgeSize-1];

    int edge0IsId = allEndIds->IsId(edgeId0);
    int edgeNIsId = allEndIds->IsId(edgeIdN);

    if (edge0IsId != -1 && edgeNIsId != -1)
    {
      // Both edges of this node already in list, delete it!
      //fprintf(stdout,"BOTH EXIST, NEED TO DELETE %d %d\n", edgeId0, edgeIdN);
      needToDelete[i] = 1;
      nodeCount[edge0IsId]++;
      nodeCount[edgeNIsId]++;
    }
    else
    {
      if (edge0IsId == -1)
      {
        allEndIds->InsertNextId(edgeId0);
        nodeCount.push_back(1);
      }
      else
        nodeCount[edge0IsId]++;
      if (edgeNIsId == -1)
      {
        allEndIds->InsertNextId(edgeIdN);
        nodeCount.push_back(1);
      }
      else
        nodeCount[edgeNIsId]++;
    }
  }
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Delete cells based on what still needs to be removed
  vtkNew(vtkIdList, pointsEdgeId);
  vtkNew(vtkIdList, edgePointIds);
  std::vector<int> isDeleted(allEdges.size(), 0);

  int done = 0;
  while (!done)
  {
    this->RemoveMarkedCells(linesPd, allEdges, needToDelete, isDeleted, allEndIds, nodeCount);

    std::vector<int> newNeedToDelete(allEdges.size(), 0);
    for (int i=0; i<needToDelete.size(); i++)
    {
      if (needToDelete[i] == 1)
      {
        int edgeSize = allEdges[i].size();

        int nodeId0 = allEndIds->IsId(allEdges[i][0]);
        int nodeIdN = allEndIds->IsId(allEdges[i][edgeSize-1]);

        if (nodeCount[nodeId0] == 1 || nodeCount[nodeIdN] == 1)
        {
          for (int j=0; j<allEdges.size(); j++)
          {
            if (isDeleted[j])
            {
              continue;
            }

            int delEdgeSize = allEdges[j].size();

            if (nodeCount[nodeId0] == 1)
            {
              if (allEdges[j][0] == allEdges[i][0] ||
                  allEdges[j][delEdgeSize-1] == allEdges[i][0])
              {
                newNeedToDelete[j] = 1;
              }
            }
            if (nodeCount[nodeIdN] == 1)
            {
              if (allEdges[j][0] == allEdges[i][edgeSize-1] ||
                  allEdges[j][delEdgeSize-1] == allEdges[i][edgeSize-1])
              {
                newNeedToDelete[j] = 1;
              }
            }
          }
        }
      }
    }

    done = 1;
    for (int i=0; i<needToDelete.size(); i++)
    {
      if (newNeedToDelete[i])
        done = 0;
      needToDelete[i] = newNeedToDelete[i];
    }
  }
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Delete very small entrance lines
  for (int i=0; i<allEdges.size(); i++)
  {
    int edgeSize = allEdges[i].size();
    int nodeId0 = allEndIds->IsId(allEdges[i][0]);
    int nodeIdN = allEndIds->IsId(allEdges[i][edgeSize-1]);

    if (edgeSize < connectThr)
    {
      if (nodeCount[nodeId0] <= 1 || nodeCount[nodeIdN] <= 1)
      {
        needToDelete[i] = 1;
      }
    }
  }

  this->RemoveMarkedCells(linesPd, allEdges, needToDelete, isDeleted, allEndIds, nodeCount);
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Delete cells based on given seeds
  if (this->SourceSeedIds)
  {
    if (this->SourceSeedIds->GetNumberOfIds() != 1)
    {
      fprintf(stdout,"Only one source seed can be provided with this method\n");
      return SV_ERROR;
    }

    linesEndPointLocator->SetDataSet(linesEndPointsPd);
    linesEndPointLocator->BuildLocator();

    std::vector<int> endPointUsed(linesEndPointIds->GetNumberOfIds(), 0);
    std::vector<int> deleteSeeds;
    for (int j=0; j<this->SourceSeedIds->GetNumberOfIds(); j++)
    {
      double sourcePt[3];
      if (this->CapCenterIds)
        input->GetPoint(this->CapCenterIds->GetId(this->SourceSeedIds->GetId(j)), sourcePt);
      else
        input->GetPoint(this->SourceSeedIds->GetId(j), sourcePt);

      int endPointId = linesEndPointLocator->FindClosestPoint(sourcePt);
      int linesPtId = linesEndPointIds->GetId(endPointId);

      if (endPointUsed[endPointId] == 1)
      {
        fprintf(stderr,"Two end lines found for different target seeds, target seeds too close\n");
        //return SV_ERROR;
      }
      else
      {
        endPointUsed[endPointId] = 1;
      }
    }

    if (this->TargetSeedIds)
    {
      int numSeeds = this->SourceSeedIds->GetNumberOfIds() + this->TargetSeedIds->GetNumberOfIds();

      if (numSeeds > linesEndPointIds->GetNumberOfIds() && this->TargetSeedIds)
      {
        fprintf(stdout,"More seeds given than found ends\n");
        vtkNew(vtkPointLocator, seedPointLocator);
        vtkNew(vtkPoints, seedPoints);
        vtkNew(vtkPolyData, seedPointsPd);  seedPointsPd->SetPoints(seedPoints);
        vtkNew(vtkIdList, seedPointIds);
        for (int j=0; j<this->TargetSeedIds->GetNumberOfIds(); j++)
        {
          double targetPt[3];
          if (this->CapCenterIds)
            input->GetPoint(this->CapCenterIds->GetId(this->TargetSeedIds->GetId(j)), targetPt);
          else
            input->GetPoint(this->TargetSeedIds->GetId(j), targetPt);

          seedPointsPd->GetPoints()->InsertNextPoint(targetPt);
          seedPointIds->InsertNextId(j);
        }
        seedPointLocator->SetDataSet(seedPointsPd);
        seedPointLocator->BuildLocator();

        for (int j=0; j<linesEndPointIds->GetNumberOfIds(); j++)
        {
          if (endPointUsed[j] == 1)
          {
            continue;
          }
          endPointUsed[j] = 1;

          int linesPtId = linesEndPointIds->GetId(j);

          double endPt[3];
          linesPd->GetPoint(linesPtId, endPt);

          int closestSeed = seedPointLocator->FindClosestPoint(endPt);
          int targetSeedId = seedPointIds->GetId(closestSeed);
        }
      }
      else
      {
        if (numSeeds == linesEndPointIds->GetNumberOfIds() && this->TargetSeedIds)
          fprintf(stdout,"Equal number of seeds and found ends\n");
        else if (numSeeds < linesEndPointIds->GetNumberOfIds() && this->TargetSeedIds)
          fprintf(stdout,"Less seeds given than found ends\n");

        for (int j=0; j<this->TargetSeedIds->GetNumberOfIds(); j++)
        {
          double targetPt[3];
          if (this->CapCenterIds)
            input->GetPoint(this->CapCenterIds->GetId(this->TargetSeedIds->GetId(j)), targetPt);
          else
            input->GetPoint(this->TargetSeedIds->GetId(j), targetPt);

          int endPointId = linesEndPointLocator->FindClosestPoint(targetPt);
          int linesPtId = linesEndPointIds->GetId(endPointId);

          if (endPointUsed[endPointId] == 1)
          {
            fprintf(stderr,"Two end lines found for different target seeds, target seeds too close\n");
            //return SV_ERROR;
          }
          else
          {
            endPointUsed[endPointId] = 1;
          }
        }
      }

      for (int i=0; i<endPointUsed.size(); i++)
      {
        if (endPointUsed[i] == 0)
        {
          fprintf(stdout,"End point %d was not used, going to remove from search list\n", i);
          int linesPtId = linesEndPointIds->GetId(i);
          deleteSeeds.push_back(linesPtId);
        }
      }
    }

    for (int i=0; i<deleteSeeds.size(); i++)
    {
      for (int j=0; j<allEdges.size(); j++)
      {
        if (isDeleted[j])
        {
          continue;
        }

        int edgeSize = allEdges[j].size();
        int nodeId0 = allEndIds->IsId(allEdges[j][0]);
        int nodeIdN = allEndIds->IsId(allEdges[j][edgeSize-1]);

        if (nodeId0 == deleteSeeds[i] || nodeIdN == deleteSeeds[i])
        {
          needToDelete[j] = 1;
        }
      }
    }

    this->RemoveMarkedCells(linesPd, allEdges, needToDelete, isDeleted, allEndIds, nodeCount);
  }

  // ------------------------------------------------------------------------
  // Remove deleted cells
  linesPd->RemoveDeletedCells();
  cleaner->SetInputData(linesPd);
  cleaner->Update();

  linesPd->DeepCopy(cleaner->GetOutput());
  linesPd->BuildLinks();
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Remove cells that may have potentially been reduced to vertices
  for (int i=0; i<linesPd->GetNumberOfCells(); i++)
  {
    if (linesPd->GetCellType(i) != VTK_LINE)
    {
      linesPd->DeleteCell(i);
    }
  }
  linesPd->RemoveDeletedCells();
  linesPd->BuildLinks();
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Check the genus of our final edges
  int nEdgeE = linesPd->GetNumberOfCells();
  int nEdgeV = linesPd->GetNumberOfPoints();
  int nFaces = nEdgeE - nEdgeV + 2;
  if (nFaces != 1)
  {
    vtkErrorMacro("After processing, centerline contains faces, or contains loops or cycles");
    return SV_ERROR;
  }
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Get end points
  this->GetLinesEndPoints(linesPd, linesEndPointIds, linesEndPoints, connectedEdgePts, firstVertex);
  linesEndPointsPd->SetPoints(linesEndPoints);
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Check if an end point found
  if (firstVertex == -1)
  {
    vtkErrorMacro("No first vertex found, lines must form loop");
    return SV_ERROR;
  }
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Get starting seed point
  if (this->SourceSeedIds)
  {
    if (this->SourceSeedIds->GetNumberOfIds() != 1)
    {
      fprintf(stdout,"Only one source seed can be provided with this method\n");
      return SV_ERROR;
    }

    linesEndPointLocator->SetDataSet(linesEndPointsPd);
    linesEndPointLocator->BuildLocator();

    double firstPt[3];
    if (this->CapCenterIds)
      input->GetPoint(this->CapCenterIds->GetId(this->SourceSeedIds->GetId(0)), firstPt);
    else
      input->GetPoint(this->SourceSeedIds->GetId(0), firstPt);

    int endPointId = linesEndPointLocator->FindClosestPoint(firstPt);
    int linesPtId = linesEndPointIds->GetId(endPointId);

    // If source seed given, use this as starting seed!
    firstVertex = linesPtId;
  }

  if (firstVertex == -1)
  {
    fprintf(stderr,"No first vertex found, lines cannot form loop\n");
    return SV_ERROR;
  }
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Re-get all the edges now that weve removed a bunch
  pointUsed.clear();
  pointUsed.resize(linesPd->GetNumberOfPoints(), 0);

  pointUsed[firstVertex] = 1;
  startVertex = connectedEdgePts[firstVertex][0];

  allEdges.clear();
  thisEdge.clear();
  thisEdge.push_back(firstVertex);

  this->RecursiveGetPolylines(linesPd, connectedEdgePts, startVertex, pointUsed, allEdges, thisEdge);
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Turn these edges into constructed lines from source end point to all other
  // end points
  std::vector<std::vector<int> > fullCenterlineEdges;

  int startEdge = 0;
  int front = allEdges[startEdge][0];
  int back  = allEdges[startEdge][allEdges[0].size()-1];

  this->RecursiveGetFullCenterlines(allEdges, fullCenterlineEdges, startEdge, front, back);
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Get the cost function
  vtkNew(vtkArrayCalculator, voronoiCostFunctionCalculator);
#if (VTK_MAJOR_VERSION <= 5)
  voronoiCostFunctionCalculator->SetInput(voronoiDiagram);
#else
  voronoiCostFunctionCalculator->SetInputData(voronoiDiagram);
#endif
  voronoiCostFunctionCalculator->SetAttributeModeToUsePointData();
  voronoiCostFunctionCalculator->AddScalarVariable("R",this->RadiusArrayName,0);
  voronoiCostFunctionCalculator->SetFunction(this->CostFunction);
  voronoiCostFunctionCalculator->SetResultArrayName(this->CostFunctionArrayName);
  voronoiCostFunctionCalculator->Update();

  surfacer->SetInputData(voronoiCostFunctionCalculator->GetOutput());
  surfacer->Update();
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // Get the actual seed points on the voronoi diagram
  std::vector<std::vector<int> > voronoiSeeds(allEdges.size());
  for (int i=0; i<allEdges.size(); i++)
  {
    int edgeSize = allEdges[i].size();
    int voronoiId0 = linesPd->GetPointData()->GetArray("TmpInternalIds")->GetTuple1(allEdges[i][0]);
    int voronoiId1 = linesPd->GetPointData()->GetArray("TmpInternalIds")->GetTuple1(allEdges[i][edgeSize-1]);
    voronoiSeeds[i].push_back(voronoiId0);
    voronoiSeeds[i].push_back(voronoiId1);
  }
  // ------------------------------------------------------------------------

  // ------------------------------------------------------------------------
  // If CapCenterIds, SourceSeedIds, or TargetSeedIds are given, we need to
  // process to get actual start points
  vtkNew(vtkIdList, voronoiCapIds);
  if (this->CapCenterIds)
  {
    this->FindVoronoiSeeds(this->DelaunayTessellation,this->CapCenterIds,surfaceNormals->GetOutput()->GetPointData()->GetNormals(),voronoiCapIds);
  }

  if (this->SourceSeedIds)
  {
    if (this->SourceSeedIds->GetNumberOfIds() != 1)
    {
      fprintf(stdout,"Only one source seed can be provided with this method\n");
      return SV_ERROR;
    }

    std::vector<int> endPointUsed(linesEndPointIds->GetNumberOfIds(), 0);
    for (int j=0; j<this->SourceSeedIds->GetNumberOfIds(); j++)
    {
      double sourcePt[3];
      if (this->CapCenterIds)
        input->GetPoint(this->CapCenterIds->GetId(this->SourceSeedIds->GetId(j)), sourcePt);
      else
        input->GetPoint(this->SourceSeedIds->GetId(j), sourcePt);

      int endPointId = linesEndPointLocator->FindClosestPoint(sourcePt);
      int linesPtId = linesEndPointIds->GetId(endPointId);
      int voronoiId = linesPd->GetPointData()->GetArray("TmpInternalIds")->GetTuple1(linesPtId);

      if (endPointUsed[endPointId] == 1)
      {
        fprintf(stderr,"Two end lines found for different target seeds, target seeds too close\n");
        //return SV_ERROR;
      }
      else
      {
        endPointUsed[endPointId] = 1;
      }

      for (int k=0; k<voronoiSeeds.size(); k++)
      {
        for (int l=0; l<voronoiSeeds[k].size(); l++)
        {
          if (voronoiSeeds[k][l] == voronoiId)
          {
            if (this->CapCenterIds)
            {
              voronoiSeeds[k][l] = voronoiCapIds->GetId(this->SourceSeedIds->GetId(j));
            }
            else
            {
              voronoiSeeds[k][l] = this->PoleIds->GetId(this->SourceSeedIds->GetId(j));
            }
          }
        }
      }
    }

    if (this->TargetSeedIds)
    {
      int numSeeds = this->SourceSeedIds->GetNumberOfIds() + this->TargetSeedIds->GetNumberOfIds();

      if (numSeeds > linesEndPointIds->GetNumberOfIds() && this->TargetSeedIds)
      {
        fprintf(stdout,"More seeds given than found ends\n");
        vtkNew(vtkPointLocator, seedPointLocator);
        vtkNew(vtkPoints, seedPoints);
        vtkNew(vtkPolyData, seedPointsPd);  seedPointsPd->SetPoints(seedPoints);
        vtkNew(vtkIdList, seedPointIds);
        for (int j=0; j<this->TargetSeedIds->GetNumberOfIds(); j++)
        {
          double targetPt[3];
          if (this->CapCenterIds)
          {
            input->GetPoint(this->CapCenterIds->GetId(this->TargetSeedIds->GetId(j)), targetPt);
          }
          else
          {
            input->GetPoint(this->TargetSeedIds->GetId(j), targetPt);
          }

          seedPointsPd->GetPoints()->InsertNextPoint(targetPt);
          seedPointIds->InsertNextId(j);
        }
        seedPointLocator->SetDataSet(seedPointsPd);
        seedPointLocator->BuildLocator();

        for (int j=0; j<linesEndPointIds->GetNumberOfIds(); j++)
        {
          if (endPointUsed[j] == 1)
          {
            continue;
          }
          endPointUsed[j] = 1;

          int linesPtId = linesEndPointIds->GetId(j);

          double endPt[3];
          linesPd->GetPoint(linesPtId, endPt);
          int voronoiId = linesPd->GetPointData()->GetArray("TmpInternalIds")->GetTuple1(linesPtId);

          int closestSeed = seedPointLocator->FindClosestPoint(endPt);
          int targetSeedId = seedPointIds->GetId(closestSeed);

          for (int k=0; k<voronoiSeeds.size(); k++)
          {
            for (int l=0; l<voronoiSeeds[k].size(); l++)
            {
              if (voronoiSeeds[k][l] == voronoiId)
              {
                if (this->CapCenterIds)
                {
                  voronoiSeeds[k][l] = voronoiCapIds->GetId(this->TargetSeedIds->GetId(targetSeedId));
                }
                else
                {
                  voronoiSeeds[k][l] = this->PoleIds->GetId(this->TargetSeedIds->GetId(targetSeedId));
                }
              }
            }
          }
        }
      }
      else
      {
        if (numSeeds == linesEndPointIds->GetNumberOfIds() && this->TargetSeedIds)
          fprintf(stdout,"Equal number of seeds and found ends\n");
        else if (numSeeds < linesEndPointIds->GetNumberOfIds() && this->TargetSeedIds)
          fprintf(stdout,"Less seeds given than found ends\n");

        for (int j=0; j<this->TargetSeedIds->GetNumberOfIds(); j++)
        {
          double targetPt[3];
          if (this->CapCenterIds)
          {
            input->GetPoint(this->CapCenterIds->GetId(this->TargetSeedIds->GetId(j)), targetPt);
          }
          else
          {
            input->GetPoint(this->TargetSeedIds->GetId(j), targetPt);
          }

          int endPointId = linesEndPointLocator->FindClosestPoint(targetPt);
          int linesPtId = linesEndPointIds->GetId(endPointId);
          int voronoiId = linesPd->GetPointData()->GetArray("TmpInternalIds")->GetTuple1(linesPtId);

          if (endPointUsed[endPointId] == 1)
          {
            fprintf(stderr,"Two end lines found for different target seeds, target seeds too close\n");
            //return SV_ERROR;
          }
          else
          {
            endPointUsed[endPointId] = 1;
          }

          for (int k=0; k<voronoiSeeds.size(); k++)
          {
            for (int l=0; l<voronoiSeeds[k].size(); l++)
            {
              if (voronoiSeeds[k][l] == voronoiId)
              {
                if (this->CapCenterIds)
                {
                  voronoiSeeds[k][l] = voronoiCapIds->GetId(this->TargetSeedIds->GetId(j));
                }
                else
                {
                  voronoiSeeds[k][l] = this->PoleIds->GetId(this->TargetSeedIds->GetId(j));
                }
              }
            }
          }
        }
      }
    }
  }
  // ------------------------------------------------------------------------

  vtkNew(vtkAppendPolyData, appender);

  vtkNew(vtkvmtkNonManifoldFastMarching, voronoiFastMarching);
  voronoiFastMarching->SetInputData(voronoiCostFunctionCalculator->GetOutput());
  voronoiFastMarching->SetCostFunctionArrayName(this->CostFunctionArrayName);
  voronoiFastMarching->SetSolutionArrayName(this->EikonalSolutionArrayName);
  voronoiFastMarching->SeedsBoundaryConditionsOn();
  for (int i=0; i<voronoiSeeds.size(); i++)
  {
    int voronoiId0 = voronoiSeeds[i][0];
    int voronoiIdN = voronoiSeeds[i][1];

    vtkNew(vtkIdList, voronoiSourceSeedIds);
    voronoiSourceSeedIds->SetNumberOfIds(1);
    voronoiSourceSeedIds->SetId(0, voronoiIdN);

    vtkNew(vtkIdList, voronoiTargetSeedIds);
    voronoiTargetSeedIds->SetNumberOfIds(1);
    voronoiTargetSeedIds->SetId(0, voronoiId0);

    fprintf(stdout,"DOING EDGE: %d %d\n", voronoiId0, voronoiIdN);
    fprintf(stdout,"MARCHING...\n");
    voronoiFastMarching->SetSeeds(voronoiSourceSeedIds);
    voronoiFastMarching->Update();
    fprintf(stdout,"Done\n");

    this->VoronoiDiagram->ShallowCopy(voronoiFastMarching->GetOutput());

    vtkNew(vtkvmtkSteepestDescentLineTracer, centerlineBacktracing);
    centerlineBacktracing->SetInputConnection(voronoiFastMarching->GetOutputPort());
    centerlineBacktracing->SetDataArrayName(this->RadiusArrayName);
    centerlineBacktracing->SetDescentArrayName(this->EikonalSolutionArrayName);
    centerlineBacktracing->SetEdgeArrayName(this->EdgeArrayName);
    centerlineBacktracing->SetEdgePCoordArrayName(this->EdgePCoordArrayName);
    centerlineBacktracing->SetSeeds(voronoiTargetSeedIds);
    centerlineBacktracing->MergePathsOff();
    centerlineBacktracing->StopOnTargetsOn();
    centerlineBacktracing->SetTargets(voronoiSourceSeedIds);
    fprintf(stdout,"BACKTRACKING...\n");
    centerlineBacktracing->Update();
    fprintf(stdout,"Done\n");

    appender->AddInputData(centerlineBacktracing->GetOutput());
  }
  appender->Update();

  vtkNew(vtkPolyData, currentLine);
  vtkNew(vtkCellArray, newCells);
  vtkNew(vtkPoints, newPoints);
  vtkNew(vtkPointData, newPointData);
  newPointData->CopyAllocate(appender->GetOutput()->GetPointData(),
                             appender->GetOutput()->GetNumberOfPoints());
  for (int i=0; i<fullCenterlineEdges.size(); i++)
  {
    vtkNew(vtkPolyLine, newLine);

    // First point if wanted
    if (this->AppendEndPointsToCenterlines)
    {
      if (this->SourceSeedIds)
      {
        // Enter source seed here
        double startPt[3];
        if (this->CapCenterIds)
          input->GetPoint(this->CapCenterIds->GetId(this->SourceSeedIds->GetId(0)), startPt);
        else
          input->GetPoint(this->SourceSeedIds->GetId(0), startPt);

        int newPointId = newPoints->InsertNextPoint(startPt);
        newLine->GetPointIds()->InsertNextId(newPointId);
        newPointData->CopyData(appender->GetInput(fullCenterlineEdges[i][0])->GetPointData(), 0, newPointId);
      }
    }

    for (int j=0; j<fullCenterlineEdges[i].size(); j++)
    {
      currentLine = appender->GetInput(fullCenterlineEdges[i][j]);

      int kStart = 1;
      if (j == 0)
        kStart = 0;
      for (int k=kStart; k<currentLine->GetNumberOfPoints(); k++)
      {
        double pt[3];
        currentLine->GetPoint(k, pt);

        int newPointId = newPoints->InsertNextPoint(pt);
        newLine->GetPointIds()->InsertNextId(newPointId);
        newPointData->CopyData(currentLine->GetPointData(), k, newPointId);
      }
    }

    // End point if wanted
    if (this->AppendEndPointsToCenterlines)
    {
      if (this->TargetSeedIds)
      {
        // Enter target seed here
        int edgeSize = fullCenterlineEdges[i].size();
        int polePointId = voronoiSeeds[fullCenterlineEdges[i][edgeSize-1]][1];
        int targetPointId;
        if (this->CapCenterIds)
          targetPointId = voronoiCapIds->IsId(polePointId);
        else
          targetPointId = this->PoleIds->IsId(polePointId);
        if (targetPointId != -1)
        {
          double endPt[3];
          if (this->CapCenterIds)
            input->GetPoint(this->CapCenterIds->GetId(targetPointId), endPt);
          else
            input->GetPoint(targetPointId, endPt);

          int numPtsInLine = appender->GetInput(fullCenterlineEdges[i][edgeSize-1])->GetNumberOfPoints();

          int newPointId = newPoints->InsertNextPoint(endPt);
          newLine->GetPointIds()->InsertNextId(newPointId);
          newPointData->CopyData(appender->GetInput(fullCenterlineEdges[i][edgeSize-1])->GetPointData(), numPtsInLine-1, newPointId);
        }
      }
    }

    newCells->InsertNextCell(newLine);
  }
  newPointData->Squeeze();

  vtkNew(vtkPolyData, finalLinesPd);
  finalLinesPd->SetPoints(newPoints);
  finalLinesPd->SetLines(newCells);
  finalLinesPd->GetPointData()->PassData(newPointData);

  output->ShallowCopy(finalLinesPd);

  if (this->CenterlineResampling)
    {
    this->ResampleCenterlines();
    }
  //this->ReverseCenterlines();
//
  return SV_OK;
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

  vtkNew(vtkPolyData, tmpTriPd); tmpTriPd->DeepCopy(inTriPd);
  vtkNew(vtkPolyData, tmpEdgePd); tmpEdgePd->DeepCopy(inEdgePd);
  outTriPd->DeepCopy(inTriPd);
  outEdgePd->DeepCopy(inEdgePd);

  int numTriCells  = tmpTriPd->GetNumberOfCells();
  int numTriPts    = tmpTriPd->GetNumberOfPoints();
  int numEdgeCells = tmpEdgePd->GetNumberOfCells();
  int numEdgePts   = tmpEdgePd->GetNumberOfPoints();

  std::vector<int> deletedCell(numTriCells, 0);
  std::vector<int> deletedEdge(numEdgeCells, 0);

  vtkNew(vtkIntArray, tmpEdgeArray);
  tmpEdgeArray->SetNumberOfTuples(numEdgeCells);
  tmpEdgeArray->SetName("TmpInternalIds");
  for (int i=0; i<numEdgeCells; i++)
    tmpEdgeArray->SetTuple1(i, i);
  tmpEdgePd->GetCellData()->AddArray(tmpEdgeArray);

  vtkNew(vtkIntArray, edgeRemoveIterArray);
  edgeRemoveIterArray->SetNumberOfTuples(numEdgeCells);
  edgeRemoveIterArray->FillComponent(0, -1);
  edgeRemoveIterArray->SetName("RemovalIteration");

  vtkNew(vtkIntArray, edgeIsolatedIterArray);
  edgeIsolatedIterArray->SetNumberOfTuples(numEdgeCells);
  edgeIsolatedIterArray->FillComponent(0, -1);
  edgeIsolatedIterArray->SetName("IsolatedIteration");

  vtkNew(vtkIntArray, endIsolatedIterArray);
  endIsolatedIterArray->SetNumberOfTuples(numEdgeCells);
  endIsolatedIterArray->FillComponent(0, -1);
  endIsolatedIterArray->SetName("IsolatedIteration");

  vtkNew(vtkIntArray, removeIterArray);
  removeIterArray->SetNumberOfTuples(numTriCells);
  removeIterArray->FillComponent(0, -1);
  removeIterArray->SetName("RemovalIteration");

  int iter = 0;
  int nDelTris = 0;
  int nDelEdges = 0;
  int nIsolated = 0;

  vtkNew(vtkIdList, ptCellIds);
  vtkNew(vtkIdList, openEdges);
  vtkNew(vtkIdList, cellNeighborIds);
  vtkNew(vtkIdList, pointIds);
  vtkNew(vtkIdList, edgeCell);
  vtkNew(vtkIdList, edgeCellIds);

  std::vector<int> tmpDeletedCells;
  std::vector<int> tmpDeletedEdges;

  vtkIdType npts, *pts;

  int loc;
  int ptId0;
  int ptId1;
  int ptId2;
  int cellId;
  int delEdge;
  int isMedEdge;
  int numDeletedNeighbors     = 0;
  int numNotDeletedNeighbors  = 0;
  int numNotDeletedNeighbors0 = 0;
  int numNotDeletedNeighbors1 = 0;

  // Set up connectivity matrices for tri pd
  std::vector<std::vector<int> > triCellPoints(numTriCells);
  for (int i=0; i<numTriCells; i++)
  {
    tmpTriPd->GetCellPoints(i, npts, pts);
    for (int j=0; j<npts; j++)
      triCellPoints[i].push_back(pts[j]);
  }

  // Set up connectivity matrices for edge pd
  std::vector<std::vector<int> > edgeCellPoints(numEdgeCells);
  std::vector<std::vector<int> > edgePointCells(numEdgePts);
  for (int i=0; i<numEdgeCells; i++)
  {
    tmpEdgePd->GetCellPoints(i, npts, pts);
    for (int j=0; j<npts; j++)
      edgeCellPoints[i].push_back(pts[j]);
  }

  for (int i=0; i<numEdgePts; i++)
  {
    tmpEdgePd->GetPointCells(i, ptCellIds);
    for (int j=0; j<ptCellIds->GetNumberOfIds(); j++)
      edgePointCells[i].push_back(ptCellIds->GetId(j));
  }

  while ( nDelTris > 0 || nDelEdges > 0 || iter == 0 )
  {
    tmpDeletedCells.clear();
    tmpDeletedEdges.clear();
    // --------------------------------------------------------------
    // Do edges before
    nDelEdges = 0;
    for (int i=0; i<numEdgeCells; i++)
    {
      if (!deletedEdge[i])
      {
        npts = edgeCellPoints[i].size();

        if (npts == 2)
        {
          ptId0 = edgeCellPoints[i][0];
          ptId1 = edgeCellPoints[i][1];

          numNotDeletedNeighbors0 = 0;
          numNotDeletedNeighbors1 = 0;
          for (int j=0; j<edgePointCells[ptId0].size(); j++)
          {
            cellId = edgePointCells[ptId0][j];
            if (!deletedEdge[cellId])
              numNotDeletedNeighbors0++;
          }
          for (int j=0; j<edgePointCells[ptId1].size(); j++)
          {
            cellId = edgePointCells[ptId1][j];
            if (!deletedEdge[cellId])
              numNotDeletedNeighbors1++;
          }

          if (numNotDeletedNeighbors0 == 1 ||
              numNotDeletedNeighbors1 == 1)
          {
            delEdge = 1;
            if (dontTouch)
            {
              isMedEdge = tmpEdgePd->GetCellData()->GetArray("MedialEdges")->GetTuple1(i);
              if (isMedEdge == 1)
              {
                delEdge = 0;
              }
            }

            if (delEdge)
            {
              nDelEdges++;
              tmpDeletedEdges.push_back(i);
              edgeRemoveIterArray->SetTuple1(i, iter);
            }
          }

        }
      }
    }
    // --------------------------------------------------------------
    nDelTris = 0;

    for (int i=0; i<numTriCells; i++)
    {
      if (!deletedCell[i])
      {
        npts = triCellPoints[i].size();

        if (npts == 3)
        {
          openEdges->Reset();
          for (int j=0; j<npts; j++)
          {
            ptId0 = triCellPoints[i][j];
            ptId1 = triCellPoints[i][(j+1)%npts];

            tmpTriPd->GetCellEdgeNeighbors(i, ptId0, ptId1, cellNeighborIds);

            if (cellNeighborIds->GetNumberOfIds() == 0)
              openEdges->InsertNextId(j);
            else
            {
              numDeletedNeighbors = 0;
              for (int k=0; k<cellNeighborIds->GetNumberOfIds(); k++)
              {
                if (deletedCell[cellNeighborIds->GetId(k)])
                  numDeletedNeighbors++;
              }
              if (numDeletedNeighbors == cellNeighborIds->GetNumberOfIds())
                openEdges->InsertNextId(j);
            }
          }
          if (openEdges->GetNumberOfIds() == 3)
          {
            nDelTris++;
            tmpDeletedCells.push_back(i);
            removeIterArray->SetTuple1(i, iter);

            // --------------------------------------------------------------
            // Remove on edge pd
            ptId0 = triCellPoints[i][0];
            ptId1 = triCellPoints[i][1];

            pointIds->Reset();
            pointIds->SetNumberOfIds(2);
            pointIds->SetId(0, ptId0);
            pointIds->SetId(1, ptId1);

            tmpEdgePd->GetCellNeighbors(-1, pointIds, edgeCell);

            if (edgeCell->GetNumberOfIds() != 1)
              fprintf(stderr,"NUMBER OF CELLS IS NOT 1, IT IS %d\n", edgeCell->GetNumberOfIds());
            else
            {
              delEdge = 1;
              if (dontTouch)
              {
                isMedEdge = tmpEdgePd->GetCellData()->GetArray("MedialEdges")->GetTuple1(edgeCell->GetId(0));
                if (isMedEdge == 1)
                {
                  delEdge = 0;
                }
              }

              if (delEdge)
              {
                nDelEdges++;
                tmpDeletedEdges.push_back(edgeCell->GetId(0));
                edgeRemoveIterArray->SetTuple1(edgeCell->GetId(0), iter);
              }
            }

            // --------------------------------------------------------------
          }
          else if (openEdges->GetNumberOfIds() == 2)
          {
            nDelTris++;
            loc;
            for (int j=0; j<npts; j++)
            {
              if (j != openEdges->GetId(0) && j != openEdges->GetId(1))
                loc = j;
            }

            ptId0 = triCellPoints[i][loc];
            ptId1 = triCellPoints[i][(loc+1)%npts];
            ptId2 = triCellPoints[i][(loc+2)%npts];

            tmpDeletedCells.push_back(i);
            removeIterArray->SetTuple1(i, iter);

            // --------------------------------------------------------------
            // Remove on edge pd
            pointIds->Reset();
            pointIds->SetNumberOfIds(2);
            pointIds->SetId(0, ptId0);
            pointIds->SetId(1, ptId2);

            tmpEdgePd->GetCellNeighbors(-1, pointIds, edgeCell);

            if (edgeCell->GetNumberOfIds() != 1)
              fprintf(stderr,"NUMBER OF CELLS IS NOT 1, IT IS %d\n", edgeCell->GetNumberOfIds());
            else
            {
              delEdge = 1;
              if (dontTouch)
              {
                isMedEdge = tmpEdgePd->GetCellData()->GetArray("MedialEdges")->GetTuple1(edgeCell->GetId(0));
                if (isMedEdge == 1)
                {
                  delEdge = 0;
                }
              }

              if (delEdge)
              {
                nDelEdges++;
                tmpDeletedEdges.push_back(edgeCell->GetId(0));
                edgeRemoveIterArray->SetTuple1(edgeCell->GetId(0), iter);
              }
            }


            // --------------------------------------------------------------
          }
          else if (openEdges->GetNumberOfIds() == 1)
          {
            nDelTris++;
            loc = openEdges->GetId(0);

            ptId0 = triCellPoints[i][loc];
            ptId1 = triCellPoints[i][(loc+1)%npts];
            ptId2 = triCellPoints[i][(loc+2)%npts];

            tmpDeletedCells.push_back(i);
            removeIterArray->SetTuple1(i, iter);

            // --------------------------------------------------------------
            // Remove on edge pd
            pointIds->Reset();
            pointIds->SetNumberOfIds(2);
            pointIds->SetId(0, ptId0);
            pointIds->SetId(1, ptId1);

            tmpEdgePd->GetCellNeighbors(-1, pointIds, edgeCell);

            if (edgeCell->GetNumberOfIds() != 1)
              fprintf(stderr,"NUMBER OF CELLS IS NOT 1, IT IS %d\n", edgeCell->GetNumberOfIds());
            else
            {
              delEdge = 1;
              if (dontTouch)
              {
                isMedEdge = tmpEdgePd->GetCellData()->GetArray("MedialEdges")->GetTuple1(edgeCell->GetId(0));
                if (isMedEdge == 1)
                {
                  delEdge = 0;
                }
              }

              if (delEdge)
              {
                nDelEdges++;
                tmpDeletedEdges.push_back(edgeCell->GetId(0));
                edgeRemoveIterArray->SetTuple1(edgeCell->GetId(0), iter);
              }
            }

            // --------------------------------------------------------------
          }
        }
      }
    }
    fprintf(stdout,"ITER %d, TRIS REMOVED: %d, EDGES REMOVED: %d\n", iter, nDelTris, nDelEdges);

    for (int i=0; i<tmpDeletedCells.size(); i++)
      deletedCell[tmpDeletedCells[i]] = 1;
    for (int i=0; i<tmpDeletedEdges.size(); i++)
      deletedEdge[tmpDeletedEdges[i]] = 1;

    // --------------------------------------------------------------
    // Now add to edge isolated list
    if (nIsolated != numEdgeCells)
    {
      for (int i=0; i<numEdgeCells; i++)
      {
        int currVal = edgeIsolatedIterArray->GetTuple1(i);

        if (currVal == -1)
        {
          npts = edgeCellPoints[i].size();

          if (npts == 2)
          {
            pointIds->Reset();
            pointIds->SetNumberOfIds(2);
            pointIds->SetId(0, edgeCellPoints[i][0]);
            pointIds->SetId(1, edgeCellPoints[i][1]);

            tmpTriPd->GetCellNeighbors(-1, pointIds, edgeCellIds);

            numDeletedNeighbors = 0;
            for (int j=0; j<edgeCellIds->GetNumberOfIds(); j++)
            {
              if (deletedCell[edgeCellIds->GetId(j)])
                numDeletedNeighbors++;
            }

            if (numDeletedNeighbors == edgeCellIds->GetNumberOfIds())
            {
              endIsolatedIterArray->SetTuple1(i, iter);
              edgeIsolatedIterArray->SetTuple1(i, iter);
              nIsolated++;
            }
          }
        }
      }
    }
    // --------------------------------------------------------------

    iter++;
  }

  outTriPd->GetCellData()->AddArray(removeIterArray);

  outEdgePd->GetCellData()->AddArray(edgeRemoveIterArray);
  outEdgePd->GetCellData()->AddArray(endIsolatedIterArray);

  vtkNew(vtkIntArray, mAbsArray);
  mAbsArray->SetNumberOfTuples(numEdgeCells);
  mAbsArray->FillComponent(0, -1);
  mAbsArray->SetName("MAbs");

  vtkNew(vtkDoubleArray, mRelArray);
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

  return SV_OK;
}

int vtkSVCenterlines::RecursiveGetPolylines(vtkPolyData *pd,
                                            std::vector<std::vector<int> > connectedEdgePts,
                                            int startVertex, std::vector<int> &pointUsed,
                                            std::vector<std::vector<int> > &allEdges,
                                            std::vector<int> &thisEdge)
{
	//int startVertex = 0;
	int i, j, index, testIndex0, testIndex1, firstVertex, prevVertex, secondVertex, countTotal;
	double tempDouble[3], tempX[3], tempY[3], tempZ[3], tempXPre[3];
  vtkNew(vtkIdList, pointCellIds);
  vtkNew(vtkIdList, edgePointIds);
	int stopCriteria;

  int numPts = pd->GetNumberOfPoints();

	stopCriteria = 0;
  firstVertex = startVertex;

	while (!stopCriteria)
	{

		prevVertex = -1;
		secondVertex = -1;

		if (connectedEdgePts[firstVertex].size() == 1)
		{
			index = connectedEdgePts[firstVertex][0];
			if (pointUsed[index] == 0)
			{
				secondVertex = index;
			}
			else
			{
				prevVertex = index;
			}
		}
		else if (connectedEdgePts[firstVertex].size() == 2)
		{
      testIndex0 = connectedEdgePts[firstVertex][0];
      testIndex1 = connectedEdgePts[firstVertex][1];

      if (pointUsed[testIndex0]  && pointUsed[testIndex1])
      {
        pointUsed[firstVertex] = 1;
        thisEdge.push_back(firstVertex);

        // Have to add the last point
        int found0 = 0, found1 = 0;
        for (int i=0; i<thisEdge.size(); i++)
        {
          if (thisEdge[i] == testIndex0)
            found0 = 1;
          if (thisEdge[i] == testIndex1)
            found1 = 1;
        }
        if (!found0)
          thisEdge.push_back(testIndex0);
        else if(!found1)
          thisEdge.push_back(testIndex1);
        else
        {
          fprintf(stderr,"BOTH ALREADY IN EDGE!\n");
          return SV_ERROR;
        }
        allEdges.push_back(thisEdge);
        //fprintf(stdout,"ALSO REACHED AN END %d\n", firstVertex);
        return SV_OK;
      }

			for (i = 0; i < connectedEdgePts[firstVertex].size(); i++)
			{
				index = connectedEdgePts[firstVertex][i];
				if (pointUsed[index] == 0)
				{
					secondVertex = index;
				}
				else
				{
					prevVertex = index;
				}
			}
		}
    else if (connectedEdgePts[firstVertex].size() > 2)
    {
      pointUsed[firstVertex] = 1;
      thisEdge.push_back(firstVertex);
      allEdges.push_back(thisEdge);
			for (i = 0; i < connectedEdgePts[firstVertex].size(); i++)
			{
				index = connectedEdgePts[firstVertex][i];
        if (pointUsed[index] == 0)
        {
          std::vector<int> newEdge;
          newEdge.push_back(firstVertex);
          this->RecursiveGetPolylines(pd, connectedEdgePts, index, pointUsed, allEdges, newEdge);
        }
        else if (connectedEdgePts[index].size() > 2)
        {
          std::vector<int> potEdge;
          potEdge.push_back(firstVertex);
          potEdge.push_back(index);
          int alreadyEdge = 0;
          for (int j=0; j<allEdges.size(); j++)
          {
            if ((allEdges[j][0] == potEdge[0] && allEdges[j][allEdges[j].size()-1] == potEdge[1]) ||
                (allEdges[j][0] == potEdge[1] && allEdges[j][allEdges[j].size()-1] == potEdge[0]))
            {
              alreadyEdge = 1;
            }
          }
          if (!alreadyEdge)
            allEdges.push_back(potEdge);
        }
      }

      return SV_OK;
    }
    else
    {
      fprintf(stderr,"Somehow point is connected to nothing\n");
      return SV_ERROR;
    }

		pointUsed[firstVertex] = 1;
    thisEdge.push_back(firstVertex);

		if (connectedEdgePts[firstVertex].size() == 1)
		{

			index = connectedEdgePts[firstVertex][0];
			if (pointUsed[index] == 1)
			{
				stopCriteria = 1;
        allEdges.push_back(thisEdge);
        return SV_OK;
			}

		}

		firstVertex = secondVertex;
	}

  return SV_OK;
}

int vtkSVCenterlines::RecursiveGetFullCenterlines(std::vector<std::vector<int> > allEdges,
                                                  std::vector<std::vector<int> > &fullCenterlineEdges,
                                                  int thisEdge, int front, int back)
{
  int stillRecursing = 0;
  for (int i=0; i<allEdges.size(); i++)
  {
    if (i == thisEdge)
      continue;

    int edgeSize = allEdges[i].size();
    int edgeId0 = allEdges[i][0];
    int edgeIdN = allEdges[i][edgeSize-1];

    std::vector<std::vector<int> > newCenterlineEdges;
    if (edgeId0 == back)
    {
      int newFront = edgeId0;
      int newBack  = edgeIdN;

      this->RecursiveGetFullCenterlines(allEdges, newCenterlineEdges, i, newFront, newBack);

    }

    if (edgeIdN == back)
    {
      int newFront = edgeIdN;
      int newBack  = edgeId0;

      this->RecursiveGetFullCenterlines(allEdges, newCenterlineEdges, i, newFront, newBack);
    }

    if (newCenterlineEdges.size() > 0)
    {
      for (int j=0; j<newCenterlineEdges.size(); j++)
      {
        std::vector<int> newEdge;
        newEdge.push_back(thisEdge);
        for (int k=0; k<newCenterlineEdges[j].size(); k++)
          newEdge.push_back(newCenterlineEdges[j][k]);

        fullCenterlineEdges.push_back(newEdge);
      }
      stillRecursing = 1;
    }
  }

  if (!stillRecursing)
  {
    std::vector<int> newEdge;
    newEdge.push_back(thisEdge);
    fullCenterlineEdges.push_back(newEdge);
  }

  return SV_OK;
}

int vtkSVCenterlines::RemoveMarkedCells(vtkPolyData *pd,
                                        std::vector<std::vector<int> > allEdges,
                                        std::vector<int> needToDelete,
                                        std::vector<int> &isDeleted,
                                        vtkIdList *allEndIds,
                                        std::vector<int> &nodeCount)
{
  int edgeSize, pointId0, pointId1, nodeId0, nodeIdN;
  vtkNew(vtkIdList, edgePointIds);
  vtkNew(vtkIdList, pointsEdgeId);
  for (int i=0; i<needToDelete.size(); i++)
  {
    if (isDeleted[i])
    {
      continue;
    }

    if (needToDelete[i] == 1)
    {
      edgeSize = allEdges[i].size();
      for (int j=0; j<edgeSize-1; j++)
      {
        pointId0 = allEdges[i][j];
        pointId1 = allEdges[i][j+1];
        edgePointIds->Reset();
        edgePointIds->InsertNextId(pointId0);
        edgePointIds->InsertNextId(pointId1);

        pd->GetCellNeighbors(-1, edgePointIds, pointsEdgeId);
        if (pointsEdgeId->GetNumberOfIds() == 1)
        {
          pd->DeleteCell(pointsEdgeId->GetId(0));
        }
      }

      nodeId0 = allEndIds->IsId(allEdges[i][0]);
      nodeIdN = allEndIds->IsId(allEdges[i][edgeSize-1]);
      nodeCount[nodeId0]--;
      nodeCount[nodeIdN]--;
      isDeleted[i] = 1;
    }
  }

  return SV_OK;
}

int vtkSVCenterlines::GetLinesEndPoints(vtkPolyData *pd,
                                        vtkIdList *endPointIds,
                                        vtkPoints *endPoints,
                                        std::vector<std::vector<int> > &connectedEdgePts,
                                        int &firstVertex)
{
	firstVertex = -1;
  connectedEdgePts.clear();
  connectedEdgePts.resize(pd->GetNumberOfPoints());

  endPointIds->Reset();
  endPoints->Reset();

  vtkIdType npts, *pts;
  int numEndPoints = 0;
  double maxRadiusValue = -1.0;
  vtkNew(vtkIdList, pointCellIds);

  for (int i=0; i<pd->GetNumberOfPoints(); i++)
  {
    pd->GetPointCells(i, pointCellIds);
    if (pointCellIds->GetNumberOfIds() == 1)
    {
      numEndPoints++;
      double radiusValue = pd->GetPointData()->
        GetArray(this->RadiusArrayName)->GetTuple1(i);

      if (radiusValue > maxRadiusValue)
      {
        maxRadiusValue = radiusValue;
        firstVertex = i;
      }

      endPointIds->InsertNextId(i);
      endPoints->InsertNextPoint(pd->GetPoint(i));
    }

    std::vector<int> connectedPts;
    for (int j=0; j<pointCellIds->GetNumberOfIds(); j++)
    {
      pd->GetCellPoints(pointCellIds->GetId(j), npts, pts);

      for (int k=0; k<npts; k++)
      {
        if (pts[k] != i)
          connectedPts.push_back(pts[k]);
      }
    }
    connectedEdgePts[i] = connectedPts;
  }

  return SV_OK;
}
