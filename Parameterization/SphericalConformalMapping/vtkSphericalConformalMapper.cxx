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

/** @file vtkSphericalConformalMapper.cxx
 *  @brief This implements the vtkSphericalConformalMapper filter as a class
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkSphericalConformalMapper.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCenterOfMass.h"
#include "vtkCleanPolyData.h"
#include "vtkConnectivityFilter.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkFeatureEdges.h"
#include "vtkFloatArray.h"
#include "vtkGradientFilter.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPlacePointsOnS2.h"
#include "vtkPointData.h"
#include "vtkPointDataToCellData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmartPointer.h"
#include "vtkTransform.h"
#include "vtkTriangle.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLPolyDataWriter.h"

#include <iostream>
#include <sstream>
#include <cmath>

//---------------------------------------------------------------------------
vtkCxxRevisionMacro(vtkSphericalConformalMapper, "$Revision: 0.0 $");
vtkStandardNewMacro(vtkSphericalConformalMapper);


//---------------------------------------------------------------------------
vtkSphericalConformalMapper::vtkSphericalConformalMapper()
{
  this->Verbose                 = 1;
  this->InitialTimeStep         = 0.01;
  this->TimeStep                = 0.1;
  this->TutteEnergyCriterion    = 0.00001;
  this->HarmonicEnergyCriterion = 0.0000001;
  this->MaxNumIterations        = 1000;
  this->NumBoundaries           = 0;
  this->CGUpdateMethod          = CG_NONE;

  this->BoundaryType            = CLOSED;
  this->BoundaryConstraintType  = 0;
  for (int i=0; i<2; i++)
  {
    this->BoundaryStart[i] = -1;
  }
  this->FirstLoopPts      = NULL;
  this->SecondLoopPts     = NULL;
  this->FirstLoopHelper   = NULL;
  this->SecondLoopHelper  = NULL;
  for (int i=0; i<2; i++)
  {
    this->CubeStart[i] = -1;
  }

  this->InitialPd      = vtkPolyData::New();
  this->EdgeTable      = vtkEdgeTable::New();
  this->EdgeWeights    = vtkFloatArray::New();
  this->PrevDescent    = vtkFloatArray::New();
  this->CurrDescent    = vtkFloatArray::New();
  this->ConjugateDir   = vtkFloatArray::New();
  this->EdgeNeighbors  = vtkIntArray::New();
  this->IsBoundary     = vtkIntArray::New();
  this->HarmonicMap[0] = vtkPolyData::New();
  this->HarmonicMap[1] = vtkPolyData::New();
  this->Boundaries     = vtkPolyData::New();
  this->SetObjectXAxis(1.0, 0.0, 0.0);
  this->SetObjectZAxis(0.0, 0.0, 1.0);

  this->IterOutputFilename = NULL;
  this->NumSaveIterations  = 100;
  this->SaveIter           = 0;
}

//---------------------------------------------------------------------------
vtkSphericalConformalMapper::~vtkSphericalConformalMapper()
{
  if (this->InitialPd != NULL)
  {
    InitialPd->Delete();
  }
  if (this->EdgeTable != NULL)
  {
    this->EdgeTable->Delete();
  }
  if (this->EdgeWeights != NULL)
  {
    this->EdgeWeights->Delete();
  }
  if (this->PrevDescent != NULL)
  {
    this->PrevDescent->Delete();
  }
  if (this->CurrDescent != NULL)
  {
    this->CurrDescent->Delete();
  }
  if (this->ConjugateDir != NULL)
  {
    this->ConjugateDir->Delete();
  }
  if (this->EdgeNeighbors != NULL)
  {
    this->EdgeNeighbors->Delete();
  }
  if (this->IsBoundary != NULL)
  {
    this->IsBoundary->Delete();
  }
  for (int i=0; i<2; i++)
  {
    if (this->HarmonicMap[i] != NULL)
    {
      this->HarmonicMap[i]->Delete();
    }
  }
  if (this->Boundaries != NULL)
  {
    this->Boundaries->Delete();
  }
}

//---------------------------------------------------------------------------
void vtkSphericalConformalMapper::PrintSelf(ostream& os, vtkIndent indent)
{
}

// Generate Separated Surfaces with Region ID Numbers
//---------------------------------------------------------------------------
int vtkSphericalConformalMapper::RequestData(
                                 vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector)
{
  // get the input and output
  vtkPolyData *input = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  //Copy the input to operate on
  this->InitialPd->DeepCopy(input);

  vtkIdType numPolys = this->InitialPd->GetNumberOfPolys();
  //Check the input to make sure it is there
  if (numPolys < 1)
  {
    vtkDebugMacro("No input!");
    return 0;
  }

  //Check the input to make sure it is manifold and a triangulated surface
  if (this->CleanAndCheckSurface(this->InitialPd) != 1)
  {
    vtkErrorMacro("Error when checking input surface");
    return 0;
  }
  //Get the number of Polys for scalar  allocation
  numPolys = this->InitialPd->GetNumberOfPolys();
  vtkIdType numPts = this->InitialPd->GetNumberOfPoints();

  //Create the edge table for the input surface
  this->CreateEdgeTable(this->InitialPd, this->EdgeTable, this->EdgeWeights,
                        this->EdgeNeighbors, this->IsBoundary);

  if (this->PerformMapping() != 1)
  {
    vtkErrorMacro("Error while doing CG Solve");
    return 0;
  }

  output->DeepCopy(this->HarmonicMap[HARMONIC]);
  output->GetPointData()->PassData(input->GetPointData());
  output->GetCellData()->PassData(input->GetCellData());
  if (this->PDCheckArrayName(output, 0 ,"Normals") == 1)
  {
    output->GetPointData()->RemoveArray("Normals");
  }
  if (this->PDCheckArrayName(output, 1,"cellNormals") == 1)
  {
    output->GetCellData()->RemoveArray("cellNormals");
  }
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::CleanAndCheckSurface(vtkPolyData *pd)
{
  vtkSmartPointer<vtkPolyDataNormals> normaler =
    vtkSmartPointer<vtkPolyDataNormals>::New();
  normaler->SetInputData(pd);
  normaler->AutoOrientNormalsOn();
  normaler->ComputeCellNormalsOff();
  normaler->SplittingOff();
  normaler->Update();

  pd->DeepCopy(normaler->GetOutput());
  pd->BuildLinks();

  int numPts = pd->GetNumberOfPoints();
  int numPolys = pd->GetNumberOfCells();

  for (int i=0; i<numPolys; i++)
  {
    vtkIdType npts, *pts;
    pd->GetCellPoints(i, npts, pts);
    if (npts != 3)
    {
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
        return 0;
      }
      //if (edgeNeighbor->GetNumberOfIds() < 1)
      //{
      //  return 0;
      //}
    }
  }
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::CreateEdgeTable(vtkPolyData *pd,
                                                  vtkEdgeTable *edgeTable,
                                                  vtkFloatArray *edgeWeights,
                                                  vtkIntArray *edgeNeighbors,
                                                  vtkIntArray *isBoundary)
{
  int numPts = pd->GetNumberOfPoints();
  int numTris = pd->GetNumberOfCells();

  edgeTable->InitEdgeInsertion(numPts, 1);
  isBoundary->SetNumberOfValues(numPts);
  for (int i=0; i<numPts; i++)
  {
    isBoundary->InsertValue(i, 0);
  }
  for (int i=0; i<numTris; i++)
  {
    vtkIdType npts, *pts;
    pd->GetCellPoints(i, npts, pts);
    double pt0[3], pt1[3], pt2[3];
    pd->GetPoint(pts[0], pt0);
    pd->GetPoint(pts[1], pt1);
    pd->GetPoint(pts[2], pt2);

    //Make sure the points are ordered counter-clockwise!
    //double area = 0.0;
    //vtkSphericalConformalMapper::ComputeArea(pt0, pt1, pt2, area);
    //if (area < 0)
    //{
    //  double tmpPoint = pts[0];
    //  pts[0] = pts[1];
    //  pts[1] = tmpPoint;
    //  pd->ReplaceCell(i, npts, pts);
    //}

    //Insert edge into table
    for (int j=0; j<npts; j++)
    {
      vtkIdType p0 = pts[j];
      vtkIdType p1 = pts[(j+1)%npts];
      vtkSmartPointer<vtkIdList> neighborCellIds =
        vtkSmartPointer<vtkIdList>::New();
      pd->GetCellEdgeNeighbors(i, p0, p1, neighborCellIds);
      vtkIdType neighborCellId = 0;
      if (neighborCellIds->GetNumberOfIds() > 0)
      {
        neighborCellId = neighborCellIds->GetId(0);
      }
      else
      {
        neighborCellId = -1;
        isBoundary->InsertValue(p0, 1);
        isBoundary->InsertValue(p1, 1);
      }
      vtkIdType checkEdge = edgeTable->IsEdge(p0, p1);
      if (checkEdge == -1)
      {
        //Compute Edge Weight
        double weight = 0.0;
        vtkSphericalConformalMapper::ComputeEdgeWeight(pd, i, neighborCellId,
                                                        p0, p1, weight);
        vtkIdType edgeId = edgeTable->InsertEdge(p0, p1);
        edgeWeights->InsertValue(edgeId, weight);
        edgeNeighbors->InsertValue(edgeId, neighborCellId);
        //if (weight < 0)
        //{
        //  fprintf(stdout,"Negative weight!: %.4f\n",weight);
        //  fprintf(stdout,"Cells: %.d and %.lld\n",i,neighborCellId);
        //}
      }
    }
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::ComputeEdgeWeight(vtkPolyData *pd,
                                                    vtkIdType cellId,
                                                    vtkIdType neighborCellId,
                                                    vtkIdType p0,
                                                    vtkIdType p1,
                                                    double &weight)
{
  //Add the edge weights based on the angle of the edge
  vtkIdType cellIds[2];
  cellIds[0] = cellId;
  cellIds[1] = neighborCellId;
  weight = 0.0;
  double v0[3]; double v1[3]; double v2[3];
  pd->GetPoint(p0, v0);
  pd->GetPoint(p1, v1);
  for (int i=0; i<2; i++)
  {
    vtkIdType npts, *pts;
    if (cellIds[i] != -1)
    {
      pd->GetCellPoints(cellIds[i], npts, pts);
      for (int k=0; k<npts; k++)
      {
	if (pts[k] != p0 && pts[k] != p1)
	{
	  pd->GetPoint(pts[k], v2);
	  double angle = 0.0;
	  vtkSphericalConformalMapper::GetEdgeCotangentAngle(v0, v1, v2, angle);

	  weight += 0.5*angle;
	}
      }
    }
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::GetEdgeCotangentAngle(double pt0[], double pt1[],
                                                        double pt2[], double &angle)
{
  double area = 0.0;
  vtkSphericalConformalMapper::ComputeArea(pt0, pt1, pt2, area);
  if (area < 0)
  {
    double tmpPoint[3];
    for (int i=0; i<3; i++)
    {
      tmpPoint[i] = pt0[i];
      pt0[i] = pt1[i];
      pt1[i] = tmpPoint[i];
    }
  }
  double vec0[3];
  double vec1[3];
  for (int i=0; i<3; i++)
  {
    vec0[i] = pt0[i] - pt2[i];
    vec1[i] = pt1[i] - pt2[i];
  }
  double numerator = vtkMath::Dot(vec0, vec1);
  double cross[3];
  vtkMath::Cross(vec0, vec1, cross);
  double denominator = vtkMath::Norm(cross);
  angle = numerator/denominator;

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::ConvertFieldToPolyData(vtkPolyData *inPd, std::string fieldName, vtkPolyData *outPd)
{
  int numCells = inPd->GetNumberOfCells();
  int numPts   = inPd->GetNumberOfPoints();
  vtkSmartPointer<vtkPoints> fieldPts = vtkSmartPointer<vtkPoints>::New();
  fieldPts->SetNumberOfPoints(numPts);
  vtkSmartPointer<vtkCellArray> fieldCells = vtkSmartPointer<vtkCellArray>::New();
  fieldCells = inPd->GetPolys();

  vtkFloatArray *fieldArray;
  fieldArray = vtkFloatArray::SafeDownCast(
    inPd->GetPointData()->GetArray(fieldName.c_str()));
  for (int i=0; i<numCells; i++)
  {
    vtkIdType npts, *pts;
    inPd->GetCellPoints(i, npts, pts);
    for (int j=0; j<npts; j++)
    {
      double pt[3];
      fieldArray->GetTuple(pts[j], pt);
      fieldPts->SetPoint(pts[j], pt);
    }
  }

  outPd->SetPolys(fieldCells);
  outPd->SetPoints(fieldPts);
  outPd->BuildLinks();

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::ProjectOntoUnitSphere(vtkPolyData *inPd,
                                                         vtkPolyData *outPd)
{
  int numPts = inPd->GetNumberOfPoints();
  double massCenter[3];

  vtkSphericalConformalMapper::ComputeMassCenter(inPd, massCenter);

  vtkSmartPointer<vtkPoints> newPts = vtkSmartPointer<vtkPoints>::New();
  newPts->SetNumberOfPoints(numPts);
  vtkSmartPointer<vtkCellArray> newCells = vtkSmartPointer<vtkCellArray>::New();
  newCells = inPd->GetPolys();

  for (int i=0; i<numPts; i++)
  {
    double pt[3];
    double newpt[3];
    inPd->GetPoint(i, pt);
    for (int j=0; j<3; j++)
    {
      newpt[j] = pt[j] - massCenter[j];
    }
    vtkMath::Normalize(newpt);
    newPts->SetPoint(i, newpt);
  }

  outPd->SetPolys(newCells);
  outPd->SetPoints(newPts);
  outPd->BuildLinks();
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::ComputeNormals(vtkPolyData *pd)
{
  vtkSmartPointer<vtkPolyDataNormals> normaler =
    vtkSmartPointer<vtkPolyDataNormals>::New();
  normaler->SetInputData(pd);
  normaler->AutoOrientNormalsOn();
  normaler->SplittingOff();
  normaler->Update();

  pd->DeepCopy(normaler->GetOutput());
  pd->BuildLinks();

  return 1;
}


//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::ComputeArea(double pt0[], double pt1[],
                                              double pt2[], double &area)
{
  area = 0.0;
  area += (pt0[0]*pt1[1])-(pt1[0]*pt0[1]);
  area += (pt1[0]*pt2[1])-(pt2[0]*pt1[1]);
  area += (pt2[0]*pt0[1])-(pt0[0]*pt2[1]);
  area *= 0.5;

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::PerformMapping()
{
  fprintf(stdout,"CG Update Method: %d\n",this->CGUpdateMethod);

  //Compute Gauss Map (i.e. get normals off of mesh)
  //TODO Check normals exists!
  //this->ConvertFieldToPolyData(this->InitialPd, "Normals", this->HarmonicMap[0]);
  vtkSmartPointer<vtkPlacePointsOnS2> initialSpot =
    vtkSmartPointer<vtkPlacePointsOnS2>::New();
  initialSpot->SetInputData(this->InitialPd);
  initialSpot->SetUseCustomAxisAlign(1);
  initialSpot->SetXAxis(this->ObjectXAxis);
  initialSpot->SetZAxis(this->ObjectZAxis);
  initialSpot->Update();
  this->HarmonicMap[0]->DeepCopy(initialSpot->GetOutput());

  if(this->SetBoundaries() != 1)
  {
    vtkErrorMacro("Error when setting boundaries");
    return 0;
  }
  //std::string filename = "/Users/adamupdegrove/Desktop/tmp/S2Placed.vtp";
  //vtkSmartPointer<vtkXMLPolyDataWriter> writer =
  //  vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  //writer->SetInputData(this->HarmonicMap[0]);
  //writer->SetFileName(filename.c_str());
  //writer->Write();

  if (this->InitiateCGArrays() != 1)
  {
    vtkErrorMacro("Error when initiating arrays");
    return 0;
  }

  //Run the Tutte Energy Step
  if (this->SphericalTutteMapping() != 1)
  {
    vtkErrorMacro("Error when computing the tutte map");
    return 0;
  }
  //Compute Initial Tutte energy
  this->HarmonicMap[1]->DeepCopy(this->HarmonicMap[0]);

  //Run the Spherical Conformal Mapping Step
  if (this->SphericalConformalMapper() != 1)
  {
    vtkErrorMacro("Error when computing the conformal map");
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
int vtkSphericalConformalMapper::InitiateCGArrays()
{
  int numPts = this->InitialPd->GetNumberOfPoints();

  this->PrevDescent->SetNumberOfComponents(3);
  this->PrevDescent->Allocate(numPts, 10000);
  this->PrevDescent->SetNumberOfTuples(numPts);

  this->CurrDescent->SetNumberOfComponents(3);
  this->CurrDescent->Allocate(numPts, 10000);
  this->CurrDescent->SetNumberOfTuples(numPts);

  this->ConjugateDir->SetNumberOfComponents(3);
  this->ConjugateDir->Allocate(numPts, 10000);
  this->ConjugateDir->SetNumberOfTuples(numPts);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::SetBoundaries()
{
  if (this->FindBoundaries() != 1)
  {
    vtkErrorMacro("Could not find boundaries");
    return 0;
  }

  if (this->NumBoundaries != 0)
  {
    //TODO: Clean up correctly on errors!!!
    vtkPolyData **boundaryLoops = new vtkPolyData*[this->NumBoundaries];
    for (int i=0; i<this->NumBoundaries; i++)
    {
      boundaryLoops[i] = vtkPolyData::New();
      boundaryLoops[i]->SetPoints(this->Boundaries->GetPoints());
      boundaryLoops[i]->GetPointData()->PassData(this->Boundaries->GetPointData());
    }
    vtkIntArray *pointIds =
      vtkIntArray::SafeDownCast(this->Boundaries->GetPointData()->GetArray("PointIds"));

    int numBoundStarts = 0;
    int boundaryStart[2];
    for (int i=0; i<2; i++)
    {
      if (this->BoundaryStart[i] != -1)
      {
        numBoundStarts++;
        boundaryStart[i] = pointIds->LookupValue(this->BoundaryStart[i]);
      }
    }
    if (numBoundStarts != this->NumBoundaries)
    {
      vtkErrorMacro("The number of boundaries does not match with the information you have provided for the boundary");
      for (int j=0; j<this->NumBoundaries; j++)
      {
        boundaryLoops[j]->Delete();
      }
      delete [] boundaryLoops;
      return 0;
    }
    if (this->SeparateLoops(this->Boundaries, boundaryLoops, numBoundStarts,
                            this->ObjectXAxis, this->ObjectZAxis, boundaryStart) != 1)
    {
      vtkErrorMacro("No separate");
      for (int j=0; j<this->NumBoundaries; j++)
      {
        boundaryLoops[j]->Delete();
      }
      delete [] boundaryLoops;
      return 0;
    }
    vtkIntArray *boundaryLoopPts[2];
    vtkIntArray *boundaryLoopHelper[2];
    boundaryLoopPts[0] = this->FirstLoopPts;
    boundaryLoopPts[1] = this->SecondLoopPts;
    boundaryLoopHelper[0] = this->FirstLoopHelper;
    boundaryLoopHelper[1] = this->SecondLoopHelper;
    //double rightval = 1 + std::sqrt(2.0)/2.0;
    double rightval = 1 + 0.58;
    double radius = std::sqrt(rightval/(2 - rightval));
    for (int i=0; i<this->NumBoundaries; i++)
    {
      int numLoopPts = boundaryLoopPts[i]->GetNumberOfTuples();
      double *lengths = new double[numLoopPts];
      if (this->CalculateSquareEdgeLengths(boundaryLoops[i], boundaryLoopPts[i], lengths) != 1)
      {
        vtkErrorMacro("Didn't work");
        delete [] lengths;
        for (int j=0; j<this->NumBoundaries; j++)
        {
          boundaryLoops[j]->Delete();
        }
        delete [] boundaryLoops;
        return 0;
      }
      double cubeStart[3];
      this->GetCubeStartPoint(this->CubeStart[i], cubeStart);
      int ok = 0;
      if (this->BoundaryConstraintType == 0)
      {
        ok = this->SetCircleBoundary(boundaryLoops[i], boundaryLoopPts[i], boundaryLoopHelper[i], cubeStart, lengths, radius);
      }
      else
      {
        ok = this->SetCubeBoundary(boundaryLoops[i], boundaryLoopPts[i], boundaryLoopHelper[i], cubeStart, lengths);
      }
      if (!ok)
      {
        delete [] lengths;
        for (int j=0; j<this->NumBoundaries; j++)
        {
          boundaryLoops[j]->Delete();
        }
        delete [] boundaryLoops;
        return 0;
      }

      //double length = 0.0;
      //if (this->CalculateCircleLength(boundaryLoops[i], length) != 1)
      //{
      //  vtkErrorMacro("Didn't work");
      //  for (int j=0; j<this->NumBoundaries; j++)
      //  {
      //    boundaryLoops[j]->Delete();
      //  }
      //  delete [] boundaryLoops;
      //  return 0;
      //}
      //if (this->SetLoopOnUnitCircle(boundaryLoops[i], length, radius) != 1)
      //{
      //  vtkErrorMacro("Didn't work");
      //  for (int j=0; j<this->NumBoundaries; j++)
      //  {
      //    boundaryLoops[j]->Delete();
      //  }
      //  delete [] boundaryLoops;
      //  return 0;
      //}
      //delete [] lengths;

      radius = std::sqrt((2 - rightval)/rightval);
    }

    for (int i=0; i<this->NumBoundaries; i++)
    {
      boundaryLoops[i]->Delete();
    }
    delete [] boundaryLoops;
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::FindBoundaries()
{
  int numPts = this->InitialPd->GetNumberOfPoints();
  vtkSmartPointer<vtkIntArray> pointIds =
    vtkSmartPointer<vtkIntArray>::New();
  pointIds->SetNumberOfValues(numPts);
  pointIds->Allocate(numPts, 10000);
  pointIds->SetName("PointIds");
  for (int i=0; i<numPts; i++)
  {
    pointIds->SetValue(i, i);
  }
  this->InitialPd->GetPointData()->AddArray(pointIds);
  vtkSmartPointer<vtkFeatureEdges> finder =
    vtkSmartPointer<vtkFeatureEdges>::New();
  finder->SetInputData(this->InitialPd);
  finder->FeatureEdgesOff();
  finder->NonManifoldEdgesOff();
  finder->BoundaryEdgesOn();
  finder->Update();

  vtkSmartPointer<vtkConnectivityFilter> connector =
    vtkSmartPointer<vtkConnectivityFilter>::New();
  connector->SetInputData(finder->GetOutput());
  connector->SetExtractionMode(VTK_EXTRACT_ALL_REGIONS);
  connector->ColorRegionsOn();
  connector->Update();

  vtkSmartPointer<vtkDataSetSurfaceFilter> surfacer =
    vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  surfacer->SetInputData(connector->GetOutput());
  surfacer->Update();

  this->NumBoundaries = connector->GetNumberOfExtractedRegions();
  this->Boundaries->DeepCopy(surfacer->GetOutput());

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::FirstStep(int map)
{
  int numPts = this->InitialPd->GetNumberOfPoints();
  vtkSmartPointer<vtkFloatArray> laplacian =
    vtkSmartPointer<vtkFloatArray>::New();
  laplacian->SetNumberOfComponents(3);
  laplacian->Allocate(numPts, 10000);
  laplacian->SetNumberOfTuples(numPts);

  if (this->ComputeMeshLaplacian(this->HarmonicMap[map], this->EdgeTable,
                                 this->EdgeWeights, this->EdgeNeighbors,
                                 laplacian, map) != 1)
  {
    vtkErrorMacro("Error when computing laplacian");
    return 0;
  }
  if (this->UpdateMap(laplacian, map, CG_NONE) != 1)
  {
    vtkErrorMacro("Error when updating tutte map");
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
int vtkSphericalConformalMapper::WolfeLineSearch(int map)
{
  int numPts = this->InitialPd->GetNumberOfPoints();

  vtkSmartPointer<vtkFloatArray> conjLaplacian =
    vtkSmartPointer<vtkFloatArray>::New();
  conjLaplacian->SetNumberOfComponents(3);
  conjLaplacian->Allocate(numPts, 10000);
  conjLaplacian->SetNumberOfTuples(numPts);
  this->ComputeDataArrayLaplacian(this->ConjugateDir, this->HarmonicMap[map],
                                  this->EdgeTable,
                                  this->EdgeWeights, this->EdgeNeighbors,
                                  conjLaplacian, map);

  double numerator[3];
  double denominator[3];
  this->VectorDotProduct(this->ConjugateDir, this->CurrDescent, numerator, numPts, 3);
  this->VectorDotProduct(this->ConjugateDir, conjLaplacian, denominator, numPts, 3);

  double maxstep = 100.0;
  for (int i=0; i<3; i++)
  {
    double alpha = numerator[i]/denominator[i];
    if (fabs(alpha) < maxstep)
    {
      maxstep = fabs(alpha);
    }
  }
  if (this->Verbose == 2)
  {
    fprintf(stdout, "New Step Size: %.5f\n", maxstep);
  }

  //Backtracking


  //vtkSmartPointer<vtkPolyData> tmpPoly =
  //  vtkSmartPointer<vtkPolyData>::New();
  //tmpPoly->DeepCopy(this->HarmonicMap[map]);
  //double E0 = 0.0;
  //this->ComputeEnergy(tmpPoly, E0, map);
  //for (int i=0; i<10; i++)
  //{
  //  for (int j=0; j<numPts; j++)
  //  {
  //    double newDescent[3], ptVal[3];
  //    this->ConjugateDir->GetTuple(j, newDescent);
  //    this->HarmonicMap[map]->GetPoint(j, ptVal);
  //    for (int k=0; k<3; k++)
  //    {
  //      ptVal[k] = ptVal[k] + maxstep * newDescent[k];
  //    }
  //    vtkMath::Normalize(ptVal);
  //    tmpPoly->GetPoints()->SetPoint(j, ptVal);
  //  }
  //  double EStep = 0.0;
  //  this->ComputeEnergy(tmpPoly, EStep, map);

  //  if( E0 - EStep < 0.0)
  //  {
  //    maxstep = 0.9*maxstep;
  //    EStep = E0;
  //    fprintf(stdout,"IT HAPPPEENEEED: %.8f\n", maxstep);
  //  }
  //  else
  //  {
  //    break;
  //  }
  //}
  this->TimeStep = maxstep;

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::SphericalTutteMapping()
{
  fprintf(stdout,"SphericalTutteMapping...\n");
  int numPts = this->InitialPd->GetNumberOfPoints();
  int numTris = this->InitialPd->GetNumberOfCells();

  double E0 = 0.0;
  this->ComputeEnergy(this->HarmonicMap[TUTTE], this->EdgeTable,
                      this->EdgeWeights, E0, TUTTE);
  double R0 = 1000.0;

  fprintf(stdout,"Starting Tutte Iterations...\n");
  if (this->Verbose)
  {
    fprintf(stdout,"Initial Tutte Energy: %.8f\n",E0);
  }
  if (this->FirstStep(TUTTE) != 1)
  {
    vtkErrorMacro("Error during initial step");
    return 0;
  }
  double EStep = 0.0;
  this->ComputeEnergy(this->HarmonicMap[TUTTE], this->EdgeTable,
                      this->EdgeWeights, EStep, TUTTE);
  double Ediff = E0-EStep;
  double RStep = 0.0;
  this->ComputeResidual(RStep);
  double Rdiff = R0-RStep;
  if (this->Verbose)
  {
    fprintf(stdout,"| Iter | 000000 | Tutte Energy | %16.8f | Res | %16.8f |\n", EStep, Ediff);
  }
  E0 = EStep;
  for (int iter=0; iter<this->MaxNumIterations; iter++)
  {
    vtkSmartPointer<vtkFloatArray> laplacian =
      vtkSmartPointer<vtkFloatArray>::New();
    laplacian->SetNumberOfComponents(3);
    laplacian->Allocate(numPts, 10000);
    laplacian->SetNumberOfTuples(numPts);
    if (this->ComputeMeshLaplacian(this->HarmonicMap[TUTTE], this->EdgeTable,
                                   this->EdgeWeights, this->EdgeNeighbors,
                                   laplacian, TUTTE) != 1)
    {
      vtkErrorMacro("Error when computing laplacian");
      return 0;
    }

    if (this->UpdateMap(laplacian, TUTTE, this->CGUpdateMethod) != 1)
    {
      vtkErrorMacro("Error when updating tutte map");
      return 0;
    }

    this->ComputeEnergy(this->HarmonicMap[TUTTE], this->EdgeTable,
                        this->EdgeWeights, EStep, TUTTE);
    Ediff = E0-EStep;
    this->ComputeResidual(RStep);
    Rdiff = R0-RStep;

    if (this->Verbose)
    {
      fprintf(stdout,"| Iter | %06d | Tutte Energy | %16.8f | Res | %16.8f |\n",iter+1, EStep, Ediff);
    }
    if (fabs(Ediff) < this->TutteEnergyCriterion)
    {
      fprintf(stdout,"Energy Criterion Met! %.5f\n",EStep);
      break;
    }
    else
    {
      E0 = EStep;
      R0 = RStep;
    }
    if (this->Verbose == 3)
    {
      if (iter%this->NumSaveIterations == 0)
      {
        if (this->IterOutputFilename != NULL)
        {
          std::stringstream iterstr;
          iterstr << this->SaveIter++;
          std::string iterName = this->IterOutputFilename;
          std::string filename =  iterName+"_"+iterstr.str()+".vtp";
          vtkSmartPointer<vtkXMLPolyDataWriter> writer =
            vtkSmartPointer<vtkXMLPolyDataWriter>::New();
          writer->SetInputData(this->HarmonicMap[TUTTE]);
          writer->SetFileName(filename.c_str());
          writer->Write();
        }
      }
    }
  }

  fprintf(stdout,"Done with SphericalTutteMapping...\n");
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::SphericalConformalMapper()
{
  fprintf(stdout,"SphericalConformalMapper...\n");
  int numPts = this->InitialPd->GetNumberOfPoints();
  int numTris = this->InitialPd->GetNumberOfCells();

  double E0 = 0.0;
  double R0 = 1000.0;
  this->ComputeEnergy(this->HarmonicMap[HARMONIC], this->EdgeTable,
                      this->EdgeWeights, E0, HARMONIC);

  fprintf(stdout,"Starting Harmonic Iterations...\n");
  if (this->Verbose)
  {
    fprintf(stdout,"Initial Harmonic Energy: %.8f\n",E0);
  }
  if (this->FirstStep(HARMONIC) != 1)
  {
    vtkErrorMacro("Error during initial step");
    return 0;
  }
  double EStep = 0.0;
  this->ComputeEnergy(this->HarmonicMap[HARMONIC], this->EdgeTable,
                      this->EdgeWeights, EStep, HARMONIC);
  double Ediff = E0-EStep;
  double RStep = 0.0;
  this->ComputeResidual(RStep);
  double Rdiff = R0-RStep;
  if (this->Verbose)
  {
    fprintf(stdout,"| Iter | 000000 | Harmonic Energy | %16.8f | Res | %16.8f |\n", EStep, Ediff);
  }
  E0 = EStep;
  for (int iter=0; iter<this->MaxNumIterations; iter++)
  {
    vtkSmartPointer<vtkFloatArray> laplacian =
      vtkSmartPointer<vtkFloatArray>::New();
    laplacian->SetNumberOfComponents(3);
    laplacian->Allocate(numPts, 10000);
    laplacian->SetNumberOfTuples(numPts);
    if (this->ComputeMeshLaplacian(this->HarmonicMap[HARMONIC], this->EdgeTable,
                                   this->EdgeWeights, this->EdgeNeighbors,
                                   laplacian, HARMONIC) != 1)
    {
      vtkErrorMacro("Error when computing laplacian");
      return 0;
    }

    if (this->UpdateMap(laplacian, HARMONIC, this->CGUpdateMethod) != 1)
    {
      vtkErrorMacro("Error when updating tutte map");
      return 0;
    }

    //Compute Mobius transformation
    if (this->NumBoundaries == 0)
    {
      if (this->ComputeMobiusTransformation() != 1)
      {
        vtkErrorMacro("Error when computing the mobius transformation");
        return 0;
      }
    }

    this->ComputeEnergy(this->HarmonicMap[HARMONIC], this->EdgeTable,
                        this->EdgeWeights, EStep, HARMONIC);
    Ediff = E0-EStep;
    this->ComputeResidual(RStep);
    Rdiff = R0-RStep;

    if (this->Verbose)
    {
      fprintf(stdout,"| Iter | %06d | Harmonic Energy | %16.8f | Res | %16.8f |\n",iter+1, EStep, Ediff);
    }
    if (fabs(Ediff) < this->HarmonicEnergyCriterion)
    {
      fprintf(stdout,"Energy Criterion Met! %.5f\n",EStep);
      break;
    }
    else
    {
      E0 = EStep;
      R0 = RStep;
    }
    if (this->Verbose == 3)
    {
      if (iter%this->NumSaveIterations == 0)
      {
        if (this->IterOutputFilename != NULL)
        {
          std::stringstream iterstr;
          iterstr << this->SaveIter++;
          std::string iterName = this->IterOutputFilename;
          std::string filename = iterName+"_"+iterstr.str()+".vtp";
          vtkSmartPointer<vtkXMLPolyDataWriter> writer =
            vtkSmartPointer<vtkXMLPolyDataWriter>::New();
          writer->SetInputData(this->HarmonicMap[HARMONIC]);
          writer->SetFileName(filename.c_str());
          writer->Write();
        }
      }
    }
  }

  fprintf(stdout,"Done with SphericalConformalMapper...\n");
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::ComputeMeshLaplacian(vtkPolyData *pd,
                                                      vtkEdgeTable *edgeTable,
                                                      vtkFloatArray *edgeWeights,
                                                      vtkIntArray *edgeNeighbors,
                                                      vtkFloatArray *laplacian, int map)
{
  int numPts = pd->GetNumberOfPoints();

  for (int i=0; i<numPts; i++)
  {
    double pointLaplacian[3];
    vtkSphericalConformalMapper::ComputePointLaplacian(i, pd,
                                                       edgeTable,
                                                       edgeWeights, edgeNeighbors,
                                                       pointLaplacian, map);
    laplacian->SetTuple(i, pointLaplacian);
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::ComputeDataArrayLaplacian(vtkFloatArray *data,
                                                           vtkPolyData *pd,
                                                           vtkEdgeTable *edgeTable,
                                                           vtkFloatArray *edgeWeights,
                                                           vtkIntArray *edgeNeighbors,
                                                           vtkFloatArray *laplacian, int map)
{
  int numPts = data->GetNumberOfTuples();

  for (int i=0; i<numPts; i++)
  {
    double pointLaplacian[3];
    vtkSphericalConformalMapper::ComputeDataLaplacian(i, data, pd, edgeTable, edgeWeights,
                                                        edgeNeighbors, pointLaplacian, map);
    laplacian->SetTuple(i, pointLaplacian);
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::ComputePointLaplacian(vtkIdType p0,
                                                         vtkPolyData *pd,
                                                         vtkEdgeTable *edgeTable,
                                                         vtkFloatArray *edgeWeights,
                                                         vtkIntArray *edgeNeighbors,
                                                         double laplacian[],
                                                         int map)
{
  vtkSmartPointer<vtkIdList> pointNeighbors = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> dummyList = vtkSmartPointer<vtkIdList>::New();
  vtkSphericalConformalMapper::GetPointNeighbors(p0, pd, pointNeighbors, dummyList);

  laplacian[0] = 0.0; laplacian[1] = 0.0; laplacian[2] = 0.0;
  for (int i=0; i<pointNeighbors->GetNumberOfIds(); i++)
  {
    vtkIdType p1 = pointNeighbors->GetId(i);
    vtkIdType edgeId = edgeTable->IsEdge(p0, p1);
    double weight = edgeWeights->GetValue(edgeId);
    if (map == TUTTE)
    {
      weight = 1.0;
    }
    int edgeNeighbor = edgeNeighbors->GetValue(edgeId);
    if (edgeNeighbor == -1)
    {
      continue;
    }
    double p0Metric[3], p1Metric[3], data0[3], data1[3];
    pd->GetPoint(p0, p0Metric);
    pd->GetPoint(p1, p1Metric);

    for (int j=0; j<3; j++)
    {
      laplacian[j] += weight * (p0Metric[j] - p1Metric[j]);
    }
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::ComputeDataLaplacian(vtkIdType p0,
                                                        vtkFloatArray *data,
                                                        vtkPolyData *pd,
                                                        vtkEdgeTable *edgeTable,
                                                        vtkFloatArray *edgeWeights,
                                                        vtkIntArray *edgeNeighbors,
                                                        double laplacian[],
                                                        int map)
{
  vtkSmartPointer<vtkIdList> pointNeighbors = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> dummyList = vtkSmartPointer<vtkIdList>::New();
  vtkSphericalConformalMapper::GetPointNeighbors(p0, pd, pointNeighbors, dummyList);

  laplacian[0] = 0.0; laplacian[1] = 0.0; laplacian[2] = 0.0;
  for (int i=0; i<pointNeighbors->GetNumberOfIds(); i++)
  {
    vtkIdType p1 = pointNeighbors->GetId(i);
    vtkIdType edgeId = edgeTable->IsEdge(p0, p1);
    double weight = edgeWeights->GetValue(edgeId);
    if (map == TUTTE)
    {
      weight = 1.0;
    }
    int edgeNeighbor = edgeNeighbors->GetValue(edgeId);
    if (edgeNeighbor == -1)
    {
      continue;
    }
    double p0Metric[3], p1Metric[3], data0[3], data1[3];
    data->GetTuple(p0, p0Metric);
    data->GetTuple(p1, p1Metric);

    for (int j=0; j<3; j++)
    {
      laplacian[j] += weight * (p0Metric[j] - p1Metric[j]);
    }
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::ComputeMassCenter(vtkPolyData *pd, double massCenter[])
{
  massCenter[0] = 0.0;
  massCenter[1] = 0.0;
  massCenter[2] = 0.0;
  vtkSmartPointer<vtkCenterOfMass> centerFinder =
    vtkSmartPointer<vtkCenterOfMass>::New();
  centerFinder->SetInputData(pd);
  centerFinder->Update();
  centerFinder->GetCenter(massCenter);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::ComputeMobiusTransformation()
{
  int numPts = this->InitialPd->GetNumberOfPoints();

  double massCenter[3];
  this->ComputeMassCenter(this->HarmonicMap[HARMONIC], massCenter);
  if (this->Verbose == 2)
  {
    fprintf(stdout,"Mass Center: %.16f, %.16f, %.16f\n",massCenter[0], massCenter[1], massCenter[2]);
  }

  for (int i =0; i<numPts; i++)
  {
    double h[3];
    this->HarmonicMap[HARMONIC]->GetPoint(i, h);
    for (int j=0; j<3; j++)
    {
      h[j] = h[j] - massCenter[j];
    }
    double hNorm = vtkMath::Norm(h);
    for (int j=0; j<3; j++)
    {
      h[j] = h[j]/hNorm;
    }
    if (this->IsBoundary->GetValue(i) == 0)
    {
      this->HarmonicMap[HARMONIC]->GetPoints()->SetPoint(i, h);
    }
    this->MassCenter[0] = massCenter[0];
    this->MassCenter[1] = massCenter[1];
    this->MassCenter[2] = massCenter[2];
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::ComputeEnergy(vtkPolyData *pd,
                                                 vtkEdgeTable *edgeTable,
                                                 vtkFloatArray *edgeWeights,
                                                 double &harmonicEnergy,
                                                 int map)
{
  harmonicEnergy = 0.0;
  double compEnergy[3];
  compEnergy[0] = 0.0; compEnergy[1] = 0.0; compEnergy[2] =0.0;
  int numEdges = edgeTable->GetNumberOfEdges();

  edgeTable->InitTraversal();
  for (int i=0; i<numEdges; i++)
  {
    vtkIdType p0, p1;
    vtkIdType edgeId = edgeTable->GetNextEdge(p0, p1);
    double weight = edgeWeights->GetValue(edgeId);
    if (map == TUTTE)
    {
      weight = 1.0;
    }

    double h0[3];
    double h1[3];
    pd->GetPoint(p0, h0);
    pd->GetPoint(p1, h1);

    //Calculate String Energy!
    double edgeEnergy[3];
    vtkSphericalConformalMapper::ComputeStringEnergy(h0, h1, weight, edgeEnergy);
    for (int j=0; j<3; j++)
    {
      compEnergy[j] += edgeEnergy[j];
    }
  }
  for (int i=0; i<3; i++)
  {
    harmonicEnergy += compEnergy[i];
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::ComputeStringEnergy(double e0[],
                                                      double e1[],
                                                      double weight,
                                                      double stringEnergy[])
{
  double edge[3];
  edge[0] = e0[0] - e1[0];
  edge[1] = e0[1] - e1[1];
  edge[2] = e0[2] - e1[2];

  for (int i=0; i<3; i++)
  {
    stringEnergy[i] = weight * pow(edge[i], 2);
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::ComputeResidual(double &residual)
{
  int numPts = this->InitialPd->GetNumberOfPoints();
  double newRes[3];
  this->VectorDotProduct(this->ConjugateDir, this->ConjugateDir, newRes,
                         numPts, 3);

  residual = 0.0;
  for (int i=0; i<3; i++)
  {
    newRes[i] = sqrt(newRes[i]);
    residual += newRes[i];
  }
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::UpdateMap(vtkFloatArray *laplacian,
                                             int map,
                                             int cg_update)
{
  int numPts = this->HarmonicMap[map]->GetNumberOfPoints();

  for (int i=0; i<numPts; i++)
  {
    double pointLaplacian[3], pointNormal[3];
    laplacian->GetTuple(i, pointLaplacian);
    this->HarmonicMap[map]->GetPoint(i, pointNormal);

    double pointScalar = vtkMath::Dot(pointLaplacian, pointNormal);
    double pointLaplacianNormal[3];
    double pointLaplacianTangential[3];
    for (int j=0; j<3; j++)
    {
      pointLaplacianNormal[j] = pointScalar * pointNormal[j];
      pointLaplacianTangential[j] = -1.0*(pointLaplacian[j] - pointLaplacianNormal[j]);
    }
    if (cg_update == CG_NONE)
    {
      this->PrevDescent->SetTuple(i, pointLaplacianTangential);
    }
    else
    {
      this->CurrDescent->SetTuple(i, pointLaplacianTangential);
    }
  }

  if (this->StepForward(map, cg_update) != 1)
  {
    vtkErrorMacro("Error when updating tutte map in CG step");
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
int vtkSphericalConformalMapper::VectorDotProduct(vtkFloatArray *v0, vtkFloatArray *v1, double product[], int numVals, int numComps)
{
  for (int i=0; i<numComps; i++)
  {
    product[i] = 0.0;
  }
  for (int i=0; i<numVals; i++)
  {
    for (int j=0; j<numComps; j++)
    {
      double val0, val1;
      val0 = v0->GetComponent(i, j);
      val1 = v1->GetComponent(i, j);
      product[j] += val0 * val1;
    }
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::VectorAdd(vtkFloatArray *v0, vtkFloatArray *v1, double scalar, vtkFloatArray *result, int numVals, int numComps)
{
  for (int i=0; i<numVals; i++)
  {
    for (int j=0; j<numComps; j++)
    {
      double val0, val1;
      val0 = v0->GetComponent(i, j);
      val1 = v1->GetComponent(i, j);
      double sum = val0 + scalar * val1;
      result->SetComponent(i, j, sum);
    }
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::StepForward(int map, int cg_update)
{

  if (cg_update == CG_NONE)
  {
    int numPts = this->InitialPd->GetNumberOfPoints();
    for (int i=0; i<numPts; i++)
    {
      double ptVal[3], descent[3];
      this->HarmonicMap[map]->GetPoint(i, ptVal);
      this->PrevDescent->GetTuple(i, descent);
      this->ConjugateDir->SetTuple(i, descent);
      if (this->IsBoundary->GetValue(i) == 0)
      {
        for (int j=0; j<3; j++)
        {
          ptVal[j] = ptVal[j] + this->InitialTimeStep * descent[j];
        }
        vtkMath::Normalize(ptVal);
      }
      this->HarmonicMap[map]->GetPoints()->SetPoint(i, ptVal);
    }
    return 1;
  }
  else if (cg_update == CG_FLETCHER_REEVES)
  {
    this->FRUpdateMap(map);
  }
  else if (cg_update == CG_POLAK_RIBIERE)
  {
    this->PRUpdateMap(map);
  }
  else if (cg_update == CG_HESTENESS_STIEFEL)
  {
    this->HSUpdateMap(map);
  }
  else if (cg_update == CG_DAI_YUAN)
  {
    this->DYUpdateMap(map);
  }
  else
  {
    fprintf(stderr,"No correct option give\n");
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
int vtkSphericalConformalMapper::FRUpdateMap(int map)
{
  int numPts = this->PrevDescent->GetNumberOfTuples();
  double numerator[3], denominator[3];
  this->VectorDotProduct(this->CurrDescent, this->CurrDescent,
                         numerator, numPts, 3);
  this->VectorDotProduct(this->PrevDescent, this->PrevDescent,
                         denominator, numPts, 3);

  double beta[3];
  for (int i=0; i<3; i++)
  {
    beta[i] = numerator[i]/denominator[i];
  }

  this->CGUpdateMap(map, beta);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::PRUpdateMap(int map)
{
  int numPts = this->PrevDescent->GetNumberOfTuples();
  double numerator[3], denominator[3];
  vtkSmartPointer<vtkFloatArray> difference =
    vtkSmartPointer<vtkFloatArray>::New();
  difference->SetNumberOfComponents(3);
  difference->Allocate(numPts, 10000);
  difference->SetNumberOfTuples(numPts);
  this->VectorAdd(this->CurrDescent, this->PrevDescent, -1.0,
                  difference, numPts, 3);
  this->VectorDotProduct(this->CurrDescent, difference,
                         numerator, numPts, 3);
  this->VectorDotProduct(this->PrevDescent, this->PrevDescent,
                         denominator, numPts, 3);

  double beta[3];
  for (int i=0; i<3; i++)
  {
    beta[i] = std::max(0.0, numerator[i]/denominator[i]);
  }

  this->CGUpdateMap(map, beta);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::HSUpdateMap(int map)
{
  int numPts = this->PrevDescent->GetNumberOfTuples();
  double numerator[3], denominator[3];
  vtkSmartPointer<vtkFloatArray> difference =
    vtkSmartPointer<vtkFloatArray>::New();
  difference->SetNumberOfComponents(3);
  difference->Allocate(numPts, 10000);
  difference->SetNumberOfTuples(numPts);
  this->VectorAdd(this->CurrDescent, this->PrevDescent, -1.0,
                  difference, numPts, 3);
  this->VectorDotProduct(this->CurrDescent, difference,
                         numerator, numPts, 3);
  this->VectorDotProduct(this->ConjugateDir, difference,
                         denominator, numPts, 3);

  double beta[3];
  for (int i=0; i<3; i++)
  {
    beta[i] = -1.0 *numerator[i]/denominator[i];
  }

  this->CGUpdateMap(map, beta);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::DYUpdateMap(int map)
{
  int numPts = this->PrevDescent->GetNumberOfTuples();
  double numerator[3], denominator[3];
  vtkSmartPointer<vtkFloatArray> difference =
    vtkSmartPointer<vtkFloatArray>::New();
  difference->SetNumberOfComponents(3);
  difference->Allocate(numPts, 10000);
  difference->SetNumberOfTuples(numPts);
  this->VectorAdd(this->CurrDescent, this->PrevDescent, -1.0,
                  difference, numPts, 3);
  this->VectorDotProduct(this->CurrDescent, this->CurrDescent,
                         numerator, numPts, 3);
  this->VectorDotProduct(this->ConjugateDir, difference,
                         denominator, numPts, 3);

  double beta[3];
  for (int i=0; i<3; i++)
  {
    beta[i] = -1.0 * numerator[i]/denominator[i];
  }

  this->CGUpdateMap(map, beta);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::CGUpdateMap(int map, double beta[])
{
  int numPts = this->InitialPd->GetNumberOfPoints();

  double descCond[3];
  this->VectorDotProduct(this->CurrDescent, this->ConjugateDir, descCond,
                         numPts, 3);

  //If descent condition isn't satisfied, restart from new steepest descent dir
  for (int i=0; i<3; i++)
  {
    if (descCond[i] <= 0.0)
    {
      beta[0] = 0.0; beta[1] = 0.0; beta[2] = 0.0;
    }
  }

  for (int i=0; i<numPts; i++)
  {
    double conjdir[3], descent[3], newDescent[3];
    this->CurrDescent->GetTuple(i, descent);
    this->ConjugateDir->GetTuple(i, conjdir);
    for (int j=0; j<3; j++)
    {
      newDescent[j] = descent[j] + beta[j]*conjdir[j];
    }
    this->PrevDescent->SetTuple(i, descent);
    this->ConjugateDir->SetTuple(i, newDescent);
  }

  this->WolfeLineSearch(map);

  for (int i=0; i<numPts; i++)
  {
    double newDescent[3], ptVal[3];
    this->ConjugateDir->GetTuple(i, newDescent);
    this->HarmonicMap[map]->GetPoint(i, ptVal);
    if (this->IsBoundary->GetValue(i) == 0)
    {
      for (int j=0; j<3; j++)
      {
        {
          ptVal[j] = ptVal[j] + this->TimeStep * newDescent[j];
        }
      }
      vtkMath::Normalize(ptVal);
    }
    this->HarmonicMap[map]->GetPoints()->SetPoint(i, ptVal);
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::GetPointNeighbors(vtkIdType p0,
                                                     vtkPolyData *pd,
						     vtkIdList *pointNeighbors,
						     vtkIdList *cellIdList)
{
  //Assuming that pointNeighbors is set with no neighbors already
  pd->GetPointCells(p0, cellIdList);

  for (int i=0; i<cellIdList->GetNumberOfIds(); i++)
  {
    vtkIdType cellId = cellIdList->GetId(i);
    vtkIdType npts, *pts;
    pd->GetCellPoints(cellId, npts, pts);

    for (int j=0; j<npts; j++)
    {
      vtkIdType neighborPoint = pts[j];
      if (neighborPoint != p0)
      {
        pointNeighbors->InsertUniqueId(neighborPoint);
      }
    }
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
int vtkSphericalConformalMapper::SeparateLoops(vtkPolyData *pd,
                                               vtkPolyData **loops,
                                               int numBoundaries,
                                               const double xvec[3],
                                               const double zvec[3],
                                               const int boundaryStart[2])
{
  vtkIdType nextCell;
  vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
  int numInterPts = pd->GetNumberOfPoints();
  int numInterLines = pd->GetNumberOfLines();
  pd->BuildLinks();

  int count = 0;
  for (int i=0;i<numBoundaries;i++)
  {
    vtkIdType startPt = boundaryStart[i];
    vtkPolyData *newloop = loops[count];
    newloop->Allocate(pd->GetNumberOfCells(), 1000);
    pd->GetPointCells(startPt,cellIds);

    nextCell = cellIds->GetId(0);
    vtkIdType npts, *pts;
    int testPt = -1;
    pd->GetCellPoints(nextCell, npts, pts);
    if (pts[0] == startPt)
      testPt = pts[1];
    else
      testPt = pts[0];

    double pt0[3], pt1[3], vec0[3], vec1[3];
    pd->GetPoint(startPt, pt0);
    pd->GetPoint(testPt, pt1);
    vtkMath::Subtract(pt1, pt0, vec0);
    vtkMath::Normalize(vec0);
    vtkMath::Cross(zvec, xvec, vec1);
    vtkMath::Normalize(vec1);
    if (vtkMath::Dot(vec0, vec1) < 0)
    {
      nextCell = cellIds->GetId(1);
    }
    //if (testPt != boundaryStart[i+2])
    //{
    //  nextCell = cellIds->GetId(1);
    //}

    //Run through intersection lines to get loops!
    vtkSphericalConformalMapper::RunLoopFind(pd, startPt, nextCell, newloop);
    loops[count++] = newloop;
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::RunLoopFind(vtkPolyData *pd,
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
int vtkSphericalConformalMapper::CalculateCircleLength(vtkPolyData *lines,
                                                       double &length)
{
  int numLines = lines->GetNumberOfLines();

  length = 0.0;
  for (int i=0; i<numLines; i++)
  {
    vtkIdType npts, *pts;
    lines->GetCellPoints(i, npts, pts);

    double pt0[3], pt1[3];
    lines->GetPoint(pts[0], pt0);
    lines->GetPoint(pts[1], pt1);

    double dist = std::sqrt(std::pow(pt0[0]-pt1[0], 2.0) +
                            std::pow(pt0[1]-pt1[1], 2.0) +
                            std::pow(pt0[2]-pt1[2], 2.0));
    length += dist;
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::CalculateSquareEdgeLengths(vtkPolyData *lines,
                                                           vtkIntArray *markerPts,
                                                           double lengths[])
{
  vtkIntArray *pointIds = vtkIntArray::SafeDownCast(lines->GetPointData()->GetArray("PointIds"));
  int numLines = lines->GetNumberOfLines();

  for (int i=0; i<markerPts->GetNumberOfTuples(); i++)
  {
    lengths[i] = 0.0;
  }
  int currCell = 0;
  for (int i=0; i<markerPts->GetNumberOfTuples(); i++)
  {
    int lastPt  = markerPts->GetValue(i);
    vtkIdType npts, *pts;
    int checkPt = -1;
    while (checkPt != lastPt)
    {
      //fprintf(stdout,"What is: %d\n", checkPt);
      lines->GetCellPoints(currCell, npts, pts);
      double pt0[3], pt1[3];
      lines->GetPoint(pts[0], pt0);
      lines->GetPoint(pts[1], pt1);
      checkPt = pointIds->GetValue(pts[1]);

      double dist = std::sqrt(std::pow(pt0[0]-pt1[0], 2.0) +
                              std::pow(pt0[1]-pt1[1], 2.0) +
                              std::pow(pt0[2]-pt1[2], 2.0));
      lengths[i] += dist;

      currCell++;
    }
  }
  //fprintf(stdout,"Curr!: %d\n", currCell);
  //fprintf(stdout,"NumLines!: %d\n", numLines);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::SetLoopOnUnitCircle(vtkPolyData *lines,
                                                    double length,
                                                    double radius)
{
  vtkIntArray *pointIds = vtkIntArray::SafeDownCast(lines->GetPointData()->GetArray("PointIds"));
  int numLines = lines->GetNumberOfLines();

  double currLength = 0.0;
  for (int i=0; i<numLines; i++)
  {
    vtkIdType npts, *pts;
    lines->GetCellPoints(i, npts, pts);

    double pt0[3], pt1[3];
    lines->GetPoint(pts[0], pt0);
    lines->GetPoint(pts[1], pt1);

    double dist = std::sqrt(std::pow(pt0[0]-pt1[0], 2.0) +
                            std::pow(pt0[1]-pt1[1], 2.0) +
                            std::pow(pt0[2]-pt1[2], 2.0));
    currLength += dist;

    double angle = currLength/length * 2.0 * M_PI;
    double x_val = (2.0 * radius * std::cos(angle))/(1.0 + std::pow(radius, 2.0));
    double y_val = (2.0 * radius * std::sin(angle))/(1.0 + std::pow(radius, 2.0));
    double z_val = (-1.0 + std::pow(radius, 2.0))/(1.0 + std::pow(radius, 2.0));

    int id = pointIds->GetValue(pts[1]);

    this->HarmonicMap[0]->GetPoints()->SetPoint(id, x_val, y_val, z_val);
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::SetCircleBoundary(vtkPolyData *lines,
                                                   vtkIntArray *markerPts,
                                                   vtkIntArray *markerDirs,
                                                   double cubeStart[],
                                                   double lengths[],
                                                   double radius)
{
  vtkIntArray *pointIds = vtkIntArray::SafeDownCast(lines->GetPointData()->GetArray("PointIds"));
  int numLines = lines->GetNumberOfLines();

  double currCoords[3];
  for (int i=0; i<3; i++)
  {
    currCoords[i] = cubeStart[i];
  }

  double unitLength = M_PI / 2.0;
  int currCell = 0;
  int checkPt = -1;
  for (int i=0; i<markerPts->GetNumberOfTuples(); i++)
  {
    double currLength = 0.0;
    int lastPt  = markerPts->GetValue(i);
    int dir     = markerDirs->GetValue(i);
    vtkIdType npts, *pts;
    while (checkPt != lastPt)
    {
      lines->GetCellPoints(currCell, npts, pts);
      double pt0[3], pt1[3];
      lines->GetPoint(pts[0], pt0);
      lines->GetPoint(pts[1], pt1);

      checkPt = pointIds->GetValue(pts[1]);

      double dist = std::sqrt(std::pow(pt0[0]-pt1[0], 2.0) +
                              std::pow(pt0[1]-pt1[1], 2.0) +
                              std::pow(pt0[2]-pt1[2], 2.0));
      currLength += dist;

      double boundaryVal[3];
      //if (dir == 7)
      //{
      //  double transPt[3];
      //  double angle = currLength/lengths[i] * unitLength + 5.0 * M_PI/4.0;
      //  transPt[0] = (2.0 * radius * std::cos(angle))/(1.0 + std::pow(radius, 2.0));
      //  transPt[1] = (2.0 * radius * std::sin(angle))/(1.0 + std::pow(radius, 2.0));
      //  transPt[2] = (-1.0 + std::pow(radius, 2.0))/(1.0 + std::pow(radius, 2.0));
      //  vtkSphericalConformalMapper::RotateByAngle(transPt, 180.0, boundaryVal);
      //}
      //else if (dir == 2)
      //{
      //  double transPt[3];
      //  double angle = currLength/lengths[i] * unitLength + 5.0 * M_PI/4.0;
      //  transPt[0] = (2.0 * radius * std::cos(angle))/(1.0 + std::pow(radius, 2.0));
      //  transPt[1] = (2.0 * radius * std::sin(angle))/(1.0 + std::pow(radius, 2.0));
      //  transPt[2] = (-1.0 + std::pow(radius, 2.0))/(1.0 + std::pow(radius, 2.0));
      //  vtkSphericalConformalMapper::RotateByAngle(transPt, 90.0, boundaryVal);
      //}
      //else if (dir == 3 && i != 2)
      //{
      //  double transPt[3];
      //  double angle = currLength/lengths[i] * unitLength + M_PI/4.0;
      //  transPt[0] = (2.0 * radius * std::cos(angle))/(1.0 + std::pow(radius, 2.0));
      //  transPt[1] = (2.0 * radius * std::sin(angle))/(1.0 + std::pow(radius, 2.0));
      //  transPt[2] = (-1.0 + std::pow(radius, 2.0))/(1.0 + std::pow(radius, 2.0));
      //  vtkSphericalConformalMapper::RotateByAngle(transPt, 180.0, boundaryVal);
      //}
      //else if (dir == 4)
      //{
      //  double transPt[3];
      //  double angle = currLength/lengths[i] * unitLength +  7.0 * M_PI/4.0;
      //  transPt[0] = (2.0 * radius * std::cos(angle))/(1.0 + std::pow(radius, 2.0));
      //  transPt[1] = (2.0 * radius * std::sin(angle))/(1.0 + std::pow(radius, 2.0));
      //  transPt[2] = (-1.0 + std::pow(radius, 2.0))/(1.0 + std::pow(radius, 2.0));
      //  vtkSphericalConformalMapper::RotateByAngle(transPt, 180.0, boundaryVal);
      //}
      //else if (dir == 5)
      //{
      //  double transPt[3];
      //  double angle = currLength/lengths[i] * unitLength +  M_PI/4.0;
      //  transPt[0] = (2.0 * radius * std::cos(angle))/(1.0 + std::pow(radius, 2.0));
      //  transPt[1] = (2.0 * radius * std::sin(angle))/(1.0 + std::pow(radius, 2.0));
      //  transPt[2] = (-1.0 + std::pow(radius, 2.0))/(1.0 + std::pow(radius, 2.0));
      //  vtkSphericalConformalMapper::RotateByAngle(transPt, 90.0, boundaryVal);
      //}
      //else
      //{
        double angle = currLength/lengths[i] * unitLength + unitLength*i;// + unitLength/2.0;
        boundaryVal[0] = (2.0 * radius * std::cos(angle))/(1.0 + std::pow(radius, 2.0));
        boundaryVal[1] = (2.0 * radius * std::sin(angle))/(1.0 + std::pow(radius, 2.0));
        boundaryVal[2] = (-1.0 + std::pow(radius, 2.0))/(1.0 + std::pow(radius, 2.0));
      //}

      int id = pointIds->GetValue(pts[1]);

      this->HarmonicMap[0]->GetPoints()->SetPoint(id, boundaryVal);

      currCell++;
    }
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::RotateByAngle(const double pt[3], const double angle,
                                               double returnPt[3])
{

  vtkSmartPointer<vtkTransform> transformer =
    vtkSmartPointer<vtkTransform>::New();
  transformer->RotateY(angle);
  transformer->TransformPoint(pt, returnPt);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::GetPointWeights(double f[], double pt0[],
                                                   double pt1[], double pt2[],
					           double &a0, double &a1, double &a2)
{
  double f0[3], f1[3], f2[3], v0[3], v1[3];
  for (int i=0; i<3; i++)
  {
    v0[i] = pt0[i] - pt1[i];
    v1[i] = pt0[i] - pt2[i];
    f0[i] = pt0[i] - f[i];
    f1[i] = pt1[i] - f[i];
    f2[i] = pt2[i] - f[i];
  }

  double vArea[3], vA0[3], vA1[3], vA2[3];
  vtkMath::Cross(v0, v1, vArea);
  double area = vtkMath::Norm(vArea);
  vtkMath::Cross(f1, f2, vA0);
  a0 = vtkMath::Norm(vA0)/area;
  vtkMath::Cross(f2, f0, vA1);
  a1 = vtkMath::Norm(vA1)/area;
  vtkMath::Cross(f0, f1, vA2);
  a2 = vtkMath::Norm(vA2)/area;

  return 1;
}

// ----------------------
// PDCheckArrayName
// ----------------------
/**
 * @brief Function to check is array with name exists in cell or point data
 * @param object this is the object to check if the array exists
 * @param datatype this is point or cell. point =0,cell=1
 * @param arrayname this is the name of the array to check
 * @reutrn this returns 1 if the array exists and zero if it doesn't
 * or the function does not return properly.
 */

int vtkSphericalConformalMapper::PDCheckArrayName(vtkPolyData *object,int datatype,std::string arrayname )
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

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::SetCubeBoundary(vtkPolyData *lines,
                                                 vtkIntArray *markerPts,
                                                 vtkIntArray *markerDirs,
                                                 double cubeStart[],
                                                 double lengths[])
{
  vtkIntArray *pointIds = vtkIntArray::SafeDownCast(lines->GetPointData()->GetArray("PointIds"));
  int numLines = lines->GetNumberOfLines();

  double currCoords[3];
  for (int i=0; i<3; i++)
  {
    currCoords[i] = cubeStart[i];
  }
  double unitLength = 2.0;
  int currCell = 0;
  int checkPt = -1;
  for (int i=0; i<markerPts->GetNumberOfTuples(); i++)
  {
    double currLength = 0.0;
    int lastPt  = markerPts->GetValue(i);
    int dir     = markerDirs->GetValue(i);
    vtkIdType npts, *pts;
    while (checkPt != lastPt)
    {
      lines->GetCellPoints(currCell, npts, pts);
      double pt0[3], pt1[3];
      lines->GetPoint(pts[0], pt0);
      lines->GetPoint(pts[1], pt1);
      checkPt = pointIds->GetValue(pts[1]);

      double dist = std::sqrt(std::pow(pt0[0]-pt1[0], 2.0) +
                              std::pow(pt0[1]-pt1[1], 2.0) +
                              std::pow(pt0[2]-pt1[2], 2.0));
      currLength += dist;


      if (dir == 0)
      {
        currCoords[0] -= dist/lengths[i] * unitLength;
      }
      else if (dir == 1)
      {
        currCoords[1] -= dist/lengths[i] * unitLength;
      }
      else if (dir == 2)
      {
        currCoords[2] -= dist/lengths[i] * unitLength;
      }
      else if (dir == 3)
      {
        currCoords[0] += dist/lengths[i] * unitLength;
      }
      else if (dir == 4)
      {
        currCoords[1] += dist/lengths[i] * unitLength;
      }
      else if (dir == 5)
      {
        currCoords[2] += dist/lengths[i] * unitLength;
      }

      double sphereCoords[3];
      this->CubeBoundaryToSphere(currCoords, sphereCoords);

      int id = pointIds->GetValue(pts[1]);
      this->HarmonicMap[0]->GetPoints()->SetPoint(id, sphereCoords);

      currCell++;
    }
  }
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::CubeBoundaryToSphere(double inCoords[], double outCoords[])
{
  double x2 = std::pow(inCoords[0], 2.0);
  double y2 = std::pow(inCoords[1], 2.0);
  double z2 = std::pow(inCoords[2], 2.0);

  outCoords[0] = inCoords[0] * std::sqrt(1 - (y2/2.0) - (z2/2.0) + (y2*z2/3.0));
  outCoords[1] = inCoords[1] * std::sqrt(1 - (x2/2.0) - (z2/2.0) + (x2*z2/3.0));
  outCoords[2] = inCoords[2] * std::sqrt(1 - (y2/2.0) - (x2/2.0) + (y2*x2/3.0));

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSphericalConformalMapper::DetermineBoundaryPlan(int &numLoops, int bBool[])
{
  int t = this->BoundaryType;

  numLoops = 0;
  int opps = 0;
  int bounds = 0;
  if (t == 0)
  {
    return 1;
  }

  //Determine number of boundaries and opposite boundaries
  if ((t & NORTH) == 0)
  {
    bounds += 1;
    bBool[0] = 1;
    if ((t & SOUTH) == 0)
    {
      opps += 1;
    }
  }
  if ((t & SOUTH) == 0)
  {
    bBool[1] = 1;
    bounds += 1;
  }
  if ((t & FRONT) == 0)
  {
    bBool[2] = 1;
    bounds += 1;
    if ((t & BACK) == 0)
    {
      opps += 1;
    }
  }
  if ((t & BACK) == 0)
  {
    bBool[3] = 1;
    bounds += 1;
  }
  if ((t & LEFT) == 0)
  {
    bBool[4] = 1;
    bounds += 1;
    if ((t & RIGHT) == 0)
    {
      opps += 1;
    }
  }
  if ((t & RIGHT) == 0)
  {
    bBool[5] = 1;
    bounds += 1;
  }

  //Set number of boundary loops and number of boundary points based
  //on information provided
  if (opps == 0)
  {
    numLoops = 1;
  }
  else if (opps == 1)
  {
    if (bounds == 2)
    {
      numLoops = 2;
    }
    else
    {
      numLoops = 1;
    }
  }
  else if (opps == 2)
  {
    if (bounds == 4)
    {
      numLoops = 2;
    }
    else
    {
      numLoops = 1;
    }
  }
  else
  {
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
int vtkSphericalConformalMapper::GetCubeStartPoint(int id, double startCoords[])
{
  if (id ==0)
  {
    startCoords[0] = 1.0; startCoords[1] = 1.0; startCoords[2] = 1.0;
  }
  else if (id == 1)
  {
    startCoords[0] = -1.0; startCoords[1] = 1.0; startCoords[2] = 1.0;
  }
  else if (id == 2)
  {
    startCoords[0] = -1.0; startCoords[1] = -1.0; startCoords[2] = 1.0;
  }
  else if (id == 3)
  {
    startCoords[0] = 1.0; startCoords[1] = -1.0; startCoords[2] = 1.0;
  }
  else if (id == 4)
  {
    startCoords[0] = 1.0; startCoords[1] = 1.0; startCoords[2] = -1.0;
  }
  else if (id == 5)
  {
    startCoords[0] = -1.0; startCoords[1] = 1.0; startCoords[2] = -1.0;
  }
  else if (id == 6)
  {
    startCoords[0] = -1.0; startCoords[1] = -1.0; startCoords[2] = -1.0;
  }
  else if (id == 7)
  {
    startCoords[0] = 1.0; startCoords[1] = -1.0; startCoords[2] = -1.0;
  }

  return 1;
}
