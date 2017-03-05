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

/** @file vtkS2LandmarkMatcher.cxx
 *  @brief This implements the vtkS2LandmarkMatcher filter as a class
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkS2LandmarkMatcher.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkCenterOfMass.h"
#include "vtkCleanPolyData.h"
#include "vtkFloatArray.h"
#include "vtkGenericCell.h"
#include "vtkGradientFilter.h"
#include "vtkLoopSubdivisionFilter.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointLocator.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmartPointer.h"
#include "vtkSphericalConformalMapper.h"
#include "vtkTriangle.h"
#include "vtkWarpVector.h"
#include "vtkXMLPolyDataWriter.h"

#include "math.h"
#include "sparse_matrix.h"

#include <iostream>
#include <cmath>

//---------------------------------------------------------------------------
//vtkCxxRevisionMacro(vtkS2LandmarkMatcher, "$Revision: 0.0 $");
vtkStandardNewMacro(vtkS2LandmarkMatcher);


//---------------------------------------------------------------------------
vtkS2LandmarkMatcher::vtkS2LandmarkMatcher()
{
  this->SetNumberOfInputPorts(2);

  this->Verbose                 = 1;
  this->LandmarkTimeStep        = 0.001;
  this->LandmarkEnergyCriterion = 0.000001;
  this->LandmarkWeighting       = 1.0;
  this->MaxNumIterations           = 1000;
  this->NumLandmarks            = 0;
  this->CGUpdateMethod          = vtkSphericalConformalMapper::CG_NONE;

  this->SourcePd = vtkPolyData::New();
  this->TargetPd = vtkPolyData::New();
  this->MappedPd = vtkPolyData::New();

  this->SourceEdgeTable     = vtkEdgeTable::New();
  this->SourceEdgeWeights   = vtkFloatArray::New();
  this->PrevDescent         = vtkFloatArray::New();
  this->CurrDescent         = vtkFloatArray::New();
  this->ConjugateDir        = vtkFloatArray::New();
  this->SourceEdgeNeighbors = vtkIntArray::New();
  this->SourceIsBoundary    = vtkIntArray::New();

  this->SourceLandmarkPtIds = NULL;
  this->IsSourceLandmark    = vtkIntArray::New();
  this->TargetLandmarkPtIds = NULL;
  this->IsTargetLandmark    = vtkIntArray::New();
}

//---------------------------------------------------------------------------
vtkS2LandmarkMatcher::~vtkS2LandmarkMatcher()
{
  if (this->SourcePd != NULL)
  {
    this->SourcePd->Delete();
  }
  if (this->TargetPd != NULL)
  {
    this->TargetPd->Delete();
  }
  if (this->MappedPd != NULL)
  {
    this->MappedPd->Delete();
  }
  if (this->SourceEdgeTable != NULL)
  {
    this->SourceEdgeTable->Delete();
  }
  if (this->SourceEdgeWeights != NULL)
  {
    this->SourceEdgeWeights->Delete();
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
  if (this->SourceEdgeNeighbors != NULL)
  {
    this->SourceEdgeNeighbors->Delete();
  }
  if (this->SourceIsBoundary != NULL)
  {
    this->SourceIsBoundary->Delete();
  }
}

//---------------------------------------------------------------------------
void vtkS2LandmarkMatcher::PrintSelf(ostream& os, vtkIndent indent)
{
}

// Generate Separated Surfaces with Region ID Numbers
//---------------------------------------------------------------------------
int vtkS2LandmarkMatcher::RequestData(
                                 vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector)
{
  // get the input and output
  vtkPolyData *input1 = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData *input2 = vtkPolyData::GetData(inputVector[1]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  //Copy the input to operate on
  this->SourcePd->DeepCopy(input1);
  this->TargetPd->DeepCopy(input2);

  vtkIdType numSourcePolys = this->SourcePd->GetNumberOfPolys();
  //Check the input to make sure it is there
  if (numSourcePolys < 1)
  {
    vtkDebugMacro("No input!");
    return 0;
  }
  //Check the input to make sure it is manifold and a triangulated surface
  if (vtkSphericalConformalMapper::CleanAndCheckSurface(this->SourcePd) != 1)
  {
    vtkErrorMacro("Error when checking source surface");
    return 0;
  }
  vtkIdType numTargetPolys = this->TargetPd->GetNumberOfPolys();
  //Check the input to make sure it is there
  if (numTargetPolys < 1)
  {
    vtkDebugMacro("No input!");
    return 0;
  }
  //Check the input to make sure it is manifold and a triangulated surface
  if (vtkSphericalConformalMapper::CleanAndCheckSurface(this->TargetPd) != 1)
  {
    vtkErrorMacro("Error when checking target surface");
    return 0;
  }

  if (vtkSphericalConformalMapper::CreateEdgeTable(this->SourcePd,
                                                     this->SourceEdgeTable,
                                                     this->SourceEdgeWeights,
                                                     this->SourceEdgeNeighbors,
                                                     this->SourceIsBoundary) != 1)
  {
    vtkErrorMacro("Error when creating edge table for source polydata");
    return 0;
  }

  if (this->RunLandmarkMatching() != 1)
  {
    vtkErrorMacro("Error during mapping");
    return 0;
  }

  output->DeepCopy(this->MappedPd);
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkS2LandmarkMatcher::RunLandmarkMatching()
{

  if (this->SetLandmarks() != 1)
  {
    vtkErrorMacro("Error when setting landmarks");
    return 0;
  }

  if (this->InitiateCGArrays() != 1)
  {
    vtkErrorMacro("Error when initiating arrays");
    return 0;
  }

  this->MappedPd->DeepCopy(this->SourcePd);
  //Run the Lanmdark Energy Step
  if (this->LandmarkMapping() != 1)
  {
    vtkErrorMacro("Error when computing the tutte map");
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
int vtkS2LandmarkMatcher::SetLandmarks()
{
  if (this->SourceLandmarkPtIds != NULL && this->TargetLandmarkPtIds != NULL)
  {
    int numSourceLandmarks = this->SourceLandmarkPtIds->GetNumberOfTuples();
    int numTargetLandmarks = this->TargetLandmarkPtIds->GetNumberOfTuples();

    if (numSourceLandmarks != numTargetLandmarks)
    {
      vtkErrorMacro("Number of source and target landmarks must be equal");
      return 0;
    }

    this->NumLandmarks = numSourceLandmarks;
    int numPoints = this->SourcePd->GetNumberOfPoints();
    this->IsSourceLandmark->Allocate(numPoints, 10000);
    this->IsSourceLandmark->SetNumberOfValues(numPoints);
    this->IsTargetLandmark->Allocate(numPoints, 10000);
    this->IsTargetLandmark->SetNumberOfValues(numPoints);

    for (int i=0; i<numPoints; i++)
    {
      this->IsSourceLandmark->SetValue(i, -1);
      this->IsTargetLandmark->SetValue(i, -1);
    }
    for (int i=0; i<this->NumLandmarks; i++)
    {
      int srcId = this->SourceLandmarkPtIds->GetValue(i);
      this->IsSourceLandmark->SetValue(srcId, i);
      int targId = this->TargetLandmarkPtIds->GetValue(i);
      this->IsTargetLandmark->SetValue(targId, i);
    }
  }
  else
  {
    this->NumLandmarks = 0;
    this->IsSourceLandmark->SetNumberOfValues(0);
    this->IsTargetLandmark->SetNumberOfValues(0);
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkS2LandmarkMatcher::InitiateCGArrays()
{
  int numPts = this->SourcePd->GetNumberOfPoints();

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
int vtkS2LandmarkMatcher::LandmarkMapping()
{
  fprintf(stdout,"LandmarkMapping...\n");
  int numPts = this->MappedPd->GetNumberOfPoints();
  int numTris = this->MappedPd->GetNumberOfCells();

  double E0 = 0.0;
  this->ComputeLandmarkEnergy(this->MappedPd, this->SourceEdgeTable,
                              this->SourceEdgeWeights, E0,
                              this->LandmarkWeighting,
                              this->TargetPd,
                              this->IsSourceLandmark,
                              this->SourceLandmarkPtIds);

  fprintf(stdout,"Starting Landmark Iterations...\n");
  if (this->Verbose)
  {
    fprintf(stdout,"Initial Landmark Energy: %.8f\n",E0);
  }
  if (this->FirstStep() != 1)
  {
    vtkErrorMacro("Error during initial step");
    return 0;
  }
  double EStep = 0.0;
  this->ComputeLandmarkEnergy(this->MappedPd, this->SourceEdgeTable,
                              this->SourceEdgeWeights, EStep,
                              this->LandmarkWeighting,
                              this->TargetPd,
                              this->IsSourceLandmark,
                              this->SourceLandmarkPtIds);
  double Ediff = E0-EStep;
  if (this->Verbose)
  {
    fprintf(stdout,"| Iter | 000000 | Landmark Energy | %16.8f | Res | %16.8f |\n", EStep, Ediff);
  }
  E0 = EStep;
  for (int iter=0; iter<this->MaxNumIterations; iter++)
  {
    vtkSmartPointer<vtkFloatArray> laplacian =
      vtkSmartPointer<vtkFloatArray>::New();
    laplacian->SetNumberOfComponents(3);
    laplacian->Allocate(numPts, 10000);
    laplacian->SetNumberOfTuples(numPts);
    if (this->ComputeMeshLandmarkLaplacian(this->MappedPd, this->SourceEdgeTable,
                                           this->SourceEdgeWeights, this->SourceEdgeNeighbors,
                                           laplacian, this->LandmarkWeighting,
                                           this->TargetPd, this->IsSourceLandmark,
                                           this->SourceLandmarkPtIds) != 1)
    {
      vtkErrorMacro("Error when computing laplacian");
      return 0;
    }

    if (this->UpdateMap(laplacian, this->CGUpdateMethod) != 1)
    {
      vtkErrorMacro("Error when updating landmark map");
      return 0;
    }

    this->ComputeLandmarkEnergy(this->MappedPd, this->SourceEdgeTable,
                                this->SourceEdgeWeights, EStep,
                                this->LandmarkWeighting,
                                this->TargetPd,
                                this->IsSourceLandmark,
                                this->SourceLandmarkPtIds);
    Ediff = E0-EStep;

    if (this->Verbose)
    {
      fprintf(stdout,"| Iter | %06d | Landmark Energy | %16.8f | Res | %16.8f |\n",iter+1, EStep, Ediff);
    }
    if (fabs(Ediff) < this->LandmarkEnergyCriterion)
    {
      fprintf(stdout,"Energy Criterion Met! %.5f\n",EStep);
      break;
    }
    else
    {
      E0 = EStep;
    }
    //if (this->Verbose == 3)
    //{
    //  if (iter%this->NumSaveIterations == 0)
    //  {
    //    if (this->IterOutputFilename != NULL)
    //    {
    //      std::stringstream iterstr;
    //      iterstr << this->SaveIter++;
    //      std::string iterName = this->IterOutputFilename;
    //      std::string filename =  iterName+"_"+iterstr.str()+".vtp";
    //      vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    //        vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    //      writer->SetInputData(this->MappedPd);
    //      writer->SetFileName(filename.c_str());
    //      writer->Write();
    //    }
    //  }
    //}
  }

  fprintf(stdout,"Done with LandmarkMapping...\n");
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkS2LandmarkMatcher::ComputeMeshLandmarkLaplacian(vtkPolyData *pd,
                                                       vtkEdgeTable *edgeTable,
                                                       vtkFloatArray *edgeWeights,
                                                       vtkIntArray *edgeNeighbors,
                                                       vtkFloatArray *laplacian,
                                                       double landmarkWeighting,
                                                       vtkPolyData *targetS2Pd,
                                                       vtkIntArray *isCurrentLandmark,
                                                       vtkIntArray *landmarkPtIds)
{
  int numPts = pd->GetNumberOfPoints();

  for (int i=0; i<numPts; i++)
  {
    double pointLaplacian[3];
    vtkS2LandmarkMatcher::ComputePointLandmarkLaplacian(i, pd, edgeTable, edgeWeights,
                                                        edgeNeighbors, pointLaplacian,
                                                        landmarkWeighting,
                                                        targetS2Pd, isCurrentLandmark,
                                                        landmarkPtIds);
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
int vtkS2LandmarkMatcher::ComputePointLandmarkLaplacian(vtkIdType p0,
                                                        vtkPolyData *pd,
                                                        vtkEdgeTable *edgeTable,
                                                        vtkFloatArray *edgeWeights,
                                                        vtkIntArray *edgeNeighbors,
                                                        double laplacian[],
                                                        double landmarkWeighting,
                                                        vtkPolyData *targetS2Pd,
                                                        vtkIntArray *isCurrentLandmark,
                                                        vtkIntArray *landmarkPtIds)
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

  if(landmarkPtIds != NULL)
  {
    pointNeighbors->InsertNextId(p0);
    for (int i=0; i<pointNeighbors->GetNumberOfIds(); i++)
    {
      vtkIdType pId = pointNeighbors->GetId(i);
      double currentPt[3], targetPt[3];
      pd->GetPoint(pId, currentPt);
      int isLandmark = isCurrentLandmark->GetValue(pId);
      if (isLandmark != -1)
      {
        int targId = landmarkPtIds->GetValue(isLandmark);
        targetS2Pd->GetPoint(targId, targetPt);
      }
      else
      {
        targetPt[0] = currentPt[0];
        targetPt[1] = currentPt[0];
        targetPt[2] = currentPt[0];
      }
      for (int j=0; j<3; j++)
      {
        laplacian[j] += landmarkWeighting * (currentPt[j] - targetPt[j]);
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
int vtkS2LandmarkMatcher::ComputeLandmarkEnergy(vtkPolyData *pd,
                                                vtkEdgeTable *edgeTable,
                                                vtkFloatArray *edgeWeights,
                                                double &landmarkEnergy,
                                                double landmarkWeighting,
                                                vtkPolyData *targetS2Pd,
                                                vtkIntArray *isCurrentLandmark,
                                                vtkIntArray *landmarkPtIds)
{
  landmarkEnergy = 0.0;
  double compEnergy[3];
  compEnergy[0] = 0.0; compEnergy[1] = 0.0; compEnergy[2] =0.0;
  int numEdges = edgeTable->GetNumberOfEdges();

  edgeTable->InitTraversal();
  for (int i=0; i<numEdges; i++)
  {
    vtkIdType p0, p1;
    vtkIdType edgeId = edgeTable->GetNextEdge(p0, p1);
    double weight = edgeWeights->GetValue(edgeId);

    double h0[3];
    double h1[3];
    pd->GetPoint(p0, h0);
    pd->GetPoint(p1, h1);

    //Calculate String Energy!
    double edgeEnergy[3];
    vtkSphericalConformalMapper::ComputeStringEnergy(h0, h1, weight, edgeEnergy);
    for (int j=0; j<3; j++)
    {
      compEnergy[j] += (1.0/2.0)*edgeEnergy[j];
    }
  }

  if (landmarkPtIds != NULL)
  {
    int numPoints = pd->GetNumberOfPoints();
    for (int i=0; i<numPoints; i++)
    {
      double edgeEnergy[3];
      double currentPt[3], targetPt[3];
      pd->GetPoint(i, currentPt);
      int isLandmark = isCurrentLandmark->GetValue(i);
      if (isLandmark != -1)
      {
        int targId = landmarkPtIds->GetValue(isLandmark);
        targetS2Pd->GetPoint(targId, targetPt);
      }
      else
      {
        targetPt[0] = 1.0;
        targetPt[1] = 0.0;
        targetPt[2] = 0.0;
      }
      vtkSphericalConformalMapper::ComputeStringEnergy(currentPt, targetPt, 1.0, edgeEnergy);
      for (int j=0; j<3; j++)
      {
        compEnergy[j] += (landmarkWeighting/2.0)*edgeEnergy[j];
      }
    }
  }
  for (int i=0; i<3; i++)
  {
    landmarkEnergy += compEnergy[i];
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkS2LandmarkMatcher::FirstStep()
{
  int numPts = this->MappedPd->GetNumberOfPoints();
  vtkSmartPointer<vtkFloatArray> laplacian =
    vtkSmartPointer<vtkFloatArray>::New();
  laplacian->SetNumberOfComponents(3);
  laplacian->Allocate(numPts, 10000);
  laplacian->SetNumberOfTuples(numPts);

  if (this->ComputeMeshLandmarkLaplacian(this->MappedPd, this->SourceEdgeTable,
                                         this->SourceEdgeWeights, this->SourceEdgeNeighbors,
                                         laplacian, this->LandmarkWeighting,
                                         this->TargetPd, this->IsSourceLandmark,
                                         this->SourceLandmarkPtIds) != 1)
  {
    vtkErrorMacro("Error when computing laplacian");
    return 0;
  }
  if (this->UpdateMap(laplacian, vtkSphericalConformalMapper::CG_NONE) != 1)
  {
    vtkErrorMacro("Error when updating landmark map");
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
int vtkS2LandmarkMatcher::UpdateMap(vtkFloatArray *laplacian,
                                    int cg_update)
{
  int numPts = this->MappedPd->GetNumberOfPoints();

  for (int i=0; i<numPts; i++)
  {
    double pointLaplacian[3], pointNormal[3];
    laplacian->GetTuple(i, pointLaplacian);
    this->MappedPd->GetPoint(i, pointNormal);

    double pointScalar = vtkMath::Dot(pointLaplacian, pointNormal);
    double pointLaplacianNormal[3];
    double pointLaplacianTangential[3];
    for (int j=0; j<3; j++)
    {
      pointLaplacianNormal[j] = pointScalar * pointNormal[j];
      pointLaplacianTangential[j] = -1.0*(pointLaplacian[j] - pointLaplacianNormal[j]);
    }
    if (cg_update == vtkSphericalConformalMapper::CG_NONE)
    {
      this->PrevDescent->SetTuple(i, pointLaplacianTangential);
    }
    else
    {
      this->CurrDescent->SetTuple(i, pointLaplacianTangential);
    }
  }

  if (this->StepForward(cg_update) != 1)
  {
    vtkErrorMacro("Error when updating landmark map in CG step");
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
int vtkS2LandmarkMatcher::StepForward(int cg_update)
{
  if (cg_update == vtkSphericalConformalMapper::CG_NONE)
  {
    int numPts = this->MappedPd->GetNumberOfPoints();
    for (int i=0; i<numPts; i++)
    {
      double ptVal[3], descent[3];
      this->MappedPd->GetPoint(i, ptVal);
      this->PrevDescent->GetTuple(i, descent);
      this->ConjugateDir->SetTuple(i, descent);
      for (int j=0; j<3; j++)
      {
        ptVal[j] = ptVal[j] + this->LandmarkTimeStep * descent[j];
      }
      vtkMath::Normalize(ptVal);
      this->MappedPd->GetPoints()->SetPoint(i, ptVal);
    }
    return 1;
  }
  else if (cg_update == vtkSphericalConformalMapper::CG_FLETCHER_REEVES)
  {
    this->FRUpdateMap();
  }
  else if (cg_update == vtkSphericalConformalMapper::CG_POLAK_RIBIERE)
  {
    this->PRUpdateMap();
  }
  else if (cg_update == vtkSphericalConformalMapper::CG_HESTENESS_STIEFEL)
  {
    this->HSUpdateMap();
  }
  else if (cg_update == vtkSphericalConformalMapper::CG_DAI_YUAN)
  {
    this->DYUpdateMap();
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
int vtkS2LandmarkMatcher::FRUpdateMap()
{
  int numPts = this->PrevDescent->GetNumberOfTuples();
  double numerator[3], denominator[3];
  vtkSphericalConformalMapper::VectorDotProduct(this->CurrDescent, this->CurrDescent,
                                                  numerator, numPts, 3);
  vtkSphericalConformalMapper::VectorDotProduct(this->PrevDescent, this->PrevDescent,
                                                  denominator, numPts, 3);

  double beta[3];
  for (int i=0; i<3; i++)
  {
    beta[i] = numerator[i]/denominator[i];
  }

  this->CGUpdateMap(beta);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkS2LandmarkMatcher::PRUpdateMap()
{
  int numPts = this->PrevDescent->GetNumberOfTuples();
  double numerator[3], denominator[3];
  vtkSmartPointer<vtkFloatArray> difference =
    vtkSmartPointer<vtkFloatArray>::New();
  difference->SetNumberOfComponents(3);
  difference->Allocate(numPts, 10000);
  difference->SetNumberOfTuples(numPts);
  vtkSphericalConformalMapper::VectorAdd(this->CurrDescent, this->PrevDescent, -1.0,
                                           difference, numPts, 3);
  vtkSphericalConformalMapper::VectorDotProduct(this->CurrDescent, difference,
                                                  numerator, numPts, 3);
  vtkSphericalConformalMapper::VectorDotProduct(this->PrevDescent, this->PrevDescent,
                                                  denominator, numPts, 3);

  double beta[3];
  for (int i=0; i<3; i++)
  {
    beta[i] = std::max(0.0, numerator[i]/denominator[i]);
  }

  this->CGUpdateMap(beta);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkS2LandmarkMatcher::HSUpdateMap()
{
  int numPts = this->PrevDescent->GetNumberOfTuples();
  double numerator[3], denominator[3];
  vtkSmartPointer<vtkFloatArray> difference =
    vtkSmartPointer<vtkFloatArray>::New();
  difference->SetNumberOfComponents(3);
  difference->Allocate(numPts, 10000);
  difference->SetNumberOfTuples(numPts);
  vtkSphericalConformalMapper::VectorAdd(this->CurrDescent, this->PrevDescent, -1.0,
                                           difference, numPts, 3);
  vtkSphericalConformalMapper::VectorDotProduct(this->CurrDescent, difference,
                                                  numerator, numPts, 3);
  vtkSphericalConformalMapper::VectorDotProduct(this->ConjugateDir, difference,
                                                  denominator, numPts, 3);

  double beta[3];
  for (int i=0; i<3; i++)
  {
    beta[i] = -1.0 *numerator[i]/denominator[i];
  }

  this->CGUpdateMap(beta);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkS2LandmarkMatcher::DYUpdateMap()
{
  int numPts = this->PrevDescent->GetNumberOfTuples();
  double numerator[3], denominator[3];
  vtkSmartPointer<vtkFloatArray> difference =
    vtkSmartPointer<vtkFloatArray>::New();
  difference->SetNumberOfComponents(3);
  difference->Allocate(numPts, 10000);
  difference->SetNumberOfTuples(numPts);
  vtkSphericalConformalMapper::VectorAdd(this->CurrDescent, this->PrevDescent, -1.0,
                                           difference, numPts, 3);
  vtkSphericalConformalMapper::VectorDotProduct(this->CurrDescent, this->CurrDescent,
                                                  numerator, numPts, 3);
  vtkSphericalConformalMapper::VectorDotProduct(this->ConjugateDir, difference,
                                                  denominator, numPts, 3);

  double beta[3];
  for (int i=0; i<3; i++)
  {
    beta[i] = -1.0 * numerator[i]/denominator[i];
  }

  this->CGUpdateMap(beta);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkS2LandmarkMatcher::CGUpdateMap(double beta[])
{
  int numPts = this->MappedPd->GetNumberOfPoints();

  double descCond[3];
  vtkSphericalConformalMapper::VectorDotProduct(this->CurrDescent, this->ConjugateDir,
                                                  descCond, numPts, 3);

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

  this->WolfeLineSearch();

  for (int i=0; i<numPts; i++)
  {
    double newDescent[3], ptVal[3];
    this->ConjugateDir->GetTuple(i, newDescent);
    this->MappedPd->GetPoint(i, ptVal);
    for (int j=0; j<3; j++)
    {
      ptVal[j] = ptVal[j] + this->LandmarkTimeStep * newDescent[j];
    }
    vtkMath::Normalize(ptVal);
    this->MappedPd->GetPoints()->SetPoint(i, ptVal);
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkS2LandmarkMatcher::WolfeLineSearch()
{
  int numPts = this->MappedPd->GetNumberOfPoints();

  vtkSmartPointer<vtkFloatArray> conjLaplacian =
    vtkSmartPointer<vtkFloatArray>::New();
  conjLaplacian->SetNumberOfComponents(3);
  conjLaplacian->Allocate(numPts, 10000);
  conjLaplacian->SetNumberOfTuples(numPts);
  vtkSphericalConformalMapper::ComputeDataArrayLaplacian(this->ConjugateDir,
                                                           this->MappedPd,
                                                           this->SourceEdgeTable,
                                                           this->SourceEdgeWeights,
                                                           this->SourceEdgeNeighbors,
                                                           conjLaplacian,
                                                           vtkSphericalConformalMapper::HARMONIC);

  double numerator[3];
  double denominator[3];
  vtkSphericalConformalMapper::VectorDotProduct(this->ConjugateDir,
                                                  this->CurrDescent, numerator,
                                                  numPts, 3);
  vtkSphericalConformalMapper::VectorDotProduct(this->ConjugateDir,
                                                  conjLaplacian, denominator,
                                                  numPts, 3);

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

  this->LandmarkTimeStep = maxstep;

  return 1;
}
