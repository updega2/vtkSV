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

/** @file vtkSurfaceToSurfaceMapper.cxx
 *  @brief This implements the vtkSurfaceToSurfaceMapper filter as a class
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkSurfaceToSurfaceMapper.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkCenterOfMass.h"
#include "vtkCleanPolyData.h"
#include "vtkFloatArray.h"
#include "vtkGenericCell.h"
#include "vtkGradientFilter.h"
#include "vtkS2LandmarkMatcher.h"
#include "vtkLoopSubdivisionFilter.h"
#include "vtkMath.h"
#include "vtkMapInterpolator.h"
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
vtkCxxRevisionMacro(vtkSurfaceToSurfaceMapper, "$Revision: 0.0 $");
vtkStandardNewMacro(vtkSurfaceToSurfaceMapper);


//---------------------------------------------------------------------------
vtkSurfaceToSurfaceMapper::vtkSurfaceToSurfaceMapper()
{
  this->SetNumberOfInputPorts(2);

  this->Verbose                 = 1;
  this->TutteTimeStep           = 0.001;
  this->ConformalTimeStep       = 0.00001;
  this->TutteEnergyCriterion    = 0.000001;
  this->HarmonicEnergyCriterion = 0.000001;
  this->MaxNumIterations        = 1000;
  this->NumSourceSubdivisions   = 0;
  this->NumLandmarks            = 0;

  this->SourcePd   = vtkPolyData::New();
  this->SourceS2Pd = vtkPolyData::New();
  this->TargetPd   = vtkPolyData::New();
  this->TargetS2Pd = vtkPolyData::New();
  this->MappedPd   = vtkPolyData::New();
  this->MappedS2Pd = vtkPolyData::New();

  this->SourceLandmarkPtIds = NULL;
  this->SourceLandmarkProj  = vtkFloatArray::New();
  this->TargetLandmarkPtIds = NULL;
  this->TargetLandmarkProj  = vtkFloatArray::New();
}

//---------------------------------------------------------------------------
vtkSurfaceToSurfaceMapper::~vtkSurfaceToSurfaceMapper()
{
  if (this->SourcePd != NULL)
  {
    this->SourcePd->Delete();
  }
  if (this->SourceS2Pd != NULL)
  {
    this->SourceS2Pd->Delete();
  }
  if (this->TargetPd != NULL)
  {
    this->TargetPd->Delete();
  }
  if (this->TargetS2Pd != NULL)
  {
    this->TargetS2Pd->Delete();
  }
  if (this->MappedPd != NULL)
  {
    this->MappedPd->Delete();
  }
  if (this->MappedS2Pd != NULL)
  {
    this->MappedS2Pd->Delete();
  }
  if (this->SourceLandmarkProj != NULL)
  {
    this->SourceLandmarkProj->Delete();
  }
  if (this->TargetLandmarkProj != NULL)
  {
    this->TargetLandmarkProj->Delete();
  }
}

//---------------------------------------------------------------------------
void vtkSurfaceToSurfaceMapper::PrintSelf(ostream& os, vtkIndent indent)
{
}

// Generate Separated Surfaces with Region ID Numbers
//---------------------------------------------------------------------------
int vtkSurfaceToSurfaceMapper::RequestData(
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

  if (this->RunSphericalConformalMappings() != 1)
  {
    vtkErrorMacro("Error during mapping");
    return 0;
  }

  //if (this->RunLandmarkMatching() != 1)
  //{
  //  vtkErrorMacro("Error during landmark matching");
  //  return 0;
  //}

  //if (this->SubdivideAndInterpolate() != 1)
  //{
  //  return 0;
  //}

  if (this->StereographicProjection(this->SourceS2Pd) != 1)
  {
    vtkErrorMacro("Error computing stereographic projection of source pd");
    return 0;
  }
  vtkSphericalConformalMapper::ComputeNormals(this->SourceS2Pd);
  this->ConvertValueToPolyData(this->SourceS2Pd, "ComplexValue",
                               this->MappedS2Pd);

  //if (this->StereographicProjection(this->TargetS2Pd) != 1)
  //{
  //  vtkErrorMacro("Error computing stereographic projection of source pd");
  //  return 0;
  //}
  //if (this->RunMobiusTransformation() != 1)
  //{
  //  vtkErrorMacro("Error running mobius transform");
  //  return 0;
  //}
  //if (this->InverseStereographicProjection(this->TargetPd) != 1)
  //{
  //  vtkErrorMacro("Inverse stereographic projection failed");
  //  return 0;
  //}

  output->DeepCopy(this->MappedS2Pd);
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSurfaceToSurfaceMapper::RunSphericalConformalMappings()
{
  vtkSmartPointer<vtkFloatArray> sourceMap =
    vtkSmartPointer<vtkFloatArray>::New();
  sourceMap->SetNumberOfComponents(3);
  sourceMap->Allocate(this->SourcePd->GetNumberOfPoints(), 10000);
  sourceMap->SetNumberOfTuples(this->SourcePd->GetNumberOfPoints());
  sourceMap->SetName("SourceMap");
  this->GetMap(this->SourcePd, this->SourceS2Pd, sourceMap);
  this->SourceS2Pd->GetPointData()->AddArray(sourceMap);

  vtkSmartPointer<vtkFloatArray> targetMap =
    vtkSmartPointer<vtkFloatArray>::New();
  targetMap->SetNumberOfComponents(3);
  targetMap->Allocate(this->TargetPd->GetNumberOfPoints(), 10000);
  targetMap->SetNumberOfTuples(this->TargetPd->GetNumberOfPoints());
  targetMap->SetName("TargetMap");
  this->GetMap(this->TargetPd, this->TargetS2Pd, targetMap);
  this->TargetS2Pd->GetPointData()->AddArray(targetMap);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSurfaceToSurfaceMapper::RunLandmarkMatching()
{
  vtkSmartPointer<vtkS2LandmarkMatcher> matcher =
    vtkSmartPointer<vtkS2LandmarkMatcher>::New();
  matcher->SetInputData(0, this->SourcePd);
  matcher->SetInputData(1, this->TargetPd);
  matcher->SetVerbose(this->Verbose);
  matcher->SetLandmarkTimeStep(this->TutteTimeStep);
  matcher->SetLandmarkEnergyCriterion(this->HarmonicEnergyCriterion);
  matcher->SetMaxNumIterations(this->MaxNumIterations);
  matcher->SetCGUpdateMethod(vtkSphericalConformalMapper::CG_NONE);
  if (this->SourceLandmarkPtIds != NULL && this->TargetLandmarkPtIds != NULL)
  {
    matcher->SetSourceLandmarkPtIds(this->SourceLandmarkPtIds);
    matcher->SetTargetLandmarkPtIds(this->TargetLandmarkPtIds);
  }
  matcher->Update();

  this->MappedS2Pd->DeepCopy(matcher->GetOutput());
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSurfaceToSurfaceMapper::SubdivideAndInterpolate()
{
  vtkSmartPointer<vtkMapInterpolator> mapInterpolator =
    vtkSmartPointer<vtkMapInterpolator>::New();
  mapInterpolator->SetVerbose(this->Verbose);
  mapInterpolator->SetInputData(0, this->SourcePd);
  mapInterpolator->SetInputData(1, this->TargetPd);
  mapInterpolator->SetInputData(2, this->SourceS2Pd);
  mapInterpolator->SetInputData(3, this->TargetS2Pd);
  mapInterpolator->SetNumSourceSubdivisions(this->NumSourceSubdivisions);
  mapInterpolator->Update();

  this->MappedPd->DeepCopy(mapInterpolator->GetOutput());

  return 1;
}
//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSurfaceToSurfaceMapper::GetMap(vtkPolyData *pd, vtkPolyData *spherePd,
                                      vtkFloatArray *map)
{
  vtkSmartPointer<vtkSphericalConformalMapper> mapper =
    vtkSmartPointer<vtkSphericalConformalMapper>::New();
  mapper->SetInputData(pd);
  mapper->SetVerbose(this->Verbose);
  mapper->SetInitialTimeStep(this->TutteTimeStep);
  mapper->SetTutteEnergyCriterion(this->TutteEnergyCriterion);
  mapper->SetHarmonicEnergyCriterion(this->HarmonicEnergyCriterion);
  mapper->SetMaxNumIterations(this->MaxNumIterations);
  mapper->SetCGUpdateMethod(vtkSphericalConformalMapper::CG_DAI_YUAN);
  mapper->Update();

  spherePd->DeepCopy(mapper->GetOutput(0));

  int numPts = pd->GetNumberOfPoints();
  for (int i=0; i< numPts; i++)
  {
    double originalPt[3];
    double mappedPt[3];
    pd->GetPoint(i, originalPt);
    spherePd->GetPoint(i, mappedPt);

    double mapVector[3];
    for (int j=0; j<3; j++)
    {
      mapVector[j] = originalPt[j] - mappedPt[j];
    }
    map->InsertTuple(i, mapVector);
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSurfaceToSurfaceMapper::StereographicProjection(vtkPolyData *pd)
{
  int numPts = pd->GetNumberOfPoints();

  vtkSmartPointer<vtkFloatArray> complexValue =
    vtkSmartPointer<vtkFloatArray>::New();
  complexValue->SetNumberOfComponents(2);
  complexValue->Allocate(numPts, 10000);
  complexValue->SetNumberOfTuples(numPts);
  complexValue->SetName("ComplexValue");

  for (int i=0; i<numPts; i++)
  {
    double pt[3];
    pd->GetPoint(i, pt);
    double x = pt[0] / (1.0 - pt[2]);
    double y = pt[1] / (1.0 - pt[2]);
    if (0.9999 <= pt[2] && pt[2] <= 1.0001)
    {
      x = 0.0;
      y = 0.0;
    }
    complexValue->SetComponent(i, 0, x);
    complexValue->SetComponent(i, 1, y);
  }

  pd->GetPointData()->AddArray(complexValue);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSurfaceToSurfaceMapper::InverseStereographicProjection(vtkPolyData *pd)
{
  int numPts = pd->GetNumberOfPoints();

  vtkFloatArray *complexValue = vtkFloatArray::SafeDownCast(
    pd->GetPointData()->GetArray("ComplexValue"));

  for (int i=0; i<numPts; i++)
  {
    double cv[2];
    complexValue->GetTuple(i, cv);
    double r2 = pow(cv[0], 2.0) + pow(cv[1], 2.0);
    double pt[3];
    pt[0] = 2.0*cv[0] / (1.0 + r2);
    pt[1] = 2.0*cv[1] / (1.0 + r2);
    pt[2] = (-1.0 + r2) / (1.0 + r2);
    pd->GetPoints()->SetPoint(i, pt);
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSurfaceToSurfaceMapper::RunMobiusTransformation()
{
  if (this->SetLandmarks() != 1)
  {
    vtkErrorMacro("Error setting landmarks on mapped surfaces");
    return 0;
  }

  std::complex<double> a(0.0,0.0), b(0.0,0.0);
  if (this->ComputeComplexMobius(a, b) != 1)
  {
    vtkErrorMacro("Mobuis transform calculation did not go well");
    return 0;
  }

  if (this->TransformComplexValue(a, b) != 1)
  {
    vtkErrorMacro("Transformation did not go well");
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
int vtkSurfaceToSurfaceMapper::SetLandmarks()
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
    this->SourceLandmarkProj->SetNumberOfComponents(2);
    this->SourceLandmarkProj->Allocate(this->NumLandmarks, 10000);
    this->SourceLandmarkProj->SetNumberOfTuples(this->NumLandmarks);
    this->TargetLandmarkProj->SetNumberOfComponents(2);
    this->TargetLandmarkProj->Allocate(this->NumLandmarks, 10000);
    this->TargetLandmarkProj->SetNumberOfTuples(this->NumLandmarks);
    vtkFloatArray *srcCV = vtkFloatArray::SafeDownCast(
      this->SourcePd->GetPointData()->GetArray("ComplexValue"));
    vtkFloatArray *targCV = vtkFloatArray::SafeDownCast(
      this->TargetPd->GetPointData()->GetArray("ComplexValue"));

    for (int i =0; i<this->NumLandmarks; i++)
    {
      double srcPt[2], targPt[2];
      int srcId = this->SourceLandmarkPtIds->GetValue(i);
      srcCV->GetTuple(srcId, srcPt);
      this->SourceLandmarkProj->SetTuple(i, srcPt);
      int targId = this->TargetLandmarkPtIds->GetValue(i);
      targCV->GetTuple(targId, targPt);
      this->TargetLandmarkProj->SetTuple(i, targPt);
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
int vtkSurfaceToSurfaceMapper::ComputeConformalFactor(std::complex<double> z,
                                                      std::complex<double> &g)
{
  std::complex<double> zconj = std::conj(z);
  g = 4.0 / (1.0 + z * zconj);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSurfaceToSurfaceMapper::ComputeComplexMobius(std::complex<double> &a,
                                                    std::complex<double> &b)
{
  std::vector<std::complex<double > > lhs;
  std::vector<std::complex<double > > rhs;

  this->GetSystem(lhs, rhs);

  double memMat[4][4], memInvMat[4][4];
  double *mat[4], *invMat[4];
  double memCol[4];

  // | real[0] -imag[0] real[1] -imag[1] |   | a  |
  // | imag[0]  real[0] imag[1]  real[1] |   | bi |
  // | real[2] -imag[2] real[3] -imag[3] |   | c  |
  // | imag[2]  real[2] imag[3]  real[3] |   | di |
  for (int i=0; i<2; i++)
  {
    for (int j=0; j<2; j++)
    {
      memMat[2*i][2*j]     = (lhs[2*i+j]).real();
      memMat[2*i][2*j+1]   = -1.0 * lhs[2*i+j].imag();
      memMat[2*i+1][2*j]   = lhs[2*i+j].imag();
      memMat[2*i+1][2*j+1] = lhs[2*i+j].real();
    }
    memCol[2*i] = rhs[i].real();
    memCol[2*i+1] = rhs[i].imag();
  }

  for (int i=0; i<4; i++)
  {
    mat[i] = &(memMat[i][0]);
    invMat[i] = &(memInvMat[i][0]);
  }
  fprintf(stdout, "Matrix is\n");
  fprintf(stdout, "| %.8f %.8f %.8f %.8f |\n", memMat[0][0], memMat[0][1], memMat[0][2], memMat[0][3]);
  fprintf(stdout, "| %.8f %.8f %.8f %.8f |\n", memMat[1][0], memMat[1][1], memMat[1][2], memMat[1][3]);
  fprintf(stdout, "| %.8f %.8f %.8f %.8f |\n", memMat[2][0], memMat[2][1], memMat[2][2], memMat[2][3]);
  fprintf(stdout, "| %.8f %.8f %.8f %.8f |\n", memMat[3][0], memMat[3][1], memMat[3][2], memMat[3][3]);

  fprintf(stdout, "Right hand side is:\n");
  fprintf(stdout, " %.8f |\n", memCol[0]);
  fprintf(stdout, " %.8f |\n", memCol[1]);
  fprintf(stdout, " %.8f |\n", memCol[2]);
  fprintf(stdout, " %.8f |\n", memCol[3]);

  if (vtkMath::InvertMatrix(mat, invMat, 4) == 0)
  {
    vtkErrorMacro("Linear system solve failed");
    return 0;
  }
  fprintf(stdout, "Inverse is\n");
  fprintf(stdout, "| %.8f %.8f %.8f %.8f |\n", memInvMat[0][0], memInvMat[0][1], memInvMat[0][2], memInvMat[0][3]);
  fprintf(stdout, "| %.8f %.8f %.8f %.8f |\n", memInvMat[1][0], memInvMat[1][1], memInvMat[1][2], memInvMat[1][3]);
  fprintf(stdout, "| %.8f %.8f %.8f %.8f |\n", memInvMat[2][0], memInvMat[2][1], memInvMat[2][2], memInvMat[2][3]);
  fprintf(stdout, "| %.8f %.8f %.8f %.8f |\n", memInvMat[3][0], memInvMat[3][1], memInvMat[3][2], memInvMat[3][3]);


  double newA = 0.0;
  double newB = 0.0;
  double newC = 0.0;
  double newD = 0.0;
  for (int i=0; i<4; i++)
  {
    for (int j=0; j<4; j++)
    {
      newA += memInvMat[i][j] * memCol[j];
      newB += memInvMat[i+1][j] * memCol[j];
      newC += memInvMat[i+2][j] * memCol[j];
      newD += memInvMat[i+3][j] * memCol[j];
    }
  }
  std::complex<double> tmpA(newA, newB);
  std::complex<double> tmpB(newC, newD);
  a = tmpA;
  b = tmpB;
  fprintf(stdout, "What is a: %.8f + %.8fi\n", newA, newB);
  fprintf(stdout, "What is b: %.8f + %.8fi\n", newC, newD);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSurfaceToSurfaceMapper::TransformComplexValue(std::complex<double> a,
                                                     std::complex<double> b)
{
  int numPts = this->TargetPd->GetNumberOfPoints();
  vtkFloatArray *complexValue = vtkFloatArray::SafeDownCast(
    this->TargetPd->GetPointData()->GetArray("ComplexValue"));
  for (int i=0; i<numPts; i++)
  {
    double cv[2];
    complexValue->GetTuple(i, cv);
    std::complex<double> tmp(cv[0], cv[1]);
    double theta = 3.14/4.0;
    std::complex<double> duh(0.0,1);
    a = std::exp(duh * theta /2.0);
    b = 0.0;
    std::complex<double> c(0.0, 0.0);
    std::complex<double> d = std::exp( -1.0 *duh * theta / 2.0);
    std::complex<double> newTmp = (a*tmp + b)/ (c*tmp + d);
    double newCV[2];
    newCV[0] = newTmp.real();
    newCV[1] = newTmp.imag();

    complexValue->SetTuple(i, newCV);
  }

  this->TargetPd->GetPointData()->AddArray(complexValue);
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSurfaceToSurfaceMapper::GetSystem(std::vector<std::complex<double> > &lhs,
                                         std::vector<std::complex<double> > &rhs)
{

  lhs.clear();
  for (int i=0; i<4; i++)
  {
    std::complex<double> zero(0.0, 0.0);
    lhs.push_back( zero);
  }
  rhs.clear();
  for (int i=0; i<2; i++)
  {
    std::complex<double> zero(0.0, 0.0);
    rhs.push_back( zero);
  }

  for (int i=0; i<this->NumLandmarks; i++)
  {
    double srcProj[2], targProj[2];
    this->SourceLandmarkProj->GetTuple(i, srcProj);
    this->TargetLandmarkProj->GetTuple(i, targProj);
    fprintf(stdout,"Source Proj Val: %.8f %.8f\n", srcProj[0], srcProj[1]);
    fprintf(stdout,"Target Proj Val: %.8f %.8f\n", targProj[0], targProj[1]);
    std::complex<double> z(srcProj[0], srcProj[1]);
    std::complex<double> t(targProj[0], targProj[1]);
    std::complex<double> g(0.0, 0.0);
    this->ComputeConformalFactor(z, g);

    fprintf(stdout,"z: %.8f + %.8fi\n", z.real(), z.imag());
    fprintf(stdout,"g: %.8f + %.8fi\n", g.real(), g.imag());
    fprintf(stdout,"t: %.8f + %.8fi\n", t.real(), t.imag());
    lhs[0] += g * z * z;
    lhs[1] += g * z;
    lhs[2] += g * z;
    lhs[3] += g;

    rhs[0] += g * z * t;
    rhs[1] += g * t;
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSurfaceToSurfaceMapper::ConvertValueToPolyData(vtkPolyData *inPd,
                                                      std::string fieldName,
                                                      vtkPolyData *outPd)
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
      double pt[2];
      fieldArray->GetTuple(pts[j], pt);

      double newPt[3];
      for (int k=0; k<2; k++)
      {
        newPt[k] = pt[k];
      }
      newPt[2] = 0.0;
      fieldPts->SetPoint(pts[j], newPt);
    }
  }

  outPd->SetPolys(fieldCells);
  outPd->SetPoints(fieldPts);
  outPd->BuildLinks();

  return 1;
}


