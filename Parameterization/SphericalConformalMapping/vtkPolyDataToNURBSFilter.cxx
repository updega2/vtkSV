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

/** @file vtkPolyDataToNURBSFilter.cxx
 *  @brief This implements the vtkPolyDataToNURBSFilter filter as a class
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkPolyDataToNURBSFilter.h"

#include "vtkAppendPolyData.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCenterOfMass.h"
#include "vtkCleanPolyData.h"
#include "vtkClipPolyData.h"
#include "vtkConnectivityFilter.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkDoubleArray.h"
#include "vtkExtractGeometry.h"
#include "vtkFloatArray.h"
#include "vtkIdFilter.h"
#include "vtkIntArray.h"
#include "vtkGradientFilter.h"
#include "vtkMapInterpolator.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPlane.h"
#include "vtkPointData.h"
#include "vtkPointDataToCellData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkPolyDataSliceAndDiceFilter.h"
#include "vtkSmartPointer.h"
#include "vtkSphericalConformalMapper.h"
#include "vtkTextureMapToSphere.h"
#include "vtkThreshold.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTriangle.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLPolyDataWriter.h"

#define vtkNew(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

#include <iostream>
#include <cmath>

//---------------------------------------------------------------------------
//vtkCxxRevisionMacro(vtkPolyDataToNURBSFilter, "$Revision: 0.0 $");
vtkStandardNewMacro(vtkPolyDataToNURBSFilter);


//---------------------------------------------------------------------------
vtkPolyDataToNURBSFilter::vtkPolyDataToNURBSFilter()
{
  this->Verbose = 1;

  this->InputPd         = vtkPolyData::New();
  this->ParameterizedPd = vtkPolyData::New();
  this->Centerlines     = NULL;
  this->CubeS2Pd        = NULL;
  this->OpenCubeS2Pd    = NULL;
  this->Polycube        = vtkGeneralizedPolycube::New();

  this->BoundaryPointsArrayName = NULL;
  this->GroupIdsArrayName       = NULL;
  this->SegmentIdsArrayName     = NULL;
  this->SliceIdsArrayName       = NULL;
  this->SphereRadiusArrayName   = NULL;
  this->InternalIdsArrayName    = NULL;
  this->DijkstraArrayName       = NULL;
}

//---------------------------------------------------------------------------
vtkPolyDataToNURBSFilter::~vtkPolyDataToNURBSFilter()
{
  if (this->InputPd != NULL)
  {
    this->InputPd->Delete();
    this->InputPd = NULL;
  }
  if (this->CubeS2Pd != NULL)
  {
    this->CubeS2Pd->Delete();
    this->CubeS2Pd = NULL;
  }
  if (this->OpenCubeS2Pd != NULL)
  {
    this->OpenCubeS2Pd->Delete();
    this->OpenCubeS2Pd = NULL;
  }
  if (this->ParameterizedPd != NULL)
  {
    this->ParameterizedPd->Delete();
    this->ParameterizedPd = NULL;
  }
  if (this->Centerlines != NULL)
  {
    this->Centerlines->Delete();
    this->Centerlines = NULL;
  }
  if (this->Polycube != NULL)
  {
    this->Polycube->Delete();
    this->Polycube = NULL;
  }
  if (this->BoundaryPointsArrayName != NULL)
  {
    delete [] this->BoundaryPointsArrayName;
    this->BoundaryPointsArrayName = NULL;
  }
  if (this->GroupIdsArrayName != NULL)
  {
    delete [] this->GroupIdsArrayName;
    this->GroupIdsArrayName = NULL;
  }
  if (this->SegmentIdsArrayName != NULL)
  {
    delete [] this->SegmentIdsArrayName;
    this->SegmentIdsArrayName = NULL;
  }
  if (this->SliceIdsArrayName != NULL)
  {
    delete [] this->SliceIdsArrayName;
    this->SliceIdsArrayName = NULL;
  }
  if (this->SphereRadiusArrayName != NULL)
  {
    delete [] this->SphereRadiusArrayName;
    this->SphereRadiusArrayName = NULL;
  }
  if (this->InternalIdsArrayName != NULL)
  {
    delete [] this->InternalIdsArrayName;
    this->InternalIdsArrayName = NULL;
  }
  if (this->DijkstraArrayName != NULL)
  {
    delete [] this->DijkstraArrayName;
    this->DijkstraArrayName = NULL;
  }
}

//---------------------------------------------------------------------------
void vtkPolyDataToNURBSFilter::PrintSelf(ostream& os, vtkIndent indent)
{
}

// Generate Separated Surfaces with Region ID Numbers
//---------------------------------------------------------------------------
int vtkPolyDataToNURBSFilter::RequestData(
                                 vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector)
{
  // get the input and output
  vtkPolyData *input1 = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  //Copy the input to operate on
  this->InputPd->DeepCopy(input1);

  if (this->Centerlines == NULL)
  {
    this->Centerlines = vtkPolyData::New();
    this->ComputeCenterlines();
    this->ExtractBranches();
  }

  if (this->SliceAndDice() != 1)
  {
    vtkErrorMacro("Error in slicing polydata\n");
    return 0;
  }

  if (this->PerformMappings() != 1)
  {
    vtkErrorMacro("Error in perform mappings\n");
    return 0;
  }

  output->DeepCopy(this->ParameterizedPd);
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataToNURBSFilter::ComputeCenterlines()
{
  return 0;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataToNURBSFilter::ExtractBranches()
{
  return 0;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataToNURBSFilter::SliceAndDice()
{
  vtkNew(vtkPolyDataSliceAndDiceFilter, slicer);
  slicer->SetInputData(this->InputPd);
  slicer->SetCenterlines(this->Centerlines);
  slicer->SetVerbose(0);
  slicer->SetSliceLength(1.5);
  slicer->SetConstructPolycube(1);
  slicer->SetBoundaryPointsArrayName(this->BoundaryPointsArrayName);
  slicer->SetGroupIdsArrayName(this->GroupIdsArrayName);
  slicer->SetSegmentIdsArrayName(this->SegmentIdsArrayName);
  slicer->SetSliceIdsArrayName(this->SliceIdsArrayName);
  slicer->SetSphereRadiusArrayName(this->SphereRadiusArrayName);
  slicer->SetInternalIdsArrayName(this->InternalIdsArrayName);
  slicer->SetDijkstraArrayName(this->DijkstraArrayName);
  slicer->Update();

  this->InputPd->DeepCopy(slicer->GetOutput());
  this->Polycube->DeepCopy(slicer->GetPolycube());
  this->StartPtIds = vtkIntArray::SafeDownCast(
    this->Polycube->GetCellData()->GetArray("CornerPtIds"));
  this->SliceRightNormals = vtkDoubleArray::SafeDownCast(
    this->Polycube->GetCellData()->GetArray("RightNormal"));
  this->SliceTopNormals = vtkDoubleArray::SafeDownCast(
    this->Polycube->GetCellData()->GetArray("TopNormal"));

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataToNURBSFilter::PerformMappings()
{
  vtkDataArray *segmentIds = this->InputPd->GetCellData()->GetArray(this->SegmentIdsArrayName);
  vtkDataArray *groupIds = this->InputPd->GetCellData()->GetArray(this->GroupIdsArrayName);
  int numIds = segmentIds->GetNumberOfTuples();
  vtkNew(vtkIdFilter, ider);
  ider->SetInputData(this->InputPd);
  ider->SetIdsArrayName(this->InternalIdsArrayName);
  ider->Update();
  this->InputPd->DeepCopy(ider->GetOutput());

  vtkNew(vtkAppendPolyData, appender);
  double segminmax[2], groupminmax[2];
  segmentIds->GetRange(segminmax);
  groupIds->GetRange(groupminmax);
  for (int i=segminmax[0]; i<=segminmax[1]; i++)
  {
    vtkNew(vtkPolyData, branchPd);
    this->GetBranch(i, branchPd);
    vtkDataArray *sliceIds = branchPd->GetCellData()->GetArray(this->SliceIdsArrayName);

    int numPoints = branchPd->GetNumberOfPoints();

    if (numPoints != 0)
    {
      if (i <= groupminmax[1])
      {
        this->MapBranch(branchPd, sliceIds, appender);
      }
      else
      {
        this->MapBifurcation(branchPd, sliceIds, appender);
      }
    }
  }
  appender->Update();
  this->ParameterizedPd->DeepCopy(appender->GetOutput());

  this->InputPd->GetCellData()->RemoveArray(this->InternalIdsArrayName);
  this->InputPd->GetPointData()->RemoveArray(this->InternalIdsArrayName);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataToNURBSFilter::GetBranch(const int branchId, vtkPolyData *branchPd)
{
  vtkNew(vtkThreshold, thresholder);
  thresholder->SetInputData(this->InputPd);
  //Set Input Array to 0 port,0 connection,1 for Cell Data, and Regions is the type name
  thresholder->SetInputArrayToProcess(0, 0, 0, 1, this->SegmentIdsArrayName);
  thresholder->ThresholdBetween(branchId, branchId);
  thresholder->Update();

  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(thresholder->GetOutput());
  surfacer->Update();

  branchPd->DeepCopy(surfacer->GetOutput());

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataToNURBSFilter::GetSlice(const int sliceId, vtkPolyData *branchPd, vtkPolyData *slicePd)
{
  vtkNew(vtkThreshold, thresholder);
  thresholder->SetInputData(branchPd);
  //Set Input Array to 0 port,0 connection,1 for Cell Data, and Regions is the type name
  thresholder->SetInputArrayToProcess(0, 0, 0, 1, this->SliceIdsArrayName);
  thresholder->ThresholdBetween(sliceId, sliceId);
  thresholder->Update();

  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(thresholder->GetOutput());
  surfacer->Update();

  slicePd->DeepCopy(surfacer->GetOutput());

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataToNURBSFilter::MapBranch(vtkPolyData *branchPd,
                                        vtkDataArray *sliceIds,
                                        vtkAppendPolyData *appender)
{
  double minmax[2];
  sliceIds->GetRange(minmax);

  for (int i=minmax[0]; i<=minmax[1]; i++)
  {
    vtkNew(vtkPolyData, slicePd);
    this->GetSlice(i, branchPd, slicePd);
    int numPoints = slicePd->GetNumberOfPoints();
    vtkDataArray *ptIds = slicePd->GetPointData()->GetArray(this->InternalIdsArrayName);

    if (numPoints != 0)
    {
      double xvec[3], zvec[3];
      vtkNew(vtkIntArray, firstLoopPts);
      vtkNew(vtkIntArray, secondLoopPts);
      //fprintf(stdout,"Double Check 1: %.2f\n", this->StartPtIds->GetComponent(i, 0));
      //fprintf(stdout,"Double Check 2: %.2f\n", this->StartPtIds->GetComponent(i, 4));
      for (int j=0; j<3; j++)
      {
        firstLoopPts->InsertNextValue(ptIds->LookupValue(this->StartPtIds->GetComponent(i, j+1)));
        //fprintf(stdout,"Double Check 1: %.2f\n", this->StartPtIds->GetComponent(i, j+1));
        secondLoopPts->InsertNextValue(ptIds->LookupValue(this->StartPtIds->GetComponent(i, j+5)));
        //fprintf(stdout,"Double Check 2: %.2f\n", this->StartPtIds->GetComponent(i, j+5));
      }
      firstLoopPts->InsertNextValue(ptIds->LookupValue(this->StartPtIds->GetComponent(i, 0)));
      secondLoopPts->InsertNextValue(ptIds->LookupValue(this->StartPtIds->GetComponent(i, 4)));
      this->SliceRightNormals->GetTuple(i, xvec);
      this->SliceTopNormals->GetTuple(i, zvec);
      vtkNew(vtkPolyData, sliceS2Pd);
      vtkNew(vtkPolyData, mappedPd);
      fprintf(stdout,"Mapping region %d...\n", i);

      this->MapSliceToS2(slicePd, sliceS2Pd, firstLoopPts, secondLoopPts, xvec, zvec);
      //this->GetCorrespondingCube(cubeS2Pd, boundary);
      this->InterpolateMapOntoTarget(this->CubeS2Pd, slicePd, sliceS2Pd, mappedPd);
      fprintf(stdout,"Done with mapping...\n");
      appender->AddInputData(mappedPd);


      std::stringstream out;
      out << i;

      std::string filename = "/Users/adamupdegrove/Desktop/tmp/S2Slice_"+out.str()+".vtp";
      vtkNew(vtkXMLPolyDataWriter, writer);
      writer->SetInputData(sliceS2Pd);
      writer->SetFileName(filename.c_str());
      writer->Write();

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
int vtkPolyDataToNURBSFilter::MapBifurcation(vtkPolyData *bifurcationPd,
                                             vtkDataArray *sliceIds,
                                             vtkAppendPolyData *appender)
{
  double minmax[2];
  sliceIds->GetRange(minmax);

  for (int i=minmax[0]; i<=minmax[1]; i++)
  {
    vtkNew(vtkPolyData, slicePd);
    this->GetSlice(i, bifurcationPd, slicePd);
    int numPoints = slicePd->GetNumberOfPoints();
    vtkDataArray *ptIds = slicePd->GetPointData()->GetArray(this->InternalIdsArrayName);

    if (numPoints != 0)
    {
      double xvec[3], zvec[3];
      vtkNew(vtkIntArray, firstLoopPts);
      //fprintf(stdout,"Double Check 1: %.2f\n", this->StartPtIds->GetComponent(i, 0));
      for (int j=0; j<7; j++)
      {
        firstLoopPts->InsertNextValue(ptIds->LookupValue(this->StartPtIds->GetComponent(i, j+1)));
        //fprintf(stdout,"Double Check 1: %.2f\n", this->StartPtIds->GetComponent(i, j+1));
      }
      firstLoopPts->InsertNextValue(ptIds->LookupValue(this->StartPtIds->GetComponent(i, 0)));
      this->SliceRightNormals->GetTuple(i, xvec);
      this->SliceTopNormals->GetTuple(i, zvec);
      vtkNew(vtkPolyData, sliceS2Pd);
      vtkNew(vtkPolyData, mappedPd);
      fprintf(stdout,"Mapping region %d...\n", i);

      this->MapOpenSliceToS2(bifurcationPd, sliceS2Pd, firstLoopPts, xvec, zvec);
      //this->GetCorrespondingCube(cubeS2Pd, boundary);
      this->InterpolateMapOntoTarget(this->OpenCubeS2Pd, bifurcationPd, sliceS2Pd, mappedPd);
      fprintf(stdout,"Done with mapping...\n");
      appender->AddInputData(mappedPd);

      std::stringstream out;
      out << i;

      std::string filename = "/Users/adamupdegrove/Desktop/tmp/S2Slice_"+out.str()+".vtp";
      vtkNew(vtkXMLPolyDataWriter, writer);
      writer->SetInputData(sliceS2Pd);
      writer->SetFileName(filename.c_str());
      writer->Write();

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
int vtkPolyDataToNURBSFilter::MapSliceToS2(vtkPolyData *slicePd,
                                       vtkPolyData *sliceS2Pd,
                                       vtkIntArray *firstCorners,
                                       vtkIntArray *secondCorners,
                                       double xvec[3],
                                       double zvec[3])
{
  int starts[2];
  starts[0] = firstCorners->GetValue(3);
  starts[1] = secondCorners->GetValue(3);
  int cubeStart[2]; cubeStart[0] = 0; cubeStart[1] = 4;;
  //fprintf(stdout,"Four diff values are: %d %d %d %d\n", firstCorners->GetValue(0),
  //                                                      firstCorners->GetValue(1),
  //                                                      firstCorners->GetValue(2),
  //                                                      firstCorners->GetValue(3));
  //fprintf(stdout,"Other values are: %d %d %d %d\n", secondCorners->GetValue(0),
  //                                                  secondCorners->GetValue(1),
  //                                                  secondCorners->GetValue(2),
  //                                                  secondCorners->GetValue(3));

  vtkNew(vtkIntArray, firstLoopHelper);
  vtkNew(vtkIntArray, secondLoopHelper);
  firstLoopHelper->InsertNextValue(0);
  firstLoopHelper->InsertNextValue(1);
  firstLoopHelper->InsertNextValue(3);
  firstLoopHelper->InsertNextValue(4);
  secondLoopHelper->InsertNextValue(0);
  secondLoopHelper->InsertNextValue(1);
  secondLoopHelper->InsertNextValue(3);
  secondLoopHelper->InsertNextValue(4);

  vtkNew(vtkSphericalConformalMapper, mapper);
  mapper->SetInputData(slicePd);
  mapper->SetVerbose(0);
  mapper->SetInitialTimeStep(0.1);
  mapper->SetTutteEnergyCriterion(1.0e-6);
  mapper->SetHarmonicEnergyCriterion(1.0e-7);
  mapper->SetMaxNumIterations(1e4);
  mapper->SetCGUpdateMethod(vtkSphericalConformalMapper::CG_NONE);
  mapper->SetBoundaryConstraintType(1);
  mapper->SetBoundaryStart(starts);
  mapper->SetFirstLoopPts(firstCorners);
  mapper->SetSecondLoopPts(secondCorners);
  mapper->SetFirstLoopHelper(firstLoopHelper);
  mapper->SetSecondLoopHelper(secondLoopHelper);
  mapper->SetCubeStart(cubeStart);
  mapper->SetObjectXAxis(xvec);
  mapper->SetObjectZAxis(zvec);
  mapper->Update();

  sliceS2Pd->DeepCopy(mapper->GetOutput());

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataToNURBSFilter::MapOpenSliceToS2(vtkPolyData *slicePd,
                                               vtkPolyData *sliceS2Pd,
                                               vtkIntArray *firstCorners,
                                               double xvec[3],
                                               double zvec[3])
{
  int starts[2];
  starts[0] = firstCorners->GetValue(7);
  starts[1] = -1;
  int cubeStart[2]; cubeStart[0] = 0; cubeStart[1] = 4;


  vtkNew(vtkIntArray, firstLoopHelper);
  firstLoopHelper->InsertNextValue(1);
  firstLoopHelper->InsertNextValue(0);
  firstLoopHelper->InsertNextValue(4);
  firstLoopHelper->InsertNextValue(2);
  firstLoopHelper->InsertNextValue(1);
  firstLoopHelper->InsertNextValue(3);
  firstLoopHelper->InsertNextValue(4);
  firstLoopHelper->InsertNextValue(5);

  vtkNew(vtkSphericalConformalMapper, mapper);
  mapper->SetInputData(slicePd);
  mapper->SetVerbose(0);
  mapper->SetInitialTimeStep(0.1);
  mapper->SetTutteEnergyCriterion(1.0e-6);
  mapper->SetHarmonicEnergyCriterion(1.0e-7);
  mapper->SetMaxNumIterations(1e4);
  mapper->SetCGUpdateMethod(vtkSphericalConformalMapper::CG_NONE);
  mapper->SetBoundaryStart(starts);
  mapper->SetObjectXAxis(xvec);
  mapper->SetObjectZAxis(zvec);
  mapper->SetFirstLoopPts(firstCorners);
  mapper->SetFirstLoopHelper(firstLoopHelper);
  mapper->SetBoundaryConstraintType(1);
  mapper->SetCubeStart(cubeStart);
  mapper->Update();

  sliceS2Pd->DeepCopy(mapper->GetOutput());

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataToNURBSFilter::InterpolateMapOntoTarget(vtkPolyData *sourceS2Pd,
                                                            vtkPolyData *targetPd,
                                                            vtkPolyData *targetS2Pd,
                                                            vtkPolyData *mappedPd)
{
  vtkNew(vtkMapInterpolator, interpolator);
  interpolator->SetInputData(0, sourceS2Pd);
  interpolator->SetInputData(1, targetPd);
  interpolator->SetInputData(2, targetS2Pd);
  interpolator->SetVerbose(0);
  interpolator->SetNumSourceSubdivisions(0);
  interpolator->Update();

  mappedPd->DeepCopy(interpolator->GetOutput());

  return 1;
}
