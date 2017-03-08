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

/** @file vtkSVPolyDataToNURBSFilter.cxx
 *  @brief This implements the vtkSVPolyDataToNURBSFilter filter as a class
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkSVPolyDataToNURBSFilter.h"

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
#include "vtkSVLoftNURBSSurface.h"
#include "vtkSVMapInterpolator.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPlane.h"
#include "vtkPointData.h"
#include "vtkPointDataToCellData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkSVGlobals.h"
#include "vtkSVPolyDataSliceAndDiceFilter.h"
#include "vtkSmartPointer.h"
#include "vtkSVPlanarMapper.h"
#include "vtkSVPullApartPolyData.h"
#include "vtkSVSuperSquareBoundaryMapper.h"
#include "vtkTextureMapToSphere.h"
#include "vtkThreshold.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTriangle.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLPolyDataWriter.h"

#include <iostream>
#include <cmath>

//---------------------------------------------------------------------------
//vtkCxxRevisionMacro(vtkSVPolyDataToNURBSFilter, "$Revision: 0.0 $");
vtkStandardNewMacro(vtkSVPolyDataToNURBSFilter);


//---------------------------------------------------------------------------
vtkSVPolyDataToNURBSFilter::vtkSVPolyDataToNURBSFilter()
{
  this->AddTextureCoordinates = 1;

  this->InputPd         = vtkPolyData::New();
  this->ParameterizedPd = vtkPolyData::New();
  this->TexturedPd      = vtkPolyData::New();
  this->LoftedPd        = vtkPolyData::New();
  this->Centerlines     = NULL;
  this->CubeS2Pd        = NULL;
  this->OpenCubeS2Pd    = NULL;
  this->SurgeryLines    = vtkPolyData::New();
  this->Polycube        = vtkSVGeneralizedPolycube::New();

  this->BoundaryPointsArrayName = NULL;
  this->GroupIdsArrayName       = NULL;
  this->SegmentIdsArrayName     = NULL;
  this->SliceIdsArrayName       = NULL;
  this->SphereRadiusArrayName   = NULL;
  this->InternalIdsArrayName    = NULL;
  this->DijkstraArrayName       = NULL;
  this->BooleanPathArrayName    = NULL;
}

//---------------------------------------------------------------------------
vtkSVPolyDataToNURBSFilter::~vtkSVPolyDataToNURBSFilter()
{
  if (this->InputPd != NULL)
  {
    this->InputPd->Delete();
    this->InputPd = NULL;
  }
  if (this->CubeS2Pd != NULL)
  {
    this->CubeS2Pd->UnRegister(this);
    this->CubeS2Pd = NULL;
  }
  if (this->OpenCubeS2Pd != NULL)
  {
    this->OpenCubeS2Pd->UnRegister(this);
    this->OpenCubeS2Pd = NULL;
  }
  if (this->ParameterizedPd != NULL)
  {
    this->ParameterizedPd->Delete();
    this->ParameterizedPd = NULL;
  }
  if (this->TexturedPd != NULL)
  {
    this->TexturedPd->Delete();
    this->TexturedPd = NULL;
  }
  if (this->LoftedPd != NULL)
  {
    this->LoftedPd->Delete();
    this->LoftedPd = NULL;
  }
  if (this->Centerlines != NULL)
  {
    this->Centerlines->UnRegister(this);
    this->Centerlines = NULL;
  }
  if (this->Polycube != NULL)
  {
    this->Polycube->Delete();
    this->Polycube = NULL;
  }
  if (this->SurgeryLines != NULL)
  {
    this->SurgeryLines->Delete();
    this->SurgeryLines = NULL;
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
  if (this->BooleanPathArrayName != NULL)
  {
    delete [] this->BooleanPathArrayName;
    this->BooleanPathArrayName = NULL;
  }
}

//---------------------------------------------------------------------------
void vtkSVPolyDataToNURBSFilter::PrintSelf(ostream& os, vtkIndent indent)
{
}

// Generate Separated Surfaces with Region ID Numbers
//---------------------------------------------------------------------------
int vtkSVPolyDataToNURBSFilter::RequestData(
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

  if (this->SliceAndDice() != SV_OK)
  {
    vtkErrorMacro("Error in slicing polydata\n");
    return SV_ERROR;
  }

  if (this->PerformMappings() != SV_OK)
  {
    vtkErrorMacro("Error in perform mappings\n");
    return SV_ERROR;
  }

  output->DeepCopy(this->ParameterizedPd);
  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataToNURBSFilter::ComputeCenterlines()
{
  return SV_ERROR;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataToNURBSFilter::ExtractBranches()
{
  return SV_ERROR;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataToNURBSFilter::SliceAndDice()
{
  vtkNew(vtkSVPolyDataSliceAndDiceFilter, slicer);
  slicer->SetInputData(this->InputPd);
  slicer->SetCenterlines(this->Centerlines);
  slicer->SetSliceLength(1.0);
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
  this->SurgeryLines->DeepCopy(slicer->GetSurgeryLines());

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataToNURBSFilter::PerformMappings()
{
  int numCubes = this->Polycube->GetNumberOfGrids();
  vtkNew(vtkIdFilter, ider);
  ider->SetInputData(this->InputPd);
  ider->SetIdsArrayName(this->InternalIdsArrayName);
  ider->Update();
  this->InputPd->DeepCopy(ider->GetOutput());

  vtkIntArray *segmentIds = vtkIntArray::SafeDownCast(
      this->Polycube->GetCellData()->GetArray(this->SegmentIdsArrayName));
  vtkIntArray *cubeType = vtkIntArray::SafeDownCast(
    this->Polycube->GetCellData()->GetArray("CubeType"));

  vtkNew(vtkAppendPolyData, appender);
  vtkNew(vtkAppendPolyData, inputAppender);
  vtkNew(vtkAppendPolyData, loftAppender);
  for (int i=0; i<numCubes; i++)
  {
    int segmentId = segmentIds->GetValue(i);
    if (cubeType->GetValue(i) == vtkSVGeneralizedPolycube::CUBE_BRANCH)
    {
      this->MapBranch(segmentId, appender, inputAppender, loftAppender);
    }
    else if (cubeType->GetValue(i) == vtkSVGeneralizedPolycube::CUBE_BIFURCATION)
    {
      this->MapBifurcation(segmentId, appender, inputAppender, loftAppender);
    }
  }
  appender->Update();
  inputAppender->Update();
  loftAppender->Update();
  this->ParameterizedPd->DeepCopy(appender->GetOutput());
  this->TexturedPd->DeepCopy(inputAppender->GetOutput());
  this->LoftedPd->DeepCopy(loftAppender->GetOutput());

  this->InputPd->GetCellData()->RemoveArray(this->InternalIdsArrayName);
  this->InputPd->GetPointData()->RemoveArray(this->InternalIdsArrayName);

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataToNURBSFilter::GetSegment(const int segmentId, vtkPolyData *segmentPd,
                                        vtkPolyData *surgeryLinePd)
{
  vtkNew(vtkThreshold, thresholder);
  thresholder->SetInputData(this->InputPd);
  //Set Input Array to 0 port,0 connection,1 for Cell Data, and Regions is the type name
  thresholder->SetInputArrayToProcess(0, 0, 0, 1, this->SegmentIdsArrayName);
  thresholder->ThresholdBetween(segmentId, segmentId);
  thresholder->Update();

  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(thresholder->GetOutput());
  surfacer->Update();

  segmentPd->DeepCopy(surfacer->GetOutput());

  thresholder->SetInputData(this->SurgeryLines);
  thresholder->Update();

  surfacer->SetInputData(thresholder->GetOutput());
  surfacer->Update();

  surgeryLinePd->DeepCopy(surfacer->GetOutput());
  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataToNURBSFilter::GetSlice(const int sliceId, vtkPolyData *segmentPd, vtkPolyData *slicePd)
{
  vtkNew(vtkThreshold, thresholder);
  thresholder->SetInputData(segmentPd);
  //Set Input Array to 0 port,0 connection,1 for Cell Data, and Regions is the type name
  thresholder->SetInputArrayToProcess(0, 0, 0, 1, this->SliceIdsArrayName);
  thresholder->ThresholdBetween(sliceId, sliceId);
  thresholder->Update();

  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(thresholder->GetOutput());
  surfacer->Update();

  slicePd->DeepCopy(surfacer->GetOutput());

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataToNURBSFilter::MapBranch(const int branchId,
                                        vtkAppendPolyData *appender,
                                        vtkAppendPolyData *inputAppender,
                                        vtkAppendPolyData *loftAppender)
{
  vtkNew(vtkPolyData, branchPd);
  vtkNew(vtkPolyData, surgeryLinePd);
  this->GetSegment(branchId, branchPd, surgeryLinePd);

  vtkIntArray *startPtIds =  vtkIntArray::SafeDownCast(
    this->Polycube->GetCellData()->GetArray("CornerPtIds"));
  vtkDoubleArray *sliceRightNormals = vtkDoubleArray::SafeDownCast(
    this->Polycube->GetCellData()->GetArray("RightNormal"));
  vtkDoubleArray *sliceTopNormals = vtkDoubleArray::SafeDownCast(
    this->Polycube->GetCellData()->GetArray("TopNormal"));
  int parentDir = this->Polycube->GetCellData()->GetArray(
    "ParentDirection")->GetTuple1(branchId);
  int childDir = this->Polycube->GetCellData()->GetArray(
    "ChildDirection")->GetTuple1(branchId);
  int cellIndices[8];
  for (int i=0; i<8; i++)
    cellIndices[i] = vtkSVPolyDataSliceAndDiceFilter::LookupIndex(parentDir, childDir, i);

  double minmax[2];
  vtkIntArray *sliceIds = vtkIntArray::SafeDownCast(
    branchPd->GetCellData()->GetArray(this->SliceIdsArrayName));
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
      for (int j=0; j<4; j++)
      {
        firstLoopPts->InsertNextValue(ptIds->LookupValue(startPtIds->GetComponent(branchId, cellIndices[j])));
        secondLoopPts->InsertNextValue(ptIds->LookupValue(startPtIds->GetComponent(branchId, cellIndices[j+4])));
      }
      sliceRightNormals->GetTuple(branchId, xvec);
      sliceTopNormals->GetTuple(branchId, zvec);
      vtkNew(vtkPolyData, sliceS2Pd);
      vtkNew(vtkPolyData, mappedPd);
      fprintf(stdout,"Mapping region %d...\n", branchId);
      fprintf(stdout,"First start vals are: %f %f %f %f\n", startPtIds->GetComponent(branchId, 0),
                                                            startPtIds->GetComponent(branchId, 1),
                                                            startPtIds->GetComponent(branchId, 2),
                                                            startPtIds->GetComponent(branchId, 3));
      fprintf(stdout,"Second start are: %f %f %f %f\n", startPtIds->GetComponent(branchId, 4),
                                                        startPtIds->GetComponent(branchId, 5),
                                                        startPtIds->GetComponent(branchId, 6),
                                                        startPtIds->GetComponent(branchId, 7));

      this->MapSliceToS2(slicePd, surgeryLinePd, sliceS2Pd, firstLoopPts, secondLoopPts, xvec, zvec);
      //this->GetCorrespondingCube(cubeS2Pd, boundary);
      this->InterpolateMapOntoTarget(this->CubeS2Pd, slicePd, sliceS2Pd, mappedPd);
      fprintf(stdout,"Done with mapping...\n");
      appender->AddInputData(mappedPd);
      if (this->AddTextureCoordinates)
      {
        this->UseMapToAddTextureCoordinates(slicePd, sliceS2Pd, 1.0, 1.0);
      }
      inputAppender->AddInputData(slicePd);

      std::stringstream out;
      out << branchId;

      std::string filename = "/Users/adamupdegrove/Desktop/tmp/S2Slice_"+out.str()+".vtp";
      vtkNew(vtkXMLPolyDataWriter, writer);
      writer->SetInputData(sliceS2Pd);
      writer->SetFileName(filename.c_str());
      writer->Write();

      std::string groupfile = "/Users/adamupdegrove/Desktop/tmp/GroupFile_"+out.str();
      this->WriteToGroupsFile(mappedPd, groupfile);
      vtkNew(vtkPolyData, loftedPd);
      //this->LoftNURBSSurface(mappedPd, loftedPd);
      loftAppender->AddInputData(loftedPd);
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
int vtkSVPolyDataToNURBSFilter::MapBifurcation(const int bifurcationId,
                                             vtkAppendPolyData *appender,
                                             vtkAppendPolyData *inputAppender,
                                             vtkAppendPolyData *loftAppender)
{
  vtkNew(vtkPolyData, bifurcationPd);
  vtkNew(vtkPolyData, surgeryLinePd);
  this->GetSegment(bifurcationId, bifurcationPd, surgeryLinePd);

  vtkIntArray *startPtIds =  vtkIntArray::SafeDownCast(
    this->Polycube->GetCellData()->GetArray("CornerPtIds"));
  vtkDoubleArray *sliceRightNormals = vtkDoubleArray::SafeDownCast(
    this->Polycube->GetCellData()->GetArray("RightNormal"));
  vtkDoubleArray *sliceTopNormals = vtkDoubleArray::SafeDownCast(
    this->Polycube->GetCellData()->GetArray("TopNormal"));
  int parentDir = this->Polycube->GetCellData()->GetArray(
    "ParentDirection")->GetTuple1(bifurcationId);
  int childDir = this->Polycube->GetCellData()->GetArray(
    "ChildDirection")->GetTuple1(bifurcationId);
  int cellIndices[8];
  for (int i=0; i<8; i++)
    cellIndices[i] = vtkSVPolyDataSliceAndDiceFilter::LookupIndex(parentDir, childDir, i);

  int numPoints = bifurcationPd->GetNumberOfPoints();
  vtkDataArray *ptIds = bifurcationPd->GetPointData()->GetArray(this->InternalIdsArrayName);

  if (numPoints != 0)
  {
    double xvec[3], zvec[3];
    vtkNew(vtkIntArray, firstLoopPts);
    vtkNew(vtkIntArray, secondLoopPts);
    for (int j=0; j<4; j++)
    {
      firstLoopPts->InsertNextValue(ptIds->LookupValue(startPtIds->GetComponent(bifurcationId, cellIndices[j])));
      secondLoopPts->InsertNextValue(ptIds->LookupValue(startPtIds->GetComponent(bifurcationId, cellIndices[j+4])));
    }
      fprintf(stdout,"First Loop Points: %f %f %f %f\n", startPtIds->GetComponent(bifurcationId, 0),
                                                         startPtIds->GetComponent(bifurcationId, 1),
                                                         startPtIds->GetComponent(bifurcationId, 2),
                                                         startPtIds->GetComponent(bifurcationId, 3));
      fprintf(stdout,"Second Loop Points: %f %f %f %f\n", startPtIds->GetComponent(bifurcationId, 4),
                                                          startPtIds->GetComponent(bifurcationId, 5),
                                                          startPtIds->GetComponent(bifurcationId, 6),
                                                          startPtIds->GetComponent(bifurcationId, 7));
    sliceRightNormals->GetTuple(bifurcationId, xvec);
    sliceTopNormals->GetTuple(bifurcationId, zvec);
    vtkNew(vtkPolyData, sliceS2Pd);
    vtkNew(vtkPolyData, mappedPd);
    fprintf(stdout,"Mapping region %d...\n", bifurcationId);

    this->MapOpenSliceToS2(bifurcationPd, sliceS2Pd, firstLoopPts, secondLoopPts, xvec, zvec);
    //this->GetCorrespondingCube(cubeS2Pd, boundary);
    this->InterpolateMapOntoTarget(this->OpenCubeS2Pd, bifurcationPd, sliceS2Pd, mappedPd);
    fprintf(stdout,"Done with mapping...\n");
    appender->AddInputData(mappedPd);

    if (this->AddTextureCoordinates)
    {
      this->UseMapToAddTextureCoordinates(bifurcationPd, sliceS2Pd, 1.0, 1.0);
    }
    inputAppender->AddInputData(bifurcationPd);

    std::stringstream out;
    out << bifurcationId;

    std::string filename = "/Users/adamupdegrove/Desktop/tmp/S2Slice_"+out.str()+".vtp";
    vtkNew(vtkXMLPolyDataWriter, writer);
    writer->SetInputData(sliceS2Pd);
    writer->SetFileName(filename.c_str());
    writer->Write();

    std::string groupfile = "/Users/adamupdegrove/Desktop/tmp/GroupFile_"+out.str();
    this->WriteToGroupsFile(mappedPd, groupfile);
    vtkNew(vtkPolyData, loftedPd);
    //this->LoftNURBSSurface(mappedPd, loftedPd);
    loftAppender->AddInputData(loftedPd);

  }

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataToNURBSFilter::MapSliceToS2(vtkPolyData *slicePd,
                                       vtkPolyData *surgeryLinePd,
                                       vtkPolyData *sliceS2Pd,
                                       vtkIntArray *firstCorners,
                                       vtkIntArray *secondCorners,
                                       double xvec[3],
                                       double zvec[3])
{
  fprintf(stdout,"Four diff values are: %d %d %d %d\n", firstCorners->GetValue(0),
                                                        firstCorners->GetValue(1),
                                                        firstCorners->GetValue(2),
                                                        firstCorners->GetValue(3));
  fprintf(stdout,"Other values are: %d %d %d %d\n", secondCorners->GetValue(0),
                                                    secondCorners->GetValue(1),
                                                    secondCorners->GetValue(2),
                                                    secondCorners->GetValue(3));

  vtkNew(vtkIntArray, ripIds);
  vtkIntArray *seamIds = vtkIntArray::SafeDownCast(surgeryLinePd->GetPointData()->GetArray(this->InternalIdsArrayName));
  vtkDataArray *pointIds = slicePd->GetPointData()->GetArray(this->InternalIdsArrayName);
  ripIds->SetNumberOfComponents(1);
  ripIds->SetNumberOfTuples(seamIds->GetNumberOfTuples());
  for (int i=0; i<seamIds->GetNumberOfTuples(); i++)
  {
    ripIds->SetValue(i, pointIds->LookupValue(seamIds->GetValue(i)));
  }

  //fprintf(stdout,"XVec: %.4f %.4f %.4f\n", xvec[0], xvec[1], xvec[2]);
  //fprintf(stdout,"ZVec: %.4f %.4f %.4f\n", zvec[0], zvec[1], zvec[2]);
  vtkNew(vtkSVPullApartPolyData, ripper);
  ripper->SetInputData(slicePd);
  ripper->SetStartPtId(firstCorners->GetValue(0));
  ripper->SetObjectXAxis(xvec);
  ripper->SetObjectZAxis(zvec);
  ripper->SetCutPointsArrayName(this->BooleanPathArrayName);
  ripper->SetSeamPointIds(ripIds);
  ripper->Update();
  std::string filename = "/Users/adamupdegrove/Desktop/tmp/Totsally.vtp";
  vtkNew(vtkXMLPolyDataWriter, writer);
  writer->SetInputData(ripper->GetOutput());
  writer->SetFileName(filename.c_str());
  writer->Write();
  slicePd->DeepCopy(ripper->GetOutput());

  vtkNew(vtkIdList, replacedPoints);
  vtkNew(vtkIdList, newPoints);
  replacedPoints->DeepCopy(ripper->GetReplacePointList());
  newPoints->DeepCopy(ripper->GetNewPointList());

  int loc0 = replacedPoints->IsId(firstCorners->GetValue(0));
  int loc1 = replacedPoints->IsId(secondCorners->GetValue(0));
  int new0 = newPoints->GetId(loc0);
  int new1 = newPoints->GetId(loc1);
  fprintf(stdout,"The final order of bpoints: %d %d %d %d %d %d %d %d %d %d\n",
    firstCorners->GetValue(0),
    firstCorners->GetValue(1),
    firstCorners->GetValue(2),
    firstCorners->GetValue(3),
    new0,
    new1,
    secondCorners->GetValue(3),
    secondCorners->GetValue(2),
    secondCorners->GetValue(1),
    secondCorners->GetValue(0));

  vtkNew(vtkIntArray, boundaryCorners);
  boundaryCorners->SetNumberOfValues(10);
  boundaryCorners->SetValue(0, firstCorners->GetValue(0));
  boundaryCorners->SetValue(1, firstCorners->GetValue(1));
  boundaryCorners->SetValue(2, firstCorners->GetValue(2));
  boundaryCorners->SetValue(3, firstCorners->GetValue(3));
  boundaryCorners->SetValue(4, new0);
  boundaryCorners->SetValue(5, new1);
  boundaryCorners->SetValue(6, secondCorners->GetValue(3));
  boundaryCorners->SetValue(7, secondCorners->GetValue(2));
  boundaryCorners->SetValue(8, secondCorners->GetValue(1));
  boundaryCorners->SetValue(9, secondCorners->GetValue(0));
  double boundaryLengths[4];
  boundaryLengths[0] = 4.0; boundaryLengths[1] = 1.0; boundaryLengths[2] = 4.0; boundaryLengths[3] = 1.0;
  int boundaryDivisions[4];
  boundaryDivisions[0] = 3; boundaryDivisions[1] = 0; boundaryDivisions[2] = 3; boundaryDivisions[3] = 0;

  vtkNew(vtkSVSuperSquareBoundaryMapper, boundaryMapper);
  boundaryMapper->SetBoundaryIds(boundaryCorners);
  boundaryMapper->SetSuperBoundaryDivisions(boundaryDivisions);
  boundaryMapper->SetSuperBoundaryLengths(boundaryLengths);
  boundaryMapper->SetObjectXAxis(xvec);
  boundaryMapper->SetObjectZAxis(zvec);

  vtkNew(vtkSVPlanarMapper, mapper);
  mapper->SetInputData(slicePd);
  //mapper->SetTutteEnergyCriterion(1.0e-6);
  //mapper->SetHarmonicEnergyCriterion(1.0e-7);
  //mapper->SetMaxNumIterations(1e4);
  mapper->SetBoundaryMapper(boundaryMapper);
  mapper->Update();

  sliceS2Pd->DeepCopy(mapper->GetOutput());

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataToNURBSFilter::MapOpenSliceToS2(vtkPolyData *slicePd,
                                               vtkPolyData *sliceS2Pd,
                                               vtkIntArray *firstCorners,
                                               vtkIntArray *secondCorners,
                                               double xvec[3],
                                               double zvec[3])
{
  fprintf(stdout,"Four diff values are: %d %d %d %d\n", firstCorners->GetValue(0),
                                                        firstCorners->GetValue(1),
                                                        firstCorners->GetValue(2),
                                                        firstCorners->GetValue(3));
  fprintf(stdout,"Other values are: %d %d %d %d\n", secondCorners->GetValue(0),
                                                    secondCorners->GetValue(1),
                                                    secondCorners->GetValue(2),
                                                    secondCorners->GetValue(3));

  vtkNew(vtkIntArray, boundaryCorners);
  boundaryCorners->SetNumberOfValues(8);
  boundaryCorners->SetValue(0, firstCorners->GetValue(0));
  boundaryCorners->SetValue(1, firstCorners->GetValue(3));
  boundaryCorners->SetValue(2, firstCorners->GetValue(2));
  boundaryCorners->SetValue(3, firstCorners->GetValue(1));
  boundaryCorners->SetValue(4, secondCorners->GetValue(1));
  boundaryCorners->SetValue(5, secondCorners->GetValue(2));
  boundaryCorners->SetValue(6, secondCorners->GetValue(3));
  boundaryCorners->SetValue(7, secondCorners->GetValue(0));
  fprintf(stdout,"The final order of bpoints: %d %d %d %d %d %d %d %d\n",
    firstCorners->GetValue(0),
    firstCorners->GetValue(3),
    firstCorners->GetValue(2),
    firstCorners->GetValue(1),
    secondCorners->GetValue(1),
    secondCorners->GetValue(2),
    secondCorners->GetValue(3),
    secondCorners->GetValue(0));
  double boundaryLengths[4];
  boundaryLengths[0] = 3.0; boundaryLengths[1] = 1.0; boundaryLengths[2] = 3.0; boundaryLengths[3] = 1.0;
  int boundaryDivisions[4];
  boundaryDivisions[0] = 2; boundaryDivisions[1] = 0; boundaryDivisions[2] = 2; boundaryDivisions[3] = 0;

  vtkNew(vtkSVSuperSquareBoundaryMapper, boundaryMapper);
  boundaryMapper->SetBoundaryIds(boundaryCorners);
  boundaryMapper->SetSuperBoundaryDivisions(boundaryDivisions);
  boundaryMapper->SetSuperBoundaryLengths(boundaryLengths);
  fprintf(stdout,"XVec: %.4f %.4f %.4f\n", xvec[0], xvec[1], xvec[2]);
  fprintf(stdout,"ZVec: %.4f %.4f %.4f\n", zvec[0], zvec[1], zvec[2]);
  boundaryMapper->SetObjectXAxis(xvec);
  boundaryMapper->SetObjectZAxis(zvec);

  vtkNew(vtkSVPlanarMapper, mapper);
  mapper->SetInputData(slicePd);
  //mapper->SetTutteEnergyCriterion(1.0e-6);
  //mapper->SetHarmonicEnergyCriterion(1.0e-7);
  //mapper->SetMaxNumIterations(1e4);
  mapper->SetBoundaryMapper(boundaryMapper);
  mapper->Update();

  sliceS2Pd->DeepCopy(mapper->GetOutput());

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataToNURBSFilter::InterpolateMapOntoTarget(vtkPolyData *sourceS2Pd,
                                                            vtkPolyData *targetPd,
                                                            vtkPolyData *targetS2Pd,
                                                            vtkPolyData *mappedPd)
{
  vtkNew(vtkSVMapInterpolator, interpolator);
  interpolator->SetInputData(0, sourceS2Pd);
  interpolator->SetInputData(1, targetPd);
  interpolator->SetInputData(2, targetS2Pd);
  interpolator->SetNumSourceSubdivisions(0);
  interpolator->Update();

  mappedPd->DeepCopy(interpolator->GetOutput());

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataToNURBSFilter::UseMapToAddTextureCoordinates(vtkPolyData *pd,
                                                            vtkPolyData *mappedPd,
                                                            const double xSize,
                                                            const double ySize)
{
  int numPoints = mappedPd->GetNumberOfPoints();
  vtkNew(vtkFloatArray, textureCoordinates);
  textureCoordinates->SetNumberOfComponents(3);
  textureCoordinates->SetNumberOfTuples(numPoints);
  textureCoordinates->SetName("TextureCoordinates");
  for (int i=0; i<numPoints; i++)
  {
    double pt[3], tc[3];
    mappedPd->GetPoint(i, pt);
    tc[0] = pt[0]/xSize;
    tc[1] = pt[1]/ySize;
    tc[2] = 0.0;
    textureCoordinates->SetTuple(i, tc);
  }

  pd->GetPointData()->SetTCoords(textureCoordinates);

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataToNURBSFilter::LoftNURBSSurface(vtkPolyData *pd, vtkPolyData *loftedPd)
{
  vtkFloatArray *tCoords = vtkFloatArray::SafeDownCast(pd->GetPointData()->GetArray("TextureCoordinates"));
  int numPts = pd->GetNumberOfPoints();

  double xSpacing, ySpacing;
  vtkSVPolyDataToNURBSFilter::GetSpacingOfTCoords(pd, xSpacing, ySpacing);

  vtkNew(vtkSVLoftNURBSSurface, lofter);

  int xNum = 1.0/xSpacing + 2;
  int yNum = 1.0/ySpacing + 2;
  fprintf(stdout,"XNum: %d\n", xNum);
  fprintf(stdout,"YNum: %d\n", yNum);
  for (int i=0; i<yNum; i++)
  {
    vtkNew(vtkPoints, newPoints);
    newPoints->SetNumberOfPoints(xNum);
    for (int j=0; j<xNum; j++)
    {
      double pt[3];
      pd->GetPoint(i*xNum+j, pt);
      newPoints->SetPoint(j,pt);

    }
    vtkNew(vtkPolyData, newPoly);
    newPoly->SetPoints(newPoints);
    lofter->AddInputData(newPoly);
  }
  lofter->SetUDegree(2);
  lofter->SetVDegree(2);
  lofter->SetPolyDataUSpacing(0.1);
  lofter->SetPolyDataVSpacing(0.1);
  //lofter->SetUKnotSpanType("average");
  //lofter->SetVKnotSpanType("average");
  lofter->SetUKnotSpanType("derivative");
  lofter->SetVKnotSpanType("derivative");
  lofter->SetUParametricSpanType("chord");
  lofter->SetVParametricSpanType("chord");
  lofter->Update();

  loftedPd->DeepCopy(lofter->GetOutput());

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataToNURBSFilter::WriteToGroupsFile(vtkPolyData *pd, std::string fileName)
{
  vtkFloatArray *tCoords = vtkFloatArray::SafeDownCast(pd->GetPointData()->GetArray("TextureCoordinates"));
  int numPts = pd->GetNumberOfPoints();

  double xSpacing, ySpacing;
  vtkSVPolyDataToNURBSFilter::GetSpacingOfTCoords(pd, xSpacing, ySpacing);

  FILE *pFile;
  pFile = fopen(fileName.c_str(), "w");
  if (pFile == NULL)
  {
    fprintf(stderr,"Error opening file\n");
    return SV_ERROR;
  }

  int xNum = 1.0/xSpacing + 2;
  int yNum = 1.0/ySpacing + 2;
  fprintf(stdout,"XNum: %d\n", xNum);
  fprintf(stdout,"YNum: %d\n", yNum);
  for (int i=0; i<yNum; i++)
  {
    fprintf(pFile, "/group/test/%d\n", i);
    fprintf(pFile, "%d\n", i);
    fprintf(pFile, "center_x 0.0\n");
    for (int j=0; j<xNum; j++)
    {
      //int ptId = newPointOrder->GetValue(i*xNum + j);
      //fprintf(stdout,"New point is: %d\n", ptId);
      double pt[3];
      pd->GetPoint(i*xNum+j, pt);
      fprintf(pFile, "%.6f %.6f %.6f\n", pt[0], pt[1], pt[2]);
    }
    fprintf(pFile, "\n");
  }

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataToNURBSFilter::GetSpacingOfTCoords(vtkPolyData *pd, double &xSpacing, double &ySpacing)
{
  vtkFloatArray *tCoords = vtkFloatArray::SafeDownCast(pd->GetPointData()->GetArray("TextureCoordinates"));
  int numPts = pd->GetNumberOfPoints();

  double xMin = 1.0e9;
  double yMin = 1.0e9;
  for (int i=0; i<numPts; i++)
  {
    double ptVal[2];
    tCoords->GetTuple(i, ptVal);
    if (ptVal[0] < xMin && ptVal[0] > 1e-8)
    {
      xMin = ptVal[0];
    }
    if (ptVal[1] < yMin && ptVal[1] > 1e-8)
    {
      yMin = ptVal[1];
    }
  }

  fprintf(stdout,"xMin: %.4f\n", xMin);
  fprintf(stdout,"yMin: %.4f\n", yMin);
  xSpacing = xMin;
  ySpacing = yMin;
  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataToNURBSFilter::GetNewPointOrder(vtkPolyData *pd, double xSpacing, double ySpacing,
                                            vtkIntArray *newPointOrder)
{
  vtkFloatArray *tCoords = vtkFloatArray::SafeDownCast(pd->GetPointData()->GetArray("TextureCoordinates"));
  int numPts = pd->GetNumberOfPoints();

  int xNum = 1.0/xSpacing + 1;
  for (int i=0; i<numPts; i++)
  {
    double ptVal[2];
    tCoords->GetTuple(i, ptVal);

    double xLoc = ptVal[0]/xSpacing;
    int loc = ceil(xLoc + xNum * ptVal[1]/ySpacing);

    double pt[3];
    pd->GetPoint(i, pt);
    fprintf(stdout,"Pt Val: %.4f %.4f, Loc: %d end\n", ptVal[0], ptVal[1], loc);
    fprintf(stdout,"Point: %.6f %.6f %.6f\n", pt[0], pt[1], pt[2]);
    newPointOrder->InsertValue(loc, i);
  }
  return SV_OK;
}
