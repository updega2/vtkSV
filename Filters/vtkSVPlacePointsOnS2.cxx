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

/** @file vtkSVPlacePointsOnS2.cxx
 *  @brief This implements the vtkSVPlacePointsOnS2 filter as a class
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkSVPlacePointsOnS2.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCenterOfMass.h"
#include "vtkCleanPolyData.h"
#include "vtkFloatArray.h"
#include "vtkGradientFilter.h"
#include "vtkMath.h"
#include "vtkMatrix4x4.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointDataToCellData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmartPointer.h"
#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"
#include "vtkTextureMapToSphere.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTriangle.h"
#include "vtkXMLPolyDataWriter.h"

#include <iostream>
#include <cmath>

//---------------------------------------------------------------------------
//vtkCxxRevisionMacro(vtkSVPlacePointsOnS2, "$Revision: 0.0 $");
vtkStandardNewMacro(vtkSVPlacePointsOnS2);


//---------------------------------------------------------------------------
vtkSVPlacePointsOnS2::vtkSVPlacePointsOnS2()
{
  this->Verbose = 1;

  this->InitialPd = vtkPolyData::New();
  this->FinalPd = vtkPolyData::New();

  this->UseCustomAxisAlign = 0;
  this->SetZAxis(0.0, 0.0, 1.0);
  this->SetXAxis(1.0, 0.0, 0.0);
}

//---------------------------------------------------------------------------
vtkSVPlacePointsOnS2::~vtkSVPlacePointsOnS2()
{
  if (this->InitialPd != NULL)
  {
    this->InitialPd->Delete();
  }
  if (this->FinalPd != NULL)
  {
    this->FinalPd->Delete();
  }
}

//---------------------------------------------------------------------------
void vtkSVPlacePointsOnS2::PrintSelf(ostream& os, vtkIndent indent)
{
}

// Generate Separated Surfaces with Region ID Numbers
//---------------------------------------------------------------------------
int vtkSVPlacePointsOnS2::RequestData(
                                 vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector)
{
  // get the input and output
  vtkPolyData *input1 = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  //Copy the input to operate on
  this->InitialPd->DeepCopy(input1);

  this->FinalPd->SetPoints(this->InitialPd->GetPoints());
  this->FinalPd->SetPolys(this->InitialPd->GetPolys());

  if (this->MoveToOrigin() != 1)
  {
    vtkErrorMacro("Couldn't move to origin\n");
    return SV_ERROR;
  }
  if (this->UseCustomAxisAlign)
  {
    if (this->RotateToCubeCenterAxis() != 1)
    {
      vtkErrorMacro("Couldn't rotate\n");
      return SV_ERROR;
    }
  }
  if (this->ScaleToUnitCube() != 1)
  {
    vtkErrorMacro("Couldn't scale\n");
    return SV_ERROR;
  }

  if (this->DumbMapToSphere() != 1)
  {
    vtkErrorMacro("Point placement failed");
    return SV_ERROR;
  }

  output->DeepCopy(this->FinalPd);
  return SV_OK;
}


//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPlacePointsOnS2::DumbMapToSphere()
{
  if (this->TextureMap() != 1)
  {
    return SV_ERROR;
  }
  if (this->ConvertTextureFieldToPolyData() != 1)
  {
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
int vtkSVPlacePointsOnS2::TextureMap()
{
  vtkNew(vtkTextureMapToSphere, texturer);
  texturer->SetInputData(this->FinalPd);
  texturer->PreventSeamOff();
  texturer->Update();

  this->FinalPd->DeepCopy(texturer->GetOutput());

  return SV_OK;
}
//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPlacePointsOnS2::ConvertTextureFieldToPolyData()
{
  vtkNew(vtkFloatArray, textureCoords);
  textureCoords = vtkFloatArray::SafeDownCast(this->FinalPd->GetPointData()->GetArray("Texture Coordinates"));

  int numPts = this->FinalPd->GetNumberOfPoints();
  for (int i=0; i< numPts; i++)
  {
    double tPt[2];
    textureCoords->GetTuple(i, tPt);

    double pt[3];
    pt[0] = std::sin(M_PI * tPt[1]) * std::cos(2.0 * M_PI * tPt[0]);
    pt[1] = std::sin(M_PI * tPt[1]) * std::sin(2.0 * M_PI * tPt[0]);
    pt[2] = std::cos(M_PI * tPt[1]);

    double realpt[3];
    this->FinalPd->GetPoint(i, realpt);

    this->FinalPd->GetPoints()->SetPoint(i, pt);
  }

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPlacePointsOnS2::RotateToCubeCenterAxis()
{
  double realY[3], realZ[3];
  realY[0] = 0.0; realY[1] = 1.0; realY[2] = 0.0;
  realZ[0] = 0.0; realZ[1] = 0.0; realZ[2] = 1.0;

  double YAxis[3];
  vtkMath::Normalize(this->XAxis);
  vtkMath::Normalize(this->ZAxis);
  vtkMath::Cross(this->ZAxis, this->XAxis, YAxis);
  vtkMath::Normalize(YAxis);
  double inZ[4], outZ[4], rotZ[3];
  for (int i=0; i<3; i++)
  {
    inZ[i] = this->ZAxis[i];
  }
  inZ[3] = 1.0;

  vtkNew(vtkMatrix4x4, rotMatrix0);
  vtkNew(vtkMatrix4x4, rotMatrix1);
  vtkSVGeneralUtils::GetRotationMatrix(YAxis, realY, rotMatrix0);
  vtkSVGeneralUtils::ApplyRotationMatrix(this->FinalPd, rotMatrix0);
  rotMatrix0->MultiplyPoint(inZ, outZ);
  for (int i=0; i<3; i++)
  {
    rotZ[i] = outZ[i];
  }

  vtkSVGeneralUtils::GetRotationMatrix(rotZ, realZ, rotMatrix1);
  vtkSVGeneralUtils::ApplyRotationMatrix(this->FinalPd, rotMatrix1);

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPlacePointsOnS2::MoveToOrigin()
{
  double massCenter[3];
  vtkSVGeneralUtils::ComputeMassCenter(this->FinalPd, massCenter);

  int numPts = this->FinalPd->GetNumberOfPoints();
  for (int i=0; i<numPts; i++)
  {
    double pt[3];
    this->FinalPd->GetPoint(i, pt);

    double movePt[3];
    for (int j=0; j<3; j++)
    {
      movePt[j] = pt[j] - massCenter[j];
    }
    this->FinalPd->GetPoints()->SetPoint(i, movePt);
  }

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPlacePointsOnS2::ScaleToUnitCube()
{
  double bounds[6];
  this->FinalPd->GetBounds(bounds);
  double xScaleFactor = 1.0/(bounds[1]-bounds[0]);
  double yScaleFactor = 1.0/(bounds[3]-bounds[2]);
  double zScaleFactor = 1.0/4*(bounds[5]-bounds[4]);

  vtkNew(vtkTransform, transformer);
  transformer->Scale(xScaleFactor, yScaleFactor, zScaleFactor);
  transformer->Update();

  vtkNew(vtkTransformPolyDataFilter, pdTransformer);
  pdTransformer->SetInputData(this->FinalPd);
  pdTransformer->SetTransform(transformer);
  pdTransformer->Update();

  this->FinalPd->DeepCopy(pdTransformer->GetOutput());
  this->FinalPd->BuildLinks();

  return SV_OK;
}
