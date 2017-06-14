/*=========================================================================

Program:   VMTK
Module:    $RCSfile: vtkSVPolyBallLine.cxx,v $
Language:  C++
Date:      $Date: 2006/04/06 16:46:43 $
Version:   $Revision: 1.5 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "vtkSVPolyBallLine.h"
#include "vtkPointData.h"
#include "vtkPointLocator.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkSVGlobals.h"
#include "vtkSVGeneralUtils.h"

// ----------------------
// StandardNewMacro
// ----------------------
vtkStandardNewMacro(vtkSVPolyBallLine);

// ----------------------
// Constructor
// ----------------------
vtkSVPolyBallLine::vtkSVPolyBallLine()
{
  this->Input = NULL;
  this->InputCellIds = NULL;
  this->InputCellId = -1;
  this->PolyBallRadiusArrayName = NULL;
  this->LocalCoordinatesArrayName = NULL;
  this->LastPolyBallCellId = -1;
  this->LastPolyBallCellSubId = -1;
  this->LastPolyBallCellPCoord = 0.0;
  this->LastPolyBallCenter[0] = this->LastPolyBallCenter[1] = this->LastPolyBallCenter[2] = 0.0;
  this->LastPolyBallCenterRadius = 0.0;
  this->UseRadiusInformation = 1;
  this->UsePointNormal = 0;
  this->UseRadiusWeighting = 0;
  this->UseLocalCoordinates = 0;
  this->RemoveEndPoints = 0;
  this->PointLocator = vtkPointLocator::New();
}

// ----------------------
// Destructor
// ----------------------
vtkSVPolyBallLine::~vtkSVPolyBallLine()
{
  if (this->Input)
    {
    this->Input->Delete();
    this->Input = NULL;
    }

  if (this->InputCellIds)
    {
    this->InputCellIds->Delete();
    this->InputCellIds = NULL;
    }

  if (this->PolyBallRadiusArrayName)
    {
    delete[] this->PolyBallRadiusArrayName;
    this->PolyBallRadiusArrayName = NULL;
    }
  if (this->PointLocator != NULL)
    {
      this->PointLocator->Delete();
    this->PointLocator = NULL;
    }
}

// ----------------------
// ComplexDot
// ----------------------
double vtkSVPolyBallLine::ComplexDot(double x[4], double y[4])
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2] - x[3]*y[3];
}

// ----------------------
// BuildLocator
// ----------------------
void vtkSVPolyBallLine::BuildLocator()
{
  this->PointLocator->SetDataSet(this->Input);
  this->PointLocator->BuildLocator();
}

// ----------------------
// EvaluateFunction
// ----------------------
double vtkSVPolyBallLine::EvaluateFunction(double x[3])
{
  vtkIdType i, k;
  vtkIdType npts, *pts;
  double polyballFunctionValue, minPolyBallFunctionValue;
  double point0[3], point1[3];
  double radius0, radius1;
  double vector0[4], vector1[4], closestPoint[4];
  double local0[3][3], local1[3][3];
  double localDiffs[3][3], finalLocal[3][3];;
  double t;
  double num, den;
  vtkDataArray *polyballRadiusArray = NULL;
  vtkDataArray *localXArray = NULL;
  vtkDataArray *localYArray = NULL;
  vtkDataArray *localZArray = NULL;

  if (!this->Input)
    {
    vtkErrorMacro(<<"No Input specified!");
    return SV_ERROR;
    }

  if (this->Input->GetNumberOfPoints()==0)
    {
    vtkWarningMacro(<<"Empty Input specified!");
    return SV_ERROR;
    }

  if (this->UseRadiusInformation)
    {
    if (!this->PolyBallRadiusArrayName)
      {
      vtkErrorMacro(<<"No PolyBallRadiusArrayName specified!");
      return SV_ERROR;
      }

    polyballRadiusArray = this->Input->GetPointData()->GetArray(this->PolyBallRadiusArrayName);

    if (polyballRadiusArray==NULL)
      {
      vtkErrorMacro(<<"PolyBallRadiusArray with name specified does not exist!");
      return SV_ERROR;
      }
    }

  if (this->UseLocalCoordinates)
    {
    if (!this->LocalCoordinatesArrayName)
      {
      vtkErrorMacro("Must provide local coordinates name if using local coordinates");
      return SV_ERROR;
      }

    std::string localXName = this->LocalCoordinatesArrayName; localXName = localXName+"X";
    localXArray = this->Input->GetPointData()->GetArray(localXName.c_str());
    std::string localYName = this->LocalCoordinatesArrayName; localYName = localYName+"Y";
    localYArray = this->Input->GetPointData()->GetArray(localYName.c_str());
    std::string localZName = this->LocalCoordinatesArrayName; localZName = localZName+"Z";
    localZArray = this->Input->GetPointData()->GetArray(localZName.c_str());
    }

  if (this->Input->GetLines()==NULL)
    {
    vtkWarningMacro(<<"No lines in Input dataset.");
    return SV_ERROR;
    }

  this->Input->BuildCells();
#if (VTK_MAJOR_VERSION <= 5)
  this->Input->Update();
#endif

  minPolyBallFunctionValue = VTK_SV_LARGE_DOUBLE;

  closestPoint[0] = closestPoint[1] = closestPoint[2] = closestPoint[3] = 0.0;

  this->LastPolyBallCellId = -1;
  this->LastPolyBallCellSubId = -1;
  this->LastPolyBallCellPCoord = 0.0;
  this->LastPolyBallCenter[0] = this->LastPolyBallCenter[1] = this->LastPolyBallCenter[2] = 0.0;
  this->LastPolyBallCenterRadius = 0.0;

  vtkIdList* cellIds = vtkIdList::New();

  if (this->InputCellIds)
    {
    cellIds->DeepCopy(this->InputCellIds);
    }
  else if (this->InputCellId != -1)
    {
    cellIds->InsertNextId(this->InputCellId);
    }
  else
    {
    cellIds->SetNumberOfIds(this->Input->GetNumberOfCells());
    for (k=0; k<this->Input->GetNumberOfCells(); k++)
      {
      cellIds->SetId(k,k);
      }
    }

  vtkNew(vtkIdList, closestPoints);
  this->PointLocator->FindClosestNPoints(10, x, closestPoints);

  for (i=0; i<closestPoints->GetNumberOfIds(); i++)
    {
      int ptId0 = closestPoints->GetId(i);
      this->Input->GetPoint(ptId0, point0);
      vtkNew(vtkIdList, tmpList);
      this->Input->GetPointCells(ptId0, tmpList);
      vtkNew(vtkIdList, pointNeighbors);
      for (int r=0; r<tmpList->GetNumberOfIds(); r++)
      {
        vtkIdType *pts, npts;
        this->Input->GetCellPoints(tmpList->GetId(r), npts, pts);
        for (int f=0; f<npts; f++)
        {
          if (pts[f] == ptId0)
          {
            if (f != 0)
              pointNeighbors->InsertNextId(pts[f-1]);
            if (f != npts-1)
              pointNeighbors->InsertNextId(pts[f+1]);
          }
        }
      }
      for (int r=0; r<pointNeighbors->GetNumberOfIds(); r++)
      {
        int ptId1 = pointNeighbors->GetId(r);
      this->Input->GetPoint(ptId1, point1);

      if (this->UseRadiusInformation)
        {
          radius0 = polyballRadiusArray->GetComponent(ptId0,0);
          radius1 = polyballRadiusArray->GetComponent(ptId1,0);
        }
      else
        {
        radius0 = 0.0;
        radius1 = 0.0;
        }

        int cellId;
        vtkNew(vtkIdList, pointCells0);
        vtkNew(vtkIdList, pointCells1);
        this->Input->GetPointCells(ptId0, pointCells0);
        this->Input->GetPointCells(ptId1, pointCells1);
        pointCells0->IntersectWith(pointCells1);
        cellId = pointCells0->GetId(0);

        vtkNew(vtkIdList, tmpList2);
        this->Input->GetPointCells(ptId1, tmpList2);

      if (this->UseLocalCoordinates)
        {
        localXArray->GetTuple(ptId0, local0[0]);
        localYArray->GetTuple(ptId0, local0[1]);
        localZArray->GetTuple(ptId0, local0[2]);
        localXArray->GetTuple(ptId1, local1[0]);
        localYArray->GetTuple(ptId1, local1[1]);
        localZArray->GetTuple(ptId1, local1[2]);
        }
      else
        {
        for (int j=0; j<3; j++)
          {
          for (int k=0; k<3; k++)
            {
            local0[j][k] = 0.0;
            local1[j][k] = 0.0;
            }
          }
        }

      vector0[0] = point1[0] - point0[0];
      vector0[1] = point1[1] - point0[1];
      vector0[2] = point1[2] - point0[2];
      vector0[3] = radius1 - radius0;
      vector1[0] = x[0] - point0[0];
      vector1[1] = x[1] - point0[1];
      vector1[2] = x[2] - point0[2];
      vector1[3] = 0.0 - radius0;
      for (int j=0; j<3; j++)
        vtkMath::Subtract(local1[j], local0[j], localDiffs[j]);

//       cout<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<point0[0]<<" "<<point0[1]<<" "<<point0[2]<<" "<<point1[0]<<" "<<point1[1]<<" "<<point1[2]<<" "<<endl;

      num = this->ComplexDot(vector0,vector1);
      den = this->ComplexDot(vector0,vector0);

      if (fabs(den)<VTK_SV_DOUBLE_TOL)
        {
        continue;
        }

      t = num / den;

      if (t<VTK_SV_DOUBLE_TOL)
        {
        t = 0.0;
        closestPoint[0] = point0[0];
        closestPoint[1] = point0[1];
        closestPoint[2] = point0[2];
        closestPoint[3] = radius0;
        for (int j=0; j<3; j++)
          {
          for (int k=0; k<3; k++)
            finalLocal[j][k] = local0[j][k];
          }
        }
      else if (1.0-t<VTK_SV_DOUBLE_TOL)
        {
        t = 1.0;
        closestPoint[0] = point1[0];
        closestPoint[1] = point1[1];
        closestPoint[2] = point1[2];
        closestPoint[3] = radius1;
        for (int j=0; j<3; j++)
          {
          for (int k=0; k<3; k++)
            finalLocal[j][k] = local1[j][k];
          }
        }
      else
        {
        closestPoint[0] = point0[0] + t * vector0[0];
        closestPoint[1] = point0[1] + t * vector0[1];
        closestPoint[2] = point0[2] + t * vector0[2];
        closestPoint[3] = radius0 + t * vector0[3];
        for (int j=0; j<3; j++)
          {
          for (int k=0; k<3; k++)
            finalLocal[j][k] = local0[j][k] + t * localDiffs[j][k];
          }
        }

      polyballFunctionValue = (x[0]-closestPoint[0])*(x[0]-closestPoint[0]) + (x[1]-closestPoint[1])*(x[1]-closestPoint[1]) + (x[2]-closestPoint[2])*(x[2]-closestPoint[2]) - closestPoint[3]*closestPoint[3];

      if (this->UsePointNormal)
        {
        double dir0[3];
        vtkMath::Subtract(x, closestPoint, dir0);
        vtkMath::Normalize(dir0);
        double align0 = vtkMath::Dot(this->PointNormal, dir0);

        if (align0 <= 0.5)
          {
          // We found a false positive
          polyballFunctionValue = VTK_SV_LARGE_DOUBLE;
          }
        }

      if (this->UseRadiusWeighting && this->UseRadiusInformation)
        {
        double factor = 0.1;
        double weight = fabs(radius1 - radius0)*fabs(radius1-radius0);
        weight = (weight*factor)/svminimum(radius0, radius1);
        polyballFunctionValue += weight;
        }

       if (polyballFunctionValue<minPolyBallFunctionValue)
        {
          if (tmpList->GetNumberOfIds() > 1)
            this->FindBifurcationCellId(ptId0, tmpList, x, closestPoint, t, cellId);
          else if (tmpList2->GetNumberOfIds() > 1)
            this->FindBifurcationCellId(ptId1, tmpList2, x, closestPoint, 1-t, cellId);
        minPolyBallFunctionValue = polyballFunctionValue;
        this->LastPolyBallCellId = cellId;
        this->LastPolyBallCellSubId = i;
        this->LastPolyBallCellPCoord = t;
        this->LastPolyBallCenter[0] = closestPoint[0];
        this->LastPolyBallCenter[1] = closestPoint[1];
        this->LastPolyBallCenter[2] = closestPoint[2];
        this->LastPolyBallCenterRadius = closestPoint[3];
        for (int j=0; j<3; j++)
          {
          this->LastLocalCoordX[j] = finalLocal[0][j];
          this->LastLocalCoordY[j] = finalLocal[1][j];
          this->LastLocalCoordZ[j] = finalLocal[2][j];
          }
        }
      }
    }

  cellIds->Delete();

  return minPolyBallFunctionValue;
}

// ----------------------
// EvaluateGradient
// ----------------------
void vtkSVPolyBallLine::EvaluateGradient(double x[3], double n[3])
{
  vtkWarningMacro("Poly ball gradient computation not yet implemented!");
  // TODO
}

// ----------------------
// PrintSelf
// ----------------------
void vtkSVPolyBallLine::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

// ----------------------
// GetExcludingCellIds
// ----------------------
int vtkSVPolyBallLine::FindBifurcationCellId(const int ptId, vtkIdList *cellIds,
                                              const double surfacePt[3],
                                              const double closestPt[3],
                                              const double parameter,
                                              int &bestCellId)
{
  vtkNew(vtkIdList, checkIds);
  vtkNew(vtkIdList, cellList);

  for (int i=0; i<cellIds->GetNumberOfIds(); i++)
  {
    vtkIdType npts, *pts;
    this->Input->GetCellPoints(cellIds->GetId(i), npts, pts);
    for (int j=0; j<npts; j++)
    {
      if (pts[j] == ptId)
      {
        if (j == 0)
        {
          checkIds->InsertNextId(pts[j+1]);
          cellList->InsertNextId(cellIds->GetId(i));
        }
        if (j == npts-1)
        {
          checkIds->InsertNextId(pts[j-1]);
          cellList->InsertNextId(cellIds->GetId(i));
        }
      }
    }
  }

  double pt0[3];
  this->Input->GetPoint(ptId, pt0);

  double bigAngle = -1.0;
  int bigPoint = 0;
  for (int i=0; i<checkIds->GetNumberOfIds(); i++)
  {
    double pt1[3];
    this->Input->GetPoint(checkIds->GetId(i), pt1);

    double vec0[3];
    vtkMath::Subtract(pt1, pt0, vec0);
    vtkMath::Normalize(vec0);

    double angleSum = 0.0;
    for (int j=0; j<checkIds->GetNumberOfIds(); j++)
    {
      if (i != j)
      {
        double pt2[3];
        this->Input->GetPoint(checkIds->GetId(j), pt2);

        double vec1[3];
        vtkMath::Subtract(pt2, pt0, vec1);
        vtkMath::Normalize(vec1);

        double crossVec[3];
        vtkMath::Cross(vec0, vec1, crossVec);

        double ang = 180.0*atan2(vtkMath::Norm(crossVec), vtkMath::Dot(vec0, vec1))/SV_PI;
        angleSum += ang;
      }
    }
    if (angleSum > bigAngle)
    {
      bigAngle = angleSum;
      bigPoint = i;
    }
  }

  std::vector<double> angles;
  std::vector<int> cells;
  std::vector<int> points;

  double pt1[3];
  this->Input->GetPoint(checkIds->GetId(bigPoint), pt1);
  angles.push_back(180.0);
  cells.push_back(cellList->GetId(bigPoint));
  points.push_back(checkIds->GetId(bigPoint));

  double vec0[3];
  vtkMath::Subtract(pt1, pt0, vec0);
  vtkMath::Normalize(vec0);

  for (int i=0; i<checkIds->GetNumberOfIds(); i++)
  {
    if (i != bigPoint)
    {
      double pt2[3];
      this->Input->GetPoint(checkIds->GetId(i), pt2);

      double vec1[3];
      vtkMath::Subtract(pt2, pt0, vec1);
      vtkMath::Normalize(vec1);

      double crossVec[3];
      vtkMath::Cross(vec0, vec1, crossVec);

      double ang = 180.0*atan2(vtkMath::Norm(crossVec), vtkMath::Dot(vec0, vec1))/SV_PI;
      angles.push_back(ang);
      cells.push_back(cellList->GetId(i));
      points.push_back(checkIds->GetId(i));
    }
  }

  std::vector<int> isAligning;
  for (int i=0; i<angles.size(); i++)
  {
    if (angles[i] < 135.0)
      isAligning.push_back(0);
    else
      isAligning.push_back(1);
  }

  double surfVec[3];
  vtkMath::Subtract(surfacePt, closestPt, surfVec);
  vtkMath::Normalize(surfVec);

  for (int i=0; i<angles.size(); i++)
  {
    if (!isAligning[i] && cells[i] == bestCellId && parameter < 0.7)
    {
      double pt2[3];
      this->Input->GetPoint(points[i], pt2);

      double vec1[3];
      vtkMath::Subtract(pt2, pt0, vec1);
      vtkMath::Normalize(vec1);

      int inMaxDot = -1.1;
      int inMaxCellId = -1;
      for (int j=0; j<angles.size(); j++)
      {
        if (isAligning[j] && i != j)
        {
          double pt3[3];
          this->Input->GetPoint(points[j], pt3);

          double vec2[3];
          vtkMath::Subtract(pt3, pt0, vec2);
          vtkMath::Normalize(vec2);

          double vec3[3];
          vtkMath::Cross(vec2, vec1, vec3);
          vtkMath::Normalize(vec3);

          double vec4[3];
          vtkMath::Cross(vec1, vec3, vec4);
          vtkMath::Normalize(vec4);

          double compare = vtkMath::Dot(surfVec, vec4);
          if (compare > inMaxDot)
          {
            inMaxDot = compare;
            inMaxCellId = cells[j];
          }
        }
      }
      bestCellId = inMaxCellId;
      break;
    }
  }

  return SV_OK;
}
