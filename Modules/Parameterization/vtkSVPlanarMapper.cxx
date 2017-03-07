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

/** @file vtkSVPlanarMapper.cxx
 *  @brief This implements the vtkSVPlanarMapper filter as a class
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkSVPlanarMapper.h"

#include "vtkSVBoundaryMapper.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkConnectivityFilter.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkFeatureEdges.h"
#include "vtkIdFilter.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"
#include "vtkTransform.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLPolyDataWriter.h"

#include <iostream>
#include <sstream>
#include <cmath>

//---------------------------------------------------------------------------
//vtkCxxRevisionMacro(vtkSVPlanarMapper, "$Revision: 0.0 $");
vtkStandardNewMacro(vtkSVPlanarMapper);


//---------------------------------------------------------------------------
vtkSVPlanarMapper::vtkSVPlanarMapper()
{
  this->RemoveInternalIds = 1;

  this->InitialPd     = vtkPolyData::New();
  this->WorkPd        = vtkPolyData::New();
  this->PlanarPd      = vtkPolyData::New();
  this->EdgeTable     = vtkEdgeTable::New();
  this->EdgeWeights   = vtkFloatArray::New();
  this->EdgeNeighbors = vtkIntArray::New();
  this->IsBoundary    = vtkIntArray::New();
  this->Boundaries    = vtkPolyData::New();
  this->BoundaryLoop  = vtkPolyData::New();

  this->BoundaryMapper = NULL;

  this->InternalIdsArrayName = NULL;

  this->Lambda = 0.5;
  this->Mu     = 0.5;

  this->ATutte = NULL;
  this->AHarm  = NULL;
}

//---------------------------------------------------------------------------
vtkSVPlanarMapper::~vtkSVPlanarMapper()
{
  if (this->InitialPd != NULL)
  {
    InitialPd->Delete();
  }
  if (this->WorkPd != NULL)
  {
    WorkPd->Delete();
  }
  if (this->PlanarPd != NULL)
  {
    PlanarPd->Delete();
  }
  if (this->EdgeTable != NULL)
  {
    this->EdgeTable->Delete();
  }
  if (this->EdgeWeights != NULL)
  {
    this->EdgeWeights->Delete();
  }
  if (this->EdgeNeighbors != NULL)
  {
    this->EdgeNeighbors->Delete();
  }
  if (this->IsBoundary != NULL)
  {
    this->IsBoundary->Delete();
  }
  if (this->Boundaries != NULL)
  {
    this->Boundaries->Delete();
  }
  if (this->BoundaryLoop != NULL)
  {
    this->BoundaryLoop->Delete();
  }

  if (this->InternalIdsArrayName)
  {
    delete [] this->InternalIdsArrayName;
    this->InternalIdsArrayName = NULL;
  }

  if (this->ATutte != NULL)
  {
    delete this->ATutte;
  }
  if (this->AHarm != NULL)
  {
    delete this->AHarm;
  }
}

//---------------------------------------------------------------------------
void vtkSVPlanarMapper::PrintSelf(ostream& os, vtkIndent indent)
{
}

// Generate Separated Surfaces with Region ID Numbers
//---------------------------------------------------------------------------
int vtkSVPlanarMapper::RequestData(
                                 vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector)
{
  // get the input and output
  vtkPolyData *input = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  //Copy the input to operate on
  this->InitialPd->DeepCopy(input);

  //Copy information to the working polydata
  this->WorkPd->DeepCopy(this->InitialPd);
  this->PlanarPd->DeepCopy(this->InitialPd);

  if (this->PrepFilter() != 1)
  {
    vtkErrorMacro("Error when mapping");
    output->DeepCopy(this->InitialPd);
    return 0;
  }

  if (this->RunFilter() != 1)
  {
    vtkErrorMacro("Error when mapping");
    output->DeepCopy(this->InitialPd);
    return 0;
  }

  if (this->RemoveInternalIds)
  {
    this->WorkPd->GetPointData()->RemoveArray(this->InternalIdsArrayName);
    this->WorkPd->GetCellData()->RemoveArray(this->InternalIdsArrayName);
  }
  output->DeepCopy(this->PlanarPd);
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPlanarMapper::CheckSurface(vtkPolyData *pd)
{
  pd->BuildLinks();

  int numPts = pd->GetNumberOfPoints();
  int numPolys = pd->GetNumberOfCells();

  for (int i=0; i<numPolys; i++)
  {
    vtkIdType npts, *pts;
    pd->GetCellPoints(i, npts, pts);
    if (npts != 3)
    {
      //vtkErrorMacro("Surface contains elements that aren't triangles");
      return 0;
    }
    for (int j=0; j<npts; j++)
    {
      vtkIdType p0, p1;
      p0 = pts[j];
      p1 = pts[(j+1)%npts];

      vtkNew(vtkIdList, edgeNeighbor);
      pd->GetCellEdgeNeighbors(i, p0, p1, edgeNeighbor);

      if (edgeNeighbor->GetNumberOfIds() > 1)
      {
        //vtkErrorMacro("Surface contains triangles with multiple neighbors, not manifold");
        return 0;
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
int vtkSVPlanarMapper::PrepFilter()
{
  vtkIdType numPolys = this->InitialPd->GetNumberOfPolys();
  vtkIdType numPoints = this->InitialPd->GetNumberOfPoints();
  //Check the input to make sure it is there
  if (numPolys < 1)
  {
    vtkDebugMacro("No input!");
    return 0;
  }

  //Check the input to make sure it is manifold and a triangulated surface
  if (this->CheckSurface(this->InitialPd) != 1)
  {
    vtkErrorMacro("Error when checking input surface");
    return 0;
  }

  // Check if internal id array name is given
  if (!this->InternalIdsArrayName)
  {
    vtkDebugMacro("Internal Ids Array Name not given, setting to InternalIds");
    this->InternalIdsArrayName = new char[strlen("InternalIds") + 1];
    strcpy(this->InternalIdsArrayName, "InternalIds");
  }
  // Check if array internal ids is already on pd
  if (this->CheckArrayExists(this->WorkPd, 0, this->InternalIdsArrayName))
  {
    this->RemoveInternalIds = 0;
  }
  else
  {
    vtkNew(vtkIdFilter, ider);
    ider->SetInputData(this->WorkPd);
    ider->SetIdsArrayName(this->InternalIdsArrayName);
    ider->Update();
    this->WorkPd->DeepCopy(ider->GetOutput());
  }

  //Create the edge table for the input surface
  this->WorkPd->BuildLinks();
  if (!this->CreateEdgeTable(this->WorkPd, this->EdgeTable, this->EdgeWeights,
                             this->EdgeNeighbors, this->IsBoundary))
  {
    vtkErrorMacro("Could not create edge table");
    return 0;
  }

  // Set the size of the matrices
  this->ATutte = new SparseMatrix(numPoints, numPoints);
  this->AHarm = new SparseMatrix(numPoints, numPoints);
  this->Xu.resize(numPoints, 0.0);
  this->Xv.resize(numPoints, 0.0);
  this->Bu.resize(numPoints, 0.0);
  this->Bv.resize(numPoints, 0.0);

  return 1;
}
//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPlanarMapper::CreateEdgeTable(vtkPolyData *pd,
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
    //Insert edge into table
    vtkIdType npts, *pts;
    pd->GetCellPoints(i, npts, pts);
    for (int j=0; j<npts; j++)
    {
      vtkIdType p0 = pts[j];
      vtkIdType p1 = pts[(j+1)%npts];
      vtkNew(vtkIdList, neighborCellIds);
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
        vtkSVPlanarMapper::ComputeEdgeWeight(pd, i, neighborCellId,
                                           p0, p1, weight);
        vtkIdType edgeId = edgeTable->InsertEdge(p0, p1);
        edgeWeights->InsertValue(edgeId, weight);
        edgeNeighbors->InsertValue(edgeId, neighborCellId);
        if (weight < 0)
        {
          //vtkWarningMacro("Negative weight on edge between cells " << i <<
          //  " and "<< neighborCellId << ": " << weight);
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
int vtkSVPlanarMapper::ComputeEdgeWeight(vtkPolyData *pd,
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
          vtkSVPlanarMapper::GetEdgeCotangentAngle(v0, v1, v2, angle);

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
int vtkSVPlanarMapper::GetEdgeCotangentAngle(double pt0[3], double pt1[3],
                                           double pt2[3], double &angle)
{
  double area = 0.0;
  vtkSVPlanarMapper::ComputeArea(pt0, pt1, pt2, area);
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
int vtkSVPlanarMapper::ComputeArea(double pt0[3], double pt1[3],
                                 double pt2[3], double &area)
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
int vtkSVPlanarMapper::RunFilter()
{
  if (this->SetBoundaries() != 1)
  {
    vtkErrorMacro("Error in mapping");
    return 0;
  }

  if (this->SetInternalNodes() != 1)
  {
    vtkErrorMacro("Error setting internal nodes");
    return 0;
  }

  if (this->SolveSystem() != 1)
  {
    vtkErrorMacro("Error solving system");
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
int vtkSVPlanarMapper::SetBoundaries()
{
  this->BoundaryMapper->SetInputData(this->WorkPd);
  this->BoundaryMapper->SetEdgeTable(this->EdgeTable);
  this->BoundaryMapper->SetInternalIdsArrayName(this->InternalIdsArrayName);
  this->BoundaryMapper->Update();

  vtkNew(vtkPolyData, boundaryPd);
  boundaryPd->DeepCopy(this->BoundaryMapper->GetOutput());
  // Check if array internal ids is already on pd
  if (this->CheckArrayExists(boundaryPd, 0, this->InternalIdsArrayName) == 0)
  {
    vtkErrorMacro("No internal ids array name on boundary pd");
    return 0;
  }
  vtkDataArray *originalIds = boundaryPd->GetPointData()->GetArray(this->InternalIdsArrayName);
  int numBoundPts = boundaryPd->GetNumberOfPoints();
  for (int i=0; i<numBoundPts; i++)
  {
    int id = originalIds->GetTuple1(i);
    double pt[3];
    boundaryPd->GetPoint(i, pt);
    // Set diagonal to be 1 for boundary points
    this->AHarm->set_element(id, id, 1.0);
    this->ATutte->set_element(id, id, 1.0);
    // Set right hand side to be point on boundary
    this->Bu[id] = pt[0];
    this->Bv[id] = pt[1];
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPlanarMapper::SetCircleBoundary()
{
//  vtkDataArray *pointIds = this->BoundaryLoop->GetPointData()->GetArray(this->InternalIdsArrayName);
//  int numLines = this->BoundaryLoop->GetNumberOfLines();
//
//  double currCoords[3];
//  for (int i=0; i<3; i++)
//  {
//    currCoords[i] = 0.0;
//  }
//
//  double unitLength = M_PI / 2.0;
//  int currCell = 0;
//  int checkPt = -1;
//  for (int i=0; i<4; i++)
//  {
//    double currLength = 0.0;
//    int lastPt = pointIds->LookupValue(this->BoundaryCorners[i]);
//    int dir = 0; //TODO: Figure out dir
//    vtkIdType npts, *pts;
//    while (checkPt != lastPt)
//    {
//      this->BoundaryLoop->GetCellPoints(currCell, npts, pts);
//      double pt0[3], pt1[3];
//      this->BoundaryLoop->GetPoint(pts[0], pt0);
//      this->BoundaryLoop->GetPoint(pts[1], pt1);
//
//      checkPt = pts[1];
//
//      double dist = std::sqrt(std::pow(pt0[0]-pt1[0], 2.0) +
//                              std::pow(pt0[1]-pt1[1], 2.0) +
//                              std::pow(pt0[2]-pt1[2], 2.0));
//      currLength += dist;
//
//      double boundaryVal[3];
//
//      double angle = currLength/this->BoundaryLengths[i] * unitLength + unitLength*i;// + unitLength/2.0;
//      //boundaryVal[0] = (2.0 * radius * std::cos(angle))/(1.0 + std::pow(radius, 2.0));
//      //boundaryVal[1] = (2.0 * radius * std::sin(angle))/(1.0 + std::pow(radius, 2.0));
//      //boundaryVal[2] = (-1.0 + std::pow(radius, 2.0))/(1.0 + std::pow(radius, 2.0));
//
//      int id = pointIds->GetTuple1(pts[1]);
//
//      //this->HarmonicMap[0]->GetPoints()->SetPoint(id, boundaryVal);
//
//      currCell++;
//    }
//  }

  return 1;
}

// ----------------------
// CheckArrayExists
// ----------------------
/**
 * @brief Function to check is array with name exists in cell or point data
 * @param object this is the object to check if the array exists
 * @param datatype this is point or cell. point =0,cell=1
 * @param arrayname this is the name of the array to check
 * @reutrn this returns 1 if the array exists and zero if it doesn't
 * or the function does not return properly.
 */

int vtkSVPlanarMapper::CheckArrayExists(vtkPolyData *object,int datatype,std::string arrayname )
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
int vtkSVPlanarMapper::SetInternalNodes()
{
  int numPoints = this->WorkPd->GetNumberOfPoints();

  for (int i=0; i<numPoints; i++)
  {
    if (this->IsBoundary->GetValue(i) == 0)
    {
      double tot_weight = 0.0;
      double tot_tutte_weight = 0.0;
      vtkNew(vtkIdList, pointNeighbors);
      vtkSVPlanarMapper::GetPointNeighbors(i, this->WorkPd, pointNeighbors);
      double weight_tot;
      for (int j=0; j<pointNeighbors->GetNumberOfIds(); j++)
      {
        int p1 = pointNeighbors->GetId(j);

        // Get edge info
        vtkIdType edgeId = this->EdgeTable->IsEdge(i, p1);
        int edgeNeighbor = this->EdgeNeighbors->GetValue(edgeId);
        double weight    = this->EdgeWeights->GetValue(edgeId);

        if (edgeNeighbor == -1)
        {
          continue;
        }

        this->AHarm->set_element(i,p1, weight);
        this->ATutte->set_element(i, p1, 1.0);
        tot_weight -= weight;
        tot_tutte_weight -= 1.0;
      }
      double pt[3];
      this->WorkPd->GetPoint(i, pt);
      this->AHarm->set_element(i, i, tot_weight);
      this->ATutte->set_element(i, i, tot_tutte_weight);
      this->Xu[i] = pt[0];
      this->Xv[i] = pt[1];
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
int vtkSVPlanarMapper::GetPointNeighbors(vtkIdType p0,
                                       vtkPolyData *pd,
						                           vtkIdList *pointNeighbors)
{
  //Assuming that pointNeighbors is set with no neighbors already
  vtkNew(vtkIdList, cellIdList);
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
int vtkSVPlanarMapper::SolveSystem()
{
  //this->PrintMatrix(this->A);
  //std::vector<std::vector<double> > Ainv;
  //if(!vtkSVPlanarMapper::InvertSystem(this->A, Ainv))
  //{
  //  vtkErrorMacro("Could not invert system");
  //  return 0;
  //}

  //if (!vtkSVPlanarMapper::MatrixVectorMultiply(Ainv, this->Bu, this->Xu))
  //{
  //  vtkErrorMacro("Error computing matrix vector multiply");
  //  return 0;
  //}

  //if (!vtkSVPlanarMapper::MatrixVectorMultiply(Ainv, this->Bv, this->Xv))
  //{
  //  vtkErrorMacro("Error computing matrix vector multiply");
  //  return 0;
  //}
  int numPoints = this->WorkPd->GetNumberOfPoints();

  svMath::conjugate_gradient(*this->ATutte, &this->Bu[0], numPoints, &this->Xu[0]);
  svMath::conjugate_gradient(*this->AHarm, &this->Bu[0], numPoints, &this->Xu[0]);
  svMath::conjugate_gradient(*this->ATutte, &this->Bv[0], numPoints, &this->Xv[0]);
  svMath::conjugate_gradient(*this->AHarm, &this->Bv[0], numPoints, &this->Xv[0]);

  for (int i=0; i<numPoints; i++)
  {
    double pt[3];
    pt[0] = this->Xu[i];
    pt[1] = this->Xv[i];
    pt[2] = 0.0;
    this->PlanarPd->GetPoints()->SetPoint(i, pt);
  }

  vtkNew(vtkPolyDataNormals, normaler);
  normaler->SetInputData(this->PlanarPd);
  normaler->SplittingOff();
  normaler->ComputePointNormalsOff();
  normaler->ComputeCellNormalsOn();
  normaler->Update();

  this->PlanarPd->DeepCopy(normaler->GetOutput());
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPlanarMapper::InvertSystem(std::vector<std::vector<double> > &mat,
                                  std::vector<std::vector<double> > &invMat)
{
  int nr = mat.size();
  int nc = mat[0].size();
  if (nr != nc)
  {
    //vtkErrorMacro("Matrix is not square");
    return 0;
  }

  double **inMat  = new double*[nr];
  double **outMat = new double*[nr];

  for (int i=0; i<nr; i++)
  {
    inMat[i]  = new double[nc];
    outMat[i] = new double[nc];
  }

  for (int i=0; i<nr; i++)
  {
    for (int j=0; j<nc; j++)
    {
      inMat[i][j] = mat[i][j];
    }
  }

  if (vtkMath::InvertMatrix(inMat, outMat, nr) == 0)
  {
    for (int i=0; i<nr; i++)
    {
      delete [] inMat[i];
      delete [] outMat[i];
    }
    delete [] inMat;
    delete [] outMat;
    //vtkErrorMacro("vtkMath could not invert matrix");
    return 0;
  }

  for (int i=0; i<nr; i++)
  {
    for (int j=0; j<nr; j++)
    {
      invMat[i][j] = outMat[i][j];
    }
  }

  for (int i=0; i<nc; i++)
  {
    delete [] inMat[i];
    delete [] outMat[i];
  }
  delete [] inMat;
  delete [] outMat;
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPlanarMapper::MatrixVectorMultiply(std::vector<std::vector<double> > &mat,
                                          std::vector<double> &inVec,
                                          std::vector<double> &outVec)
{
  int nrM = mat.size();;
  int ncM = mat[0].size();
  int nrV = inVec.size();
  if (ncM != nrV)
  {
    //vtkErrorMacro("Matrix and vector do not have consistent dimension");
    return 0;
  }

  outVec.resize(nrM, 0.0);
  for (int i=0; i<nrM; i++)
  {
    double updateVal = 0.0;
    for (int j=0; j<ncM; j++)
    {
      updateVal = mat[i][j]*inVec[j];
    }
    outVec[i] = updateVal;
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPlanarMapper::PrintMatrix(std::vector<std::vector<double> > &mat)
{
  int nr = mat.size();
  int nc = mat[0].size();
  fprintf(stdout,"Matrix: %d by %d\n", nr, nc);
  fprintf(stdout,"----------------------------------------------------------\n");
  for (int i=0; i<nr; i++)
  {
    for (int j=0; j<nc; j++)
    {
      if (mat[i][j] != 0.0)
      {
        fprintf(stdout,"| Row: %d Col: %d ", i, j);
        fprintf(stdout,"%.4f ", mat[i][j]);
        fprintf(stdout,"|");
      }
    }
    fprintf(stdout,"\n");
  }
  fprintf(stdout,"----------------------------------------------------------\n");

  return 1;
}

