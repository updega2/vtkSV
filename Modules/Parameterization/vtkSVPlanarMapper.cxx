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
#include "vtkSVGeneralUtils.h"
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
    return SV_ERROR;
  }

  if (this->RunFilter() != 1)
  {
    vtkErrorMacro("Error when mapping");
    output->DeepCopy(this->InitialPd);
    return SV_ERROR;
  }

  if (this->RemoveInternalIds)
  {
    this->WorkPd->GetPointData()->RemoveArray(this->InternalIdsArrayName);
    this->WorkPd->GetCellData()->RemoveArray(this->InternalIdsArrayName);
  }
  output->DeepCopy(this->PlanarPd);
  return SV_OK;
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
    return SV_ERROR;
  }

  //Check the input to make sure it is manifold and a triangulated surface
  if (vtkSVGeneralUtils::CheckSurface(this->InitialPd) != 1)
  {
    vtkErrorMacro("Error when checking input surface");
    return SV_ERROR;
  }

  // Check if internal id array name is given
  if (!this->InternalIdsArrayName)
  {
    vtkDebugMacro("Internal Ids Array Name not given, setting to InternalIds");
    this->InternalIdsArrayName = new char[strlen("InternalIds") + 1];
    strcpy(this->InternalIdsArrayName, "InternalIds");
  }
  // Check if array internal ids is already on pd
  if (vtkSVGeneralUtils::CheckArrayExists(this->WorkPd, 0, this->InternalIdsArrayName))
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
  if (!vtkSVGeneralUtils::CreateEdgeTable(this->WorkPd, this->EdgeTable, this->EdgeWeights,
                             this->EdgeNeighbors, this->IsBoundary))
  {
    vtkErrorMacro("Could not create edge table");
    return SV_ERROR;
  }

  // Set the size of the matrices
  this->ATutte = new SparseMatrix(numPoints, numPoints);
  this->AHarm = new SparseMatrix(numPoints, numPoints);
  this->Xu.resize(numPoints, 0.0);
  this->Xv.resize(numPoints, 0.0);
  this->Bu.resize(numPoints, 0.0);
  this->Bv.resize(numPoints, 0.0);

  return SV_OK;
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
    return SV_ERROR;
  }

  if (this->SetInternalNodes() != 1)
  {
    vtkErrorMacro("Error setting internal nodes");
    return SV_ERROR;
  }

  if (this->SolveSystem() != 1)
  {
    vtkErrorMacro("Error solving system");
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
int vtkSVPlanarMapper::SetBoundaries()
{
  this->BoundaryMapper->SetInputData(this->WorkPd);
  this->BoundaryMapper->SetEdgeTable(this->EdgeTable);
  this->BoundaryMapper->SetInternalIdsArrayName(this->InternalIdsArrayName);
  this->BoundaryMapper->Update();

  vtkNew(vtkPolyData, boundaryPd);
  boundaryPd->DeepCopy(this->BoundaryMapper->GetOutput());
  // Check if array internal ids is already on pd
  if (vtkSVGeneralUtils::CheckArrayExists(boundaryPd, 0, this->InternalIdsArrayName) == 0)
  {
    vtkErrorMacro("No internal ids array name on boundary pd");
    return SV_ERROR;
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

  return SV_OK;
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

  return SV_OK;
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
      vtkSVGeneralUtils::GetPointNeighbors(i, this->WorkPd, pointNeighbors);
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

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPlanarMapper::SolveSystem()
{
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
  return SV_OK;
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
    return SV_ERROR;
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
    return SV_ERROR;
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
  return SV_OK;
}
