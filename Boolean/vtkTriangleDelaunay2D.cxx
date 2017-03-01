/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTriangleDelaunay2D.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkTriangleDelaunay2D.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkExecutive.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkTransform.h>

extern "C"
{
#define ANSI_DECLARATORS
#ifndef VOID
# define VOID void
#endif
#ifndef REAL
# define REAL double
#endif
#include "triangle.h"
}
#include "predicates.h"

#include <algorithm>


vtkStandardNewMacro(vtkTriangleDelaunay2D);

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
vtkTriangleDelaunay2D::vtkTriangleDelaunay2D()
{
  // This filters needs generators and boundary
  this->SetNumberOfInputPorts(2);

  // This filters produces a constrained Delaunay triangulation
  this->SetNumberOfOutputPorts(1);
}

//-----------------------------------------------------------------------------
vtkTriangleDelaunay2D::~vtkTriangleDelaunay2D()
{
}

//-----------------------------------------------------------------------------
void vtkTriangleDelaunay2D::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//-----------------------------------------------------------------------------
int vtkTriangleDelaunay2D::RequestData(vtkInformation* vtkNotUsed(request),
                                       vtkInformationVector** inputVector,
                                       vtkInformationVector* outputVector)
{
  // Extract inputs
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0], 0);
  vtkPolyData* source = vtkPolyData::GetData(inputVector[1], 0);

  // Extract outputs
  vtkPolyData* cdt = vtkPolyData::GetData(outputVector, 0);

  // Initialize; check input
  //
  vtkPoints* inPoints = input->GetPoints();
  if (inPoints == NULL)
    {
    vtkDebugMacro("Cannot triangulate; no input points");
    return 1;
    }

  vtkIdType numPoints = inPoints->GetNumberOfPoints();
  if (numPoints < 3)
    {
    vtkDebugMacro("Cannot triangulate; need at least 3 input points");
    return 1;
    }

  vtkPoints *points = inPoints;
  vtkPoints *tPoints = NULL;
  // If the user specified a transform, apply it to the input data.
  //
  // Only the input points are transformed.  We do not bother
  // transforming the source points (if specified).  The reason is
  // that only the topology of the Source is used during the constrain
  // operation.  The point ids in the Source topology are assumed to
  // reference points in the input. So, when an input transform is
  // used, only the input points are transformed.  We do not bother
  // with transforming the Source points since they are never
  // referenced.
  if (this->Transform)
    {
    tPoints = vtkPoints::New();
    this->Transform->TransformPoints(inPoints, tPoints);
    points = tPoints;
    }
  else
    {
    // If the user asked this filter to compute the best fitting plane,
    // proceed to compute the plane and generate a transform that will
    // map the input points into that plane.
    if (this->ProjectionPlaneMode == VTK_BEST_FITTING_PLANE)
      {
      this->SetTransform(this->ComputeBestFittingPlane(input));
      tPoints = vtkPoints::New();
      this->Transform->TransformPoints(inPoints, tPoints);
      points = tPoints;
      }
    }

  struct triangulateio in, out;
  memset(&in, 0, sizeof(struct triangulateio));
  memset(&out, 0, sizeof(struct triangulateio));

  // Copy point coordinates to triangle input data structure
  std::size_t nbPts = (std::size_t)points->GetNumberOfPoints();
  in.numberofpoints = static_cast<int>(nbPts);
  in.pointlist = new REAL[in.numberofpoints * 2];
  if (points->GetDataType() == VTK_FLOAT)
    {
    float* ptsPtr = static_cast<float*>(points->GetVoidPointer(0));
    for (std::size_t i = 0, ptr = 0; i < nbPts; i++, ptsPtr++)
      {
      in.pointlist[ptr++] = *ptsPtr++;
      in.pointlist[ptr++] = *ptsPtr++;
      }
    }
  else if (points->GetDataType() == VTK_DOUBLE)
    {
    double* ptsPtr = static_cast<double*>(points->GetVoidPointer(0));
    for (std::size_t i = 0, ptr = 0; i < nbPts; i++, ptsPtr++)
      {
      in.pointlist[ptr++] = *ptsPtr++;
      in.pointlist[ptr++] = *ptsPtr++;
      }
    }

  unsigned int nbSegments = 0;

  std::vector<bool> polyIsInternal;
  vtkIdType nbOfInternalPolys = 0;

  // Process polys then lines
  if (source)
    {
    polyIsInternal.reserve(source->GetNumberOfPolys());
    vtkCellArray* cell = source->GetPolys();
    for (int sourceCnt = 0; sourceCnt < 2; sourceCnt++)
      {
      if (cell)
        {
        cell->InitTraversal();
        vtkIdType nbpts, *pts;
        while (cell->GetNextCell(nbpts, pts))
          {
          nbSegments += nbpts;
          if (sourceCnt == 0)
            {
            bool isInternal = !this->IsCellCCW(nbpts, pts, points);
            polyIsInternal.push_back(isInternal);
            if (isInternal)
              {
              nbOfInternalPolys++;
              }
            }
          }
        }
      cell = source->GetLines();
      }
    }

  // add MA edges, if any
  in.numberofsegments = nbSegments;
  in.segmentlist = nbSegments > 0 ? new int[in.numberofsegments * 2] : 0;
  in.segmentmarkerlist = nbSegments > 0 ? new int[in.numberofsegments] : 0;
  in.numberofholes = nbOfInternalPolys;

  if (nbOfInternalPolys > 0)
    {
    in.holelist = new REAL[nbOfInternalPolys * 2 + 1];
    }

  vtkIdType cnt = 0;
  if (source)
    {
    unsigned int segCnt = 0, holecnt = 0;
    vtkCellArray* cell = source->GetPolys();
    for (int sourceCnt = 0; sourceCnt < 2; sourceCnt++)
      {
      if (cell)
        {
        cell->InitTraversal();
        vtkIdType nbpts, *pts;
        while (cell->GetNextCell(nbpts, pts))
          {
          nbpts -= (sourceCnt == 0) ? 0 : 1;
          unsigned int initSeg = segCnt;
          for (int j = 0; j < nbpts; j++)
            {
            in.segmentlist[2 * segCnt + 0] = pts[j];
            in.segmentlist[2 * segCnt + 1] =
              pts[(sourceCnt == 0 && j == nbpts - 1) ? 0 : (j+1)];
            in.segmentmarkerlist[segCnt] = cnt + 1;
            segCnt++;
            }
          // do not consider external boundary & lines
          if (sourceCnt == 0 && polyIsInternal[cnt])
            {
            // Compute hole position
            this->GetPointInsidePolygon(in.pointlist,
              &in.segmentlist[2 * initSeg], nbpts,
              &in.holelist[holecnt++ * 2]);
            }
          cnt++;
          }
        }
      cell = source->GetLines();
      }
    }

  triangulate((char*)(in.numberofsegments > 0 ? "pzQ" : "zQ"),
    &in, &out, (struct triangulateio*)NULL);

  vtkDebugMacro(
    << "Triangle CDT has " << out.numberoftriangles << " triangles");

  // Fetch output triangles
  vtkNew<vtkCellArray> tri;
  tri->Allocate(4 * out.numberoftriangles);

  for (int i = 0; i < out.numberoftriangles; i++)
    {
    vtkIdType ids[3] =
      {
      out.trianglelist[i * 3 + 0],
      out.trianglelist[i * 3 + 1],
      out.trianglelist[i * 3 + 2]
      };
    tri->InsertNextCell(3, ids);
    }

  cdt->SetPoints(points);
  cdt->SetPolys(tri.Get());
  cdt->GetPointData()->ShallowCopy(input->GetPointData());

  // Free all allocated arrays, including the one allocated by Triangle
  delete [] in.pointlist;
  delete [] in.segmentlist;
  delete [] in.segmentmarkerlist;
  delete [] in.holelist;
  trifree(out.trianglelist);
  trifree(out.pointlist);
  trifree(out.pointmarkerlist);
  trifree(out.segmentlist);
  trifree(out.segmentmarkerlist);

  // If the best fitting option was ON, then the current transform
  // is the one that was computed internally. We must now destroy it.
  if (this->ProjectionPlaneMode == VTK_BEST_FITTING_PLANE)
    {
    if (this->Transform)
      {
      this->Transform->UnRegister(this);
      this->Transform = NULL;
      }
    }
  if (tPoints)
    {
    tPoints->Delete();
    }

  return 1;
}

//-----------------------------------------------------------------------------
// Return a point inside the polygon described by the provided segment list
// made of nbptsinseg point pairs. pts is the 2D point coordinate array.
void vtkTriangleDelaunay2D::GetPointInsidePolygon(double* pts, int* segs,
  vtkIdType nbptsinseg, double* P)
{
  assert(pts && segs && P);
  double V[2], A[2], B[2], Q[2];
  double minDist = VTK_DOUBLE_MAX;
  int iv = -1;

  // Search the first convex point in the segment strip
  for (vtkIdType i = 1; i < nbptsinseg - 1; i++)
    {
    int a = segs[(i-1) * 2];
    int v = segs[ i    * 2];
    int b = segs[(i+1) * 2];
    // AVB in clockwise order
    if (p_orient2d(&pts[a*2], &pts[v*2], &pts[b*2]) < 0.)
      {
      memcpy(V, &pts[v*2], 2 * sizeof(double));
      memcpy(A, &pts[a*2], 2 * sizeof(double));
      memcpy(B, &pts[b*2], 2 * sizeof(double));
      iv = i;
      break;
      }
    }
  if (iv == -1) return;

  bool found = false;
  for (vtkIdType i = 1; i < nbptsinseg; i++)
    {
    if (i >= iv -1 && i <= iv + 1)
      {
      continue;
      }
    double q[2];
    int pt = segs[i*2];
    memcpy(q, &pts[pt * 2], 2 * sizeof(double));

    // point inside triangle test
    double a1 = p_orient2d(V, B, q);
    double a2 = p_orient2d(B, A, q);
    double a3 = p_orient2d(A, V, q);
    if ((a1 < 0. && a2 < 0. && a3 < 0.) ||
      (a1 > 0. && a2 > 0. && a3 > 0.))
      {
      double d = (V[0] - q[0]) * (V[0] - q[0]) +
        (V[1] - q[1]) * (V[1] - q[1]);
      // is point nearest from V than the previous one?
      if (d < minDist)
        {
        // if so, remember it
        minDist = d;
        Q[0] = q[0];
        Q[1] = q[1];
        found = true;
        }
      }
    }

  if (!found)
    {
    // return the triangle center
    P[0] = (V[0] + A[0] + B[0]) / 3.;
    P[1] = (V[1] + A[1] + B[1]) / 3.;
    }
  else
    {
    // return V-Q midpoint
    P[0] = (V[0] + Q[0]) / 2.;
    P[1] = (V[1] + Q[1]) / 2.;
    }
}

//-----------------------------------------------------------------------------
// Compute if cell is CCW or CW using cell area.
bool vtkTriangleDelaunay2D::IsCellCCW(
  vtkIdType npts, vtkIdType* pts, vtkPoints* points)
{
  assert(pts && points);
  double area = 0.0;
  if (npts == 0)
    {
    return false;
    }
  double p1[3], p2[3];
  double* pa = p1, *pb = p2;
  points->GetPoint(pts[0], pa);
  for (int i = 0; i < npts; i++)
    {
    points->GetPoint(pts[(i + 1) % npts], pb);
    area += pa[0] * pb[1] - pb[0] * pa[1];
    // Swap pa & pb
    double *tmp = pa;
    pa = pb;
    pb = tmp;
    }
  return area * 0.5 >= 0.0;
}
