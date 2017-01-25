/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkIntersectionPolyDataFilter2.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkIntersectionPolyDataFilter2
// .SECTION Description
//
// vtkIntersectionPolyDataFilter2 computes the intersection between two
// vtkPolyData objects. The first output is a set of lines that marks
// the intersection of the input vtkPolyData objects. The second and
// third outputs are the first and second input vtkPolyData,
// respectively. Optionally, the two output vtkPolyData can be split
// along the intersection lines.
//
// This code was contributed in the Insight Journal paper:
// "Boolean Operations on Surfaces in VTK Without External Libraries"
// by Cory Quammen, Chris Weigle C., Russ Taylor
// http://hdl.handle.net/10380/3262
// http://www.insight-journal.org/browse/publication/797

#ifndef __vtkIntersectionPolyDataFilter2_h
#define __vtkIntersectionPolyDataFilter2_h

#include "vtkFiltersGeneralModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"

#define REAL double 

class VTKFILTERSGENERAL_EXPORT vtkIntersectionPolyDataFilter2 : public vtkPolyDataAlgorithm
{
public:
  static vtkIntersectionPolyDataFilter2 *New();
  vtkTypeMacro(vtkIntersectionPolyDataFilter2, vtkPolyDataAlgorithm);
  virtual void PrintSelf(ostream &os, vtkIndent indent);

  // Description:
  // Integer describing the number of intersection points and lines
  vtkGetMacro(NumberOfIntersectionPoints, int);
  vtkSetMacro(NumberOfIntersectionPoints, int);
  vtkBooleanMacro(NumberOfIntersectionPoints, int);
  vtkGetMacro(NumberOfIntersectionLines, int);
  vtkSetMacro(NumberOfIntersectionLines, int);
  vtkBooleanMacro(NumberOfIntersectionLines, int);

  // Description:
  // If on, the second output will be the first input mesh split by the
  // intersection with the second input mesh. Defaults to on.
  vtkGetMacro(SplitFirstOutput, int);
  vtkSetMacro(SplitFirstOutput, int);
  vtkBooleanMacro(SplitFirstOutput, int);

  // Description:
  // If on, the third output will be the second input mesh split by the
  // intersection with the first input mesh. Defaults to on.
  vtkGetMacro(SplitSecondOutput, int);
  vtkSetMacro(SplitSecondOutput, int);
  vtkBooleanMacro(SplitSecondOutput, int);

  // Description:
  // If on, an array containing the boundary points will be attached to the 
  // outputs
  vtkGetMacro(ApplyBoundaryPointArray, int);
  vtkSetMacro(ApplyBoundaryPointArray, int);
  vtkBooleanMacro(ApplyBoundaryPointArray, int);

  // Description:
  // If on, the normals of the input will be checked.
  vtkGetMacro(CheckInput, int);
  vtkSetMacro(CheckInput, int);
  vtkBooleanMacro(CheckInput, int);

  // Description:
  // If on, the third output will be the second input mesh split by the
  // intersection with the first input mesh. Defaults to on.
  vtkGetMacro(CheckMesh, int);
  vtkSetMacro(CheckMesh, int);
  vtkBooleanMacro(CheckMesh, int);

  // Description:
  vtkGetMacro(Tolerance, double);
  vtkSetMacro(Tolerance, double);

  // Description:
  // Given two triangles defined by points (p1, q1, r1) and (p2, q2,
  // r2), returns whether the two triangles intersect. If they do,
  // the endpoints of the line forming the intersection are returned
  // in pt1 and pt2. The parameter coplanar is set to 1 if the
  // triangles are coplanar and 0 otherwise.
  static int TriangleTriangleIntersection(double p1[3], double q1[3], double r1[3],
                                          double p2[3], double q2[3], double r2[3],
                                          int &coplanar, double pt1[3], double pt2[3],
					  double surfaceid[2],double tolerance,
					  int verbose);

  static void CleanAndCheckSurface(vtkPolyData *pd,double tol);
  static void CleanAndCheckInput(vtkPolyData *pd);


protected:
  vtkIntersectionPolyDataFilter2();
  ~vtkIntersectionPolyDataFilter2();

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
  int FillInputPortInformation(int, vtkInformation*);

private:
  vtkIntersectionPolyDataFilter2(const vtkIntersectionPolyDataFilter2&); // Not implemented
  void operator=(const vtkIntersectionPolyDataFilter2&); // Not implemented

  int NumberOfIntersectionPoints;
  int NumberOfIntersectionLines;
  int SplitFirstOutput;
  int SplitSecondOutput;
  int ApplyBoundaryPointArray;
  int CheckMesh;
  int CheckInput;
  double Tolerance;

  class Impl;
};

REAL orient2d(REAL *pa, REAL *pb, REAL *pc);

#endif // __vtkIntersectionPolyDataFilter2_h
