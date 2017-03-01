/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTriangleDelaunay2D.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// VTK class that wraps Shewchuck's Triangle delaunay implementation.
// Author: Joachim Pouderoux, Kitware SAS (2016)

#ifndef __vtkTriangleDelaunay2D_h
#define __vtkTriangleDelaunay2D_h

#include <vtkDelaunay2D.h>

class vtkTriangleDelaunay2D : public vtkDelaunay2D
{
public:
  vtkTypeMacro(vtkTriangleDelaunay2D, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkTriangleDelaunay2D* New();

protected:
  vtkTriangleDelaunay2D();
  ~vtkTriangleDelaunay2D();

  virtual int RequestData(vtkInformation*,
    vtkInformationVector**, vtkInformationVector*);

private:
  vtkTriangleDelaunay2D(const vtkTriangleDelaunay2D&);  // Not implemented.
  void operator=(const vtkTriangleDelaunay2D&);  // Not implemented.

  static void GetPointInsidePolygon(double* pts, int* segs,
                           vtkIdType nbptsinseg, double* P);
  static bool IsCellCCW(vtkIdType npts, vtkIdType* pts, vtkPoints* points);
};

#endif
