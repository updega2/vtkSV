/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkGeneralizedPolyCube.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkGeneralizedPolyCube - topologically regular array of data
// .SECTION Description
//
// Inherets from vtkDataObject

#ifndef vtkGeneralizedPolyCube_h
#define vtkGeneralizedPolyCube_h

#include "vtkUnstructuredGrid.h"

#include "vtkStructuredGrid.h"
#include "vtkDenseArray.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkPolyData.h"

class vtkGeneralizedPolyCube : public vtkUnstructuredGrid
{
public:
  static vtkGeneralizedPolyCube *New();
  vtkGeneralizedPolyCube(int m, vtkPoints *controlPoints, int n, vtkDoubleArray *knotPoints, int deg) {;}
  vtkGeneralizedPolyCube(int m, vtkPoints *controlPoints, vtkDoubleArray *knotPoints, vtkIntArray *knotMultiplicity, int deg) {;}

  vtkTypeMacro(vtkGeneralizedPolyCube,vtkUnstructuredGrid);
  void PrintSelf(ostream& os, vtkIndent indent);

  void Initialize();

  int GenerateGridsFromCenterLines() {return 0;}
  int AddGridWithOrigin(int cellId, double origin[], int boundary) {return 0;}
  int AddGridWithCenter(int cellId, double center[], int boundary) {return 0;}
  int AddGrid(int cellId, int boundary) {return 0;}
  int GenerateFullRepresentation() {return 0;}
  int GetGrid(int cellId, vtkStructuredGrid *gridRepresentation) {return 0;}


protected:
  vtkGeneralizedPolyCube();
  ~vtkGeneralizedPolyCube();

  vtkUnstructuredGrid *FullRepresentation;
  vtkIntArray         *Boundaries;

  vtkPolyData *Centerlines;

private:
  vtkGeneralizedPolyCube(const vtkGeneralizedPolyCube&);  // Not implemented.
  void operator=(const vtkGeneralizedPolyCube&);  // Not implemented.
};

#endif
