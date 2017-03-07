/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSVControlGrid.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSVControlGrid - topologically regular array of data
// .SECTION Description
//
// Inherets from vtkStructuredGrid

#ifndef vtkSVControlGrid_h
#define vtkSVControlGrid_h

#include "vtkStructuredGrid.h"

#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"

class vtkSVControlGrid : public vtkStructuredGrid
{
public:
  static vtkSVControlGrid *New();

  vtkTypeMacro(vtkSVControlGrid,vtkStructuredGrid);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Copy the geometric and topological structure of an input poly data object.
  void CopyStructure(vtkDataSet *ds);

  void Initialize();

  int SetControlPoint(const int i, const int j, const int k, const double p[3], const double w);
  int SetControlPoint(const int i, const int j, const int k, const double pw[4]);
  int InsertControlPoint(const int i, const int j, const int k, const double p[3], const double w);
  int InsertControlPoint(const int i, const int j, const int k, const double pw[4]);
  int GetControlPoint(const int i, const int j, const int k, double p[3], double &w);
  int GetControlPoint(const int i, const int j, const int k, double pw[4]);
  int GetPointId(const int i, const int j, const int k, int &ptId);

  // Description:
  // Get dimensions of this structured points dataset.
  virtual int *GetDimensions () {return vtkStructuredGrid::GetDimensions();}
  virtual void GetDimensions (int dim[3]) {vtkStructuredGrid::GetDimensions(dim);}

  // Description:
  // Retrieve an instance of this class from an information object.
  static vtkSVControlGrid* GetData(vtkInformation* info);
  static vtkSVControlGrid* GetData(vtkInformationVector* v, int i=0);

protected:
  vtkSVControlGrid();
  ~vtkSVControlGrid();

  virtual void ComputeScalarRange() {vtkStructuredGrid::GetScalarRange();}

  vtkPoints      *InternalPoints;
  vtkDoubleArray *InternalWeights;

private:
  vtkSVControlGrid(const vtkSVControlGrid&);  // Not implemented.
  void operator=(const vtkSVControlGrid&);  // Not implemented.
};

#endif






