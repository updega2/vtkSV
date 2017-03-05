/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkNURBSSurface.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkNURBSSurface - topologically regular array of data
// .SECTION Description
//
// Inherets from vtkDataObject

#ifndef vtkNURBSSurface_h
#define vtkNURBSSurface_h

#include "vtkDataObject.h"

#include "vtkControlGrid.h"
#include "vtkDenseArray.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkPolyData.h"

class vtkNURBSSurface : public vtkDataObject
{
public:
  static vtkNURBSSurface *New();
  vtkNURBSSurface(int m, vtkPoints *controlPoints, int n, vtkDoubleArray *knotPoints, int deg) {;}
  vtkNURBSSurface(int m, vtkPoints *controlPoints, vtkDoubleArray *knotPoints, vtkIntArray *knotMultiplicity, int deg) {;}

  vtkTypeMacro(vtkNURBSSurface,vtkDataObject);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Get and set the number of control points for curve
  vtkGetMacro(NumberOfUControlPoints, int);
  vtkSetMacro(NumberOfUControlPoints, int);
  vtkGetMacro(NumberOfVControlPoints, int);
  vtkSetMacro(NumberOfVControlPoints, int);


  // Description:
  // Get and set the number of knot points for curve
  vtkGetMacro(NumberOfUKnotPoints, int);
  vtkSetMacro(NumberOfUKnotPoints, int);
  vtkGetMacro(NumberOfVKnotPoints, int);
  vtkSetMacro(NumberOfVKnotPoints, int);

  // Description:
  // Get and set the knot vector object
  vtkGetObjectMacro(ControlPointGrid, vtkControlGrid);

  // Description:
  // Get and set the knot vector object
  vtkGetObjectMacro(UKnotVector, vtkDoubleArray);
  vtkGetObjectMacro(VKnotVector, vtkDoubleArray);

  // Description:
  // Get the PolyData Representation
  vtkGetObjectMacro(SurfaceRepresentation, vtkPolyData);

  void Initialize();

  //PolyData representation functions
  int GeneratePolyDataRepresentation(const double uSpacing, const double vSpacing);

  //Functions to set control points/knots/etc.
  void SetControlPointGrid(vtkControlGrid *controlPoints) {}
  void SetControlPoints(vtkStructuredGrid *points2d);
  void SetKnotVector(vtkDoubleArray *knotVector, const int dim);

  //Functions to manipulate the geometry
  void UpdateCurve() {}
  int IncreaseDegree(const int degree, const int dim) {return 0;}
  int InsertKnot(const double newKnot, const int dim, const double tolerance) {return 0;}
  int InsertKnots(vtkDoubleArray *newKnots, const int dim, const double tolerance) {return 0;}
  int RemoveKnot(const int index, const int dim, const double tolerance) {return 0;}
  int SetKnot(const int index, const int dim, const double newKnot) {return 0;}
  int SetKnots(vtkIntArray *indices, const int dim, vtkDoubleArray *newKnots) {return 0;}
  int GetKnot(const int index, const int dim, double &knotVal) {return 0;}
  int GetKNots(const int indices, const int dim, vtkDoubleArray *knotVals) {return 0;}

  int SetControlPoint(const int index, const int dim, const double coordinate[3], const double weight) {return 0;}
  int SetControlPoints(vtkIntArray *indices, const int dim, vtkPoints *coordinates, vtkDoubleArray *weights);
  int GetControlPoint(const int index, const int dim, double coordinates[3], double &weight) {return 0;}
  int GetControlPoints(vtkIntArray *indices, const int dim, vtkPoints *coordinates, vtkDoubleArray *weights) {return 0;}

  int SetWeight(const int index, const int dim, const double weight) {return 0;}
  int GetWeight(const int index, const int dim, double &weight) {return 0;}

  void SetClosed(const int closed, const int dim) {;}
  void SetClamped(const int clamped, const int dim) {;}
  int MakePeriodic(const int continuity, const int dim) {return 0;}

  int GetMultiplicity(vtkIntArray *multiplicity, vtkDoubleArray *singleKnots);
  int GetStructuredGridConnectivity(const int numXPoints, const int numYPoints, vtkCellArray *connectivity);

  // Description:
  // Retrieve an instance of this class from an information object.
  static vtkNURBSSurface* GetData(vtkInformation* info);
  static vtkNURBSSurface* GetData(vtkInformationVector* v, int i=0);

protected:
  vtkNURBSSurface();
  ~vtkNURBSSurface();

  int NumberOfUControlPoints;
  int NumberOfVControlPoints;
  int NumberOfUKnotPoints;
  int NumberOfVKnotPoints;
  int UDegree;
  int VDegree;
  int UClamped;
  int VClamped;
  int UClosed;
  int VClosed;

  vtkControlGrid *ControlPointGrid;
  vtkDoubleArray *UKnotVector;
  vtkDoubleArray *VKnotVector;
  vtkDoubleArray *UVKnotVectors[2];
  vtkDoubleArray *UWeights;
  vtkDoubleArray *VWeights;
  vtkDoubleArray *UVWeights[2];

  vtkPolyData *SurfaceRepresentation;

private:
  vtkNURBSSurface(const vtkNURBSSurface&);  // Not implemented.
  void operator=(const vtkNURBSSurface&);  // Not implemented.
};

#endif
