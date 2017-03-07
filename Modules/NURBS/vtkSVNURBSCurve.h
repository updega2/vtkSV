/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSVNURBSCurve.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSVNURBSCurve - topologically regular array of data
// .SECTION Description
//
// Inherets from vtkDataObject

#ifndef vtkSVNURBSCurve_h
#define vtkSVNURBSCurve_h

#include "vtkDataObject.h"

#include "vtkCellArray.h"
#include "vtkSVControlGrid.h"
#include "vtkDenseArray.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkPolyData.h"

class vtkSVNURBSCurve : public vtkDataObject
{
public:
  static vtkSVNURBSCurve *New();
  vtkSVNURBSCurve(int m, vtkPoints *controlPoints, int n, vtkDoubleArray *knotPoints, int deg) {;}
  vtkSVNURBSCurve(int m, vtkPoints *controlPoints, vtkDoubleArray *knotPoints, vtkIntArray *knotMultiplicity, int deg) {;}

  vtkTypeMacro(vtkSVNURBSCurve,vtkDataObject);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Get and set the number of control points for curve
  vtkGetMacro(NumberOfControlPoints, int);
  vtkSetMacro(NumberOfControlPoints, int);

  // Description:
  // Get and set the number of knot points for curve
  vtkGetMacro(NumberOfKnotPoints, int);
  vtkSetMacro(NumberOfKnotPoints, int);

  // Description:
  // Get and set the knot vector object
  vtkGetObjectMacro(ControlPointGrid, vtkSVControlGrid);

  // Description:
  // Get and set the knot vector object
  vtkGetObjectMacro(KnotVector, vtkDoubleArray);

  // Description:
  // Get the PolyData Representation
  vtkGetObjectMacro(CurveRepresentation, vtkPolyData);

  void Initialize();

  //PolyData representation functions
  int GeneratePolyDataRepresentation(const double spacing);
  int GetStructuredGridConnectivity(const int numPoints, vtkCellArray *connectivity);

  //Functions to set control points/knots/etc.
  void SetControlPointGrid(vtkSVControlGrid *controlPoints) {}
  void SetControlPoints(vtkPoints *points1d);
  void SetControlPoints(vtkDenseArray<double[3]> *points2d) {}
  void SetKnotVector(vtkDoubleArray *knotVector);

  //Functions to manipulate the geometry
  void UpdateCurve() {}
  int IncreaseDegree(const int degree) {return 0;}
  int InsertKnot(const double newKnot, const double tolerance) {return 0;}
  int InsertKnots(vtkDoubleArray *newKnots, const double tolerance) {return 0;}
  int RemoveKnot(const int index, const double tolerance) {return 0;}
  int SetKnot(const int index, const double newKnot) {return 0;}
  int SetKnots(vtkIntArray *indices, vtkDoubleArray *newKnots) {return 0;}
  int GetKnot(const int index, double &knotVal) {return 0;}
  int GetKNots(const int indices, vtkDoubleArray *knotVals) {return 0;}

  int SetControlPoint(const int index, const double coordinate[3], const double weight) {return 0;}
  int SetControlPoints(vtkIntArray *indices, vtkPoints *coordinates, vtkDoubleArray *weights) {return 0;}
  int GetControlPoint(const int index, double coordinates[3], double &weight) {return 0;}
  int GetControlPoints(vtkIntArray *indices, vtkPoints *coordinates, vtkDoubleArray *weights) {return 0;}

  int SetWeight(const int index, const double weight) {return 0;}
  int GetWeight(const int index, double &weight) {return 0;}

  void SetClosed(const int closed) {this->Closed = closed;}
  void SetClamped(const int clamped) {this->Clamped = clamped;}
  int MakePeriodic(const int continuity) {return 0;}

  int GetMultiplicity(vtkIntArray *multiplicity, vtkDoubleArray *singleKnots);

  // Description:
  // Retrieve an instance of this class from an information object.
  static vtkSVNURBSCurve* GetData(vtkInformation* info);
  static vtkSVNURBSCurve* GetData(vtkInformationVector* v, int i=0);

protected:
  vtkSVNURBSCurve();
  ~vtkSVNURBSCurve();

  int NumberOfControlPoints;
  int NumberOfKnotPoints;
  int Clamped;
  int Closed;
  int Degree;

  vtkSVControlGrid *ControlPointGrid;
  vtkDoubleArray *KnotVector;
  vtkDoubleArray *Weights;

  vtkPolyData *CurveRepresentation;

private:
  vtkSVNURBSCurve(const vtkSVNURBSCurve&);  // Not implemented.
  void operator=(const vtkSVNURBSCurve&);  // Not implemented.
};

#endif
