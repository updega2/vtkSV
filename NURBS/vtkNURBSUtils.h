/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkNURBSUtils.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================
  Copyright 2011 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

  Contact: pppebay@sandia.gov,dcthomp@sandia.gov

=========================================================================*/
// .NAME vtkNURBSUtils - performs common math operations
// .SECTION Description
// vtkNURBSUtils provides methods to perform common math operations. These
// include providing constants such as Pi; conversion from degrees to
// radians; vector operations such as dot and cross products and vector
// norm; matrix determinant for 2x2 and 3x3 matrices; univariate polynomial
// solvers; and for random number generation (for backward compatibility only).
// .SECTION See Also
// vtkMinimalStandardRandomSequence, vtkBoxMuellerRandomSequence,
// vtkQuaternion

#ifndef vtkNURBSUtils_h
#define vtkNURBSUtils_h

#include "vtkObject.h"

#include "vtkControlGrid.h"
#include "vtkDenseArray.h"
#include "vtkIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkStructuredGrid.h"
#include "vtkTypedArray.h"

#include <cassert> // assert() in inline implementations.

class vtkNURBSUtils : public vtkObject
{
public:
  static vtkNURBSUtils *New();
  vtkTypeMacro(vtkNURBSUtils,vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);

  static int GetUs(vtkPoints *xyz, std::string type, vtkDoubleArray *U);
  static int LinSpace(double min, double max, int num, vtkDoubleArray *U);
  static int LinSpaceClamp(double min, double max, int num, int p, vtkDoubleArray *U);
  static int GetAvgKnots(double min, double max, int num, int p, vtkDoubleArray *U,
                         vtkDoubleArray *knots);
  static int GetEndDerivKnots(double min, double max, int num, int p, vtkDoubleArray *U,
                         vtkDoubleArray *knots);
  static int GetChordSpacedUs(vtkPoints *xyz, int num, vtkDoubleArray *U);
  static int GetCentripetalSpacedUs(vtkPoints *xyz, int num, vtkDoubleArray *U);
  static int GetZeroBasisFunctions(vtkDoubleArray *U, vtkDoubleArray *knots,
                                   vtkTypedArray<double> *N0);
  static int GetPBasisFunctions(vtkDoubleArray *u, vtkDoubleArray *knots,
                                const int p,
                                vtkTypedArray<double> *N);
  static int GetControlPointsOfCurve(vtkPoints *points, vtkDoubleArray *U,
                                     vtkDoubleArray *weights, vtkDoubleArray *knots,
                                     const int p,
                                     std::string ktype,
                                     const double D0[3], const double DN[3],
                                     vtkPoints *cPoints);
  static int GetControlPointsOfSurface(vtkStructuredGrid *points, vtkDoubleArray *U,
                                       vtkDoubleArray *V, vtkDoubleArray *uWeights,
                                       vtkDoubleArray *vWeights, vtkDoubleArray *uKnots,
                                       vtkDoubleArray *vKnots, const int p, const int q,
                                       std::string kutype, std::string kvtype,
                                       vtkDoubleArray *DU0, vtkDoubleArray *DUN,
                                       vtkDoubleArray *DV0, vtkDoubleArray *DVN,
                                       vtkStructuredGrid *cPoints);
  static int SetCurveEndDerivatives(vtkTypedArray<double> *NP, vtkTypedArray<double> *points,
		                                const int p, const double D0[3], const double DN[3],
                                    vtkDoubleArray *U, vtkDoubleArray *knots,
                                    vtkTypedArray<double> *newNP, vtkTypedArray<double> *newPoints);
  static int SetSurfaceEndDerivatives(vtkTypedArray<double> *NPU, vtkTypedArray<double> *NPV,
                                      vtkTypedArray<double> *points,
		                                  const int p, const int q,
                                      std::string kutype, std::string kvtype,
                                      vtkTypedArray<double> *DU0, vtkTypedArray<double> *DUN,
                                      vtkTypedArray<double> *DV0, vtkTypedArray<double> *DVN,
                                      vtkDoubleArray *U, vtkDoubleArray *V,
                                      vtkDoubleArray *uKnots, vtkDoubleArray *vKnots,
                                      vtkTypedArray<double> *newNPU, vtkTypedArray<double> *newNPV,
                                      vtkTypedArray<double> *newPoints);
  static int AddDerivativeRows(vtkTypedArray<double> *NP, vtkTypedArray<double> *newNP,
                               const int p, vtkDoubleArray *knots);
  static int AddDerivativePoints(vtkTypedArray<double> *points,
		                             const int p, const double D0[3],
                                 const double DN[3], vtkDoubleArray *U, vtkDoubleArray *knots,
                                 vtkTypedArray<double> *newPoints);
  static int GetKnots(vtkDoubleArray *u, int p, std::string type,
                      vtkDoubleArray *knots);
  static int InvertSystem(vtkTypedArray<double> *NP, vtkTypedArray<double> *NPinv);

  static int BasisEvaluation(vtkDoubleArray *knots, int p, int kEval, double uEval,
                             vtkDoubleArray *Nu);
  static int BasisEvaluationVec(vtkDoubleArray *knots, int p, int kEval, vtkDoubleArray *uEvals,
                             vtkTypedArray<double> *Nus);
  static int FindSpan(int p, double u, vtkDoubleArray *knots, int &span);

  //Conversion functions
  static int PolyDatasToStructuredGrid(vtkPolyData **inputs, const int numInputs, vtkStructuredGrid *points);
  static int StructuredGridToTypedArray(vtkStructuredGrid *grid, vtkTypedArray<double> *output);
  static int TypedArrayToStructuredGrid(vtkTypedArray<double> *array, vtkStructuredGrid *output);
  static int PointsToTypedArray(vtkPoints *points, vtkTypedArray<double> *output);
  static int TypedArrayToPoints(vtkTypedArray<double> *array, vtkPoints *output);
  static int DoubleArrayToTypedArray(vtkDoubleArray *input, vtkTypedArray<double> *output);
  static int MatrixToVector(vtkTypedArray<double> *mat, double *matVec);
  static int VectorToMatrix(double *matVec, const int nr, const int nc, vtkTypedArray<double> *mat);
  static int PointMatrixToVectors(vtkTypedArray<double> *mat, double *matVecs[3]);
  static int VectorsToPointMatrix(double *matVecs[3], const int nr, const int nc, vtkTypedArray<double> *mat);
  static int DeepCopy(vtkTypedArray<double> *input, vtkTypedArray<double> *output);

  //Matrix and vector math
  static int MatrixPointsMultiply(vtkTypedArray<double> *mat, vtkPoints *pointVec, vtkPoints *output);
  static int MatrixVecMultiply(vtkTypedArray<double> *mat, const int matIsPoints,
                               vtkTypedArray<double> *vec, const int vecIsPoints,
                               vtkTypedArray<double> *output);
  static int MatrixMatrixMultiply(vtkTypedArray<double> *mat0, const int mat0IsPoints,
                                  vtkTypedArray<double> *mat1, const int mat1IsPoints,
                                  vtkTypedArray<double> *output);
  static int MatrixMatrixForDGEMM(vtkTypedArray<double> *mat0,
                                  vtkTypedArray<double> *mat1,
                                  vtkTypedArray<double> *output);
  static int PointMatrixPointMatrixForDGEMM(vtkTypedArray<double> *mat0,
                                            vtkTypedArray<double> *mat1,
                                            vtkTypedArray<double> *output);
  static int MatrixPointMatrixForDGEMM(vtkTypedArray<double> *mat0,
                                       vtkTypedArray<double> *mat1,
                                       vtkTypedArray<double> *output);
  static int PointMatrixMatrixForDGEMM(vtkTypedArray<double> *mat0,
                                       vtkTypedArray<double> *mat1,
                                       vtkTypedArray<double> *output);
  static int DGEMM(const double *A, const int nrA, const int ncA,
                   const double *B, const int nrB, const int ncB,
                   double *C);
  static int GetMatrixComp(vtkTypedArray<double> *mat,  const int loc, const int comp,
                           const int matIsPoints, vtkTypedArray<double> *vec);
  static int SetMatrixComp(vtkTypedArray<double> *vec,  const int loc, const int comp,
                           const int matIsPoints, vtkTypedArray<double> *mat);
  static int StructuredGridTranspose(vtkStructuredGrid *sg, vtkStructuredGrid *newSg);
  static int MatrixTranspose(vtkTypedArray<double> *mat, const int matIsPoints,
                             vtkTypedArray<double> *newMat);

  //Simple index-wise vector operations
  static int Add1D(vtkDoubleArray *v0, vtkDoubleArray *v1, double scalar, vtkDoubleArray *result);
  static int AddVal1D(vtkDoubleArray *v0, double val, double scalar, vtkDoubleArray *result);
  static int AddVal1D(double val, vtkDoubleArray *v0, double scalar, vtkDoubleArray *result);
  static int MultiplyVal1D(vtkDoubleArray *v0, double val, vtkDoubleArray *result);
  static int Intersect1D(vtkIntArray *v0, vtkIntArray *v1, vtkIntArray *result);
  static int WhereGreaterEqual(double val, vtkDoubleArray *in, vtkIntArray *out);
  static int WhereGreater(double val, vtkDoubleArray *in, vtkIntArray *out);
  static int WhereLessEqual(double val, vtkDoubleArray *in, vtkIntArray *out);
  static int WhereLess(double val, vtkDoubleArray *in, vtkIntArray *out);
  static int WhereEqual(double val, vtkDoubleArray *in, vtkIntArray *out);
  static int WhereNotEqual(double val, vtkDoubleArray *in, vtkIntArray *out);

  //Print operations
  static int PrintArray(vtkIntArray *arr);
  static int PrintArray(vtkDoubleArray *arr);
  static int PrintVector(vtkTypedArray<double> *vec);
  static int PrintMatrix(vtkTypedArray<double> *mat);
  static int PrintStructuredGrid(vtkStructuredGrid *mat);
  static int PrintPoints(vtkPoints *points);
  static int Print2DArray(const double *arr, const int nr, const int nc);

protected:
  vtkNURBSUtils();
  ~vtkNURBSUtils();

private:
  vtkNURBSUtils(const vtkNURBSUtils&);  // Not implemented.
  void operator=(const vtkNURBSUtils&);  // Not implemented.
};

#endif
