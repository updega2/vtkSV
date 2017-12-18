/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSVBestFitPlane.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSVBestFitPlane - Find the least squares plane through a set of points
// .SECTION Description
// vtkSVBestFitPlane finds the least squares plane through a set of points.

#ifndef vtkSVBestFitPlane_h
#define vtkSVBestFitPlane_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkSVFiltersModule.h" // For export

class VTKSVFILTERS_EXPORT vtkSVBestFitPlane : public vtkPolyDataAlgorithm
{
public:
  static vtkSVBestFitPlane *New();
  vtkTypeMacro(vtkSVBestFitPlane,vtkPolyDataAlgorithm);

protected:

  vtkSVBestFitPlane();

  int RequestData(vtkInformation* request,
                  vtkInformationVector** inputVector,
                  vtkInformationVector* outputVector );
  int FillInputPortInformation(int port, vtkInformation* info);

  enum WeightEnum {UNIFORMWEIGHT, GAUSSIANWEIGHT};

  // Description:
  // Compute a weight based on the WeightMode.
  double WeightFunction(double distance);

  // Description:
  // Set the WeightMode to UNIFORMWEIGHT
  void SetWeightModeToUniform(){this->WeightMode = UNIFORMWEIGHT;}

  // Description:
  // Set the WeightMode to GAUSSIANWEIGHT
  void SetWeightModeToGaussian(){this->WeightMode = GAUSSIANWEIGHT;}

  // Description:
  // A flag specifying how to weight points.
  int WeightMode;

  // Description:
  // Accessor/mutator for GaussianVariance
  vtkSetMacro(GaussianVariance, double);
  vtkGetMacro(GaussianVariance, double);

  // Description:
  // The variance of the Gaussian function to use if
  // WeightMode == GAUSSIANWEIGHT
  double GaussianVariance;
private:

  vtkSVBestFitPlane(const vtkSVBestFitPlane&);  // Not implemented.
  void operator=(const vtkSVBestFitPlane&);  // Not implemented.
};

#endif
