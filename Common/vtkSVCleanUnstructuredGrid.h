/*=========================================================================

  Program:   ParaView
  Module:    vtkSVCleanUnstructuredGrid.h

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/**
 * @class   vtkSVCleanUnstructuredGrid
 * @brief   merge duplicate points
 *
 *
 * vtkSVCleanUnstructuredGrid is a filter that takes unstructured grid data as
 * input and generates unstructured grid data as output. vtkSVCleanUnstructuredGrid can
 * merge duplicate points (with coincident coordinates) using the vtkMergePoints object
 * to merge points.
 *
 * @sa
 * vtkCleanPolyData
*/

#ifndef vtkSVCleanUnstructuredGrid_h
#define vtkSVCleanUnstructuredGrid_h

#include "vtkSVCommonModule.h" //needed for exports
#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkIncrementalPointLocator.h"

class VTKSVCOMMON_EXPORT vtkSVCleanUnstructuredGrid
  : public vtkUnstructuredGridAlgorithm
{
public:
  static vtkSVCleanUnstructuredGrid* New();

  vtkTypeMacro(vtkSVCleanUnstructuredGrid, vtkUnstructuredGridAlgorithm);
  //@{
  /**
   * Set/Get a spatial locator for speeding the search process. By
   * default an instance of vtkMergePoints is used.
   */
  vtkSetObjectMacro(Locator, vtkIncrementalPointLocator);
  vtkGetObjectMacro(Locator, vtkIncrementalPointLocator);
  //@}

  //@{
  /**
   * By default ToleranceIsAbsolute is false and Tolerance is
   * a fraction of Bounding box diagonal, if true, AbsoluteTolerance is
   * used when adding points to locator (merging)
   */
  vtkSetMacro(ToleranceIsAbsolute,int);
  vtkBooleanMacro(ToleranceIsAbsolute,int);
  vtkGetMacro(ToleranceIsAbsolute,int);
  //@}

  //@{
  /**
   * Specify tolerance in terms of fraction of bounding box length.
   * Default is 0.0.
   */
  vtkSetClampMacro(Tolerance,double,0.0,1.0);
  vtkGetMacro(Tolerance,double);
  //@}

  //@{
  /**
   * Specify tolerance in absolute terms. Default is 1.0.
   */
  vtkSetClampMacro(AbsoluteTolerance,double,0.0,VTK_DOUBLE_MAX);
  vtkGetMacro(AbsoluteTolerance,double);
  //@}

  /**
   * Create default locator. Used to create one when none is specified.
   */
  void CreateDefaultLocator(vtkDataSet *input = 0);

  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkSVCleanUnstructuredGrid();
  ~vtkSVCleanUnstructuredGrid();

  double Tolerance;
  double AbsoluteTolerance;
  int ToleranceIsAbsolute;

  vtkIncrementalPointLocator *Locator;

  virtual int RequestData(
    vtkInformation*, vtkInformationVector**, vtkInformationVector*);
  virtual int FillInputPortInformation(int port, vtkInformation* info);

private:
  vtkSVCleanUnstructuredGrid(const vtkSVCleanUnstructuredGrid&);
  void operator=(const vtkSVCleanUnstructuredGrid&);
};
#endif
