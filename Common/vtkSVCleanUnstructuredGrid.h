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

class vtkPointLocator;

class VTKSVCOMMON_EXPORT vtkSVCleanUnstructuredGrid
  : public vtkUnstructuredGridAlgorithm
{
public:
  static vtkSVCleanUnstructuredGrid* New();

  vtkTypeMacro(vtkSVCleanUnstructuredGrid, vtkUnstructuredGridAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

protected:
  vtkSVCleanUnstructuredGrid();
  ~vtkSVCleanUnstructuredGrid();

  vtkPointLocator* Locator;

  virtual int RequestData(
    vtkInformation*, vtkInformationVector**, vtkInformationVector*) VTK_OVERRIDE;
  virtual int FillInputPortInformation(int port, vtkInformation* info) VTK_OVERRIDE;

private:
  vtkSVCleanUnstructuredGrid(const vtkSVCleanUnstructuredGrid&) VTK_DELETE_FUNCTION;
  void operator=(const vtkSVCleanUnstructuredGrid&) VTK_DELETE_FUNCTION;
};
#endif
