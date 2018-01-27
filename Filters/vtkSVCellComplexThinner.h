/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSVCellComplexThinner.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#ifndef vtkSVCellComplexThinner_h
#define vtkSVCellComplexThinner_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkSVFiltersModule.h" // For export

class VTKSVFILTERS_EXPORT vtkSVCellComplexThinner : public vtkPolyDataAlgorithm
{
public:
  static vtkSVCellComplexThinner *New();
  vtkTypeMacro(vtkSVCellComplexThinner,vtkPolyDataAlgorithm);

  //@{
  /// \brief Get/set name for array on edge pd if preserving
  vtkGetStringMacro(PreserveEdgeCellsArrayName);
  vtkSetStringMacro(PreserveEdgeCellsArrayName);
  //@}

  //@{
  /// \brief Get/set name for edge pds
  vtkSetObjectMacro(InputEdgePd, vtkPolyData);
  vtkGetObjectMacro(InputEdgePd, vtkPolyData);

  vtkSetObjectMacro(OutputEdgePd, vtkPolyData);
  vtkGetObjectMacro(OutputEdgePd, vtkPolyData);
  //@}

protected:

  vtkSVCellComplexThinner();
  ~vtkSVCellComplexThinner();

  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *);

private:

  vtkPolyData *InputEdgePd;
  vtkPolyData *OutputEdgePd;

  vtkPolyData *WorkTriPd;
  vtkPolyData *WorkEdgePd;

  char *PreserveEdgeCellsArrayName;

  int PrepFilter();
  int RunFilter();

  vtkSVCellComplexThinner(const vtkSVCellComplexThinner&);  // Not implemented.
  void operator=(const vtkSVCellComplexThinner&);  // Not implemented.
};

#endif
