/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSuperSquareBoundaryMapper.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkSuperSquareBoundaryMapper
 * @brief   base class for boundary mapping filters
 *
 * vtkSuperSquareBoundaryMapper is an abstract class that defines
 * the protocol for suboundary mapping surface filters.
 *
*/

#ifndef vtkSuperSquareBoundaryMapper_h
#define vtkSuperSquareBoundaryMapper_h

#include "vtkFiltersGeneralModule.h" // For export macro
#include "vtkBoundaryMapper.h"

#include "vtkIntArray.h"
#include "vtkDoubleArray.h"

class vtkSuperSquareBoundaryMapper : public vtkBoundaryMapper
{
public:
  static vtkSuperSquareBoundaryMapper* New();
  vtkTypeMacro(vtkSuperSquareBoundaryMapper,vtkBoundaryMapper);
  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

  // Vector describing how many divisions on each of the four boundaries
  // The total number of divisions should be equal to the total number of
  // boundary ids minus four (the four corners)
  vtkSetVector4Macro(SuperBoundaryDivisions, int);
  vtkGetVector4Macro(SuperBoundaryDivisions, int);

  // The length of each boundary. Default is one for each side
  vtkSetVector4Macro(SuperBoundaryLengths, double);
  vtkGetVector4Macro(SuperBoundaryLengths, double);
protected:
  vtkSuperSquareBoundaryMapper();
  ~vtkSuperSquareBoundaryMapper();

  int SetBoundaries() VTK_OVERRIDE;
  int CalculateSquareEdgeLengths(vtkIntArray *actualIds);
  int SetSquareBoundary(vtkIntArray *actualIds);

  vtkDoubleArray *BoundaryLengths;
  int SuperBoundaryDivisions[4];
  double SuperBoundaryLengths[4];

private:
  vtkSuperSquareBoundaryMapper(const vtkSuperSquareBoundaryMapper&) VTK_DELETE_FUNCTION;
  void operator=(const vtkSuperSquareBoundaryMapper&) VTK_DELETE_FUNCTION;
};

#endif
