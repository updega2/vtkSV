/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSquareBoundaryMapper.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkSquareBoundaryMapper
 * @brief   base class for boundary mapping filters
 *
 * vtkSquareBoundaryMapper is an abstract class that defines
 * the protocol for suboundary mapping surface filters.
 *
*/

#ifndef vtkSquareBoundaryMapper_h
#define vtkSquareBoundaryMapper_h

#include "vtkFiltersGeneralModule.h" // For export macro
#include "vtkBoundaryMapper.h"

#include "vtkIntArray.h"

class vtkSquareBoundaryMapper : public vtkBoundaryMapper
{
public:
  static vtkSquareBoundaryMapper* New();
  vtkTypeMacro(vtkSquareBoundaryMapper,vtkBoundaryMapper);
  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

protected:
  vtkSquareBoundaryMapper() {}

  int SetBoundaries() VTK_OVERRIDE;
  int CalculateSquareEdgeLengths();
  int SetSquareBoundary();

  double BoundaryLengths[4];

private:
  vtkSquareBoundaryMapper(const vtkSquareBoundaryMapper&) VTK_DELETE_FUNCTION;
  void operator=(const vtkSquareBoundaryMapper&) VTK_DELETE_FUNCTION;
};

#endif
