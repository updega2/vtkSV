/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSVSquareBoundaryMapper.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkSVSquareBoundaryMapper
 * @brief   base class for boundary mapping filters
 *
 * vtkSVSquareBoundaryMapper is an abstract class that defines
 * the protocol for suboundary mapping surface filters.
 *
*/

#ifndef vtkSVSquareBoundaryMapper_h
#define vtkSVSquareBoundaryMapper_h

#include "vtkFiltersGeneralModule.h" // For export macro
#include "vtkSVBoundaryMapper.h"

#include "vtkIntArray.h"

class vtkSVSquareBoundaryMapper : public vtkSVBoundaryMapper
{
public:
  static vtkSVSquareBoundaryMapper* New();
  vtkTypeMacro(vtkSVSquareBoundaryMapper,vtkSVBoundaryMapper);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkSVSquareBoundaryMapper() {}

  int SetBoundaries();
  int CalculateSquareEdgeLengths();
  int SetSquareBoundary();

  double BoundaryLengths[4];

private:
  vtkSVSquareBoundaryMapper(const vtkSVSquareBoundaryMapper&);
  void operator=(const vtkSVSquareBoundaryMapper&);
};

#endif
