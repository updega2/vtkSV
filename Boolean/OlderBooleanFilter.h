/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkBooleanOperationPolyDataFilter2.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkBooleanOperationPolyDataFilter2
// .SECTION Description
//
// Computes the boundary of the union, intersection, or difference
// volume computed from the volumes defined by two input surfaces. The
// two surfaces do not need to be manifold, but if they are not,
// unexpected results may be obtained. The resulting surface is
// available in the first output of the filter. The second output
// contains a set of polylines that represent the intersection between
// the two input surfaces.
//
// This code was contributed in the VTK Journal paper:
// "Boolean Operations on Surfaces in VTK Without External Libraries"
// by Cory Quammen, Chris Weigle C., Russ Taylor
// http://hdl.handle.net/10380/3262
// http://www.midasjournal.org/browse/publication/797

#ifndef __vtkBooleanOperationPolyDataFilter2_h
#define __vtkBooleanOperationPolyDataFilter2_h

#include "vtkFiltersGeneralModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"

#include "vtkDataSetAttributes.h" // Needed for CopyCells() method

class vtkIdList;

class VTKFILTERSGENERAL_EXPORT vtkBooleanOperationPolyDataFilter2 : public vtkPolyDataAlgorithm
{
public:
  // Description:
  // Construct object that computes the boolean surface.
  static vtkBooleanOperationPolyDataFilter2 *New();

  vtkTypeMacro(vtkBooleanOperationPolyDataFilter2,
               vtkPolyDataAlgorithm);

  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Integer describing the number of intersection points and lines
  vtkGetMacro(NumberOfIntersectionPoints, int);
  vtkSetMacro(NumberOfIntersectionPoints, int);
  vtkBooleanMacro(NumberOfIntersectionPoints, int);

  vtkGetMacro(NumberOfIntersectionLines, int);
  vtkSetMacro(NumberOfIntersectionLines, int);
  vtkBooleanMacro(NumberOfIntersectionLines, int);

  // Description:
  // Variable to determine what is output if no intersection occurs.
  // 0 = neither (default), 1 = first, 2 = second, 3 = both
  vtkGetMacro(NoIntersectionOutput, int);
  vtkSetMacro(NoIntersectionOutput, int);
  vtkBooleanMacro(NoIntersectionOutput, int);
  enum OperationType
  {
    VTK_UNION=0,
    VTK_INTERSECTION,
    VTK_DIFFERENCE
  };
  enum NoIntersectionOutputType
  {
    VTK_FIRST=0,
    VTK_SECOND,
    VTK_BOTH,
    VTK_NEITHER
  };

  // Description:
  // Set the boolean operation to perform. Defaults to union.
  vtkSetClampMacro( Operation, int, VTK_UNION, VTK_DIFFERENCE );
  vtkGetMacro( Operation, int );
  void SetOperationToUnion()
  { this->SetOperation( VTK_UNION ); }
  void SetOperationToIntersection()
  { this->SetOperation( VTK_INTERSECTION ); }
  void SetOperationToDifference()
  { this->SetOperation( VTK_DIFFERENCE ); }

  // Description:
  // Variable to determine what is output if no intersection occurs.
  // 0 = neither (default), 1 = first, 2 = second, 3 = both
  vtkGetMacro(Tolerance, double);
  vtkSetMacro(Tolerance, double);

protected:
  vtkBooleanOperationPolyDataFilter2();
  ~vtkBooleanOperationPolyDataFilter2();

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
  int FillInputPortInformation(int, vtkInformation*);

private:
  vtkBooleanOperationPolyDataFilter2(const vtkBooleanOperationPolyDataFilter2&); // Not implemented
  void operator=(const vtkBooleanOperationPolyDataFilter2&); // Not implemented

  // Description:
  // PolyDatas for each surface out of intersection and also the intersection
  // lines
  vtkPolyData *OutputSurface;
  // Description:
  // Which operation to perform.
  // Can be VTK_UNION, VTK_INTERSECTION, or VTK_DIFFERENCE.
  int Operation;
  int NoIntersectionOutput;
  int NumberOfIntersectionPoints;
  int NumberOfIntersectionLines;
  double Tolerance;

  class Impl;

};

#endif
