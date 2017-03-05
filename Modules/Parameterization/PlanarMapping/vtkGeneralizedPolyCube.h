/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkGeneralizedPolycube.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkGeneralizedPolycube - topologically regular array of data
// .SECTION Description
//
// Inherets from vtkDataObject

#ifndef vtkGeneralizedPolycube_h
#define vtkGeneralizedPolycube_h

#include "vtkUnstructuredGrid.h"

#include "vtkStructuredGrid.h"
#include "vtkDenseArray.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkPolyData.h"

class vtkGeneralizedPolycube : public vtkUnstructuredGrid
{
public:
  static vtkGeneralizedPolycube *New();
  vtkGeneralizedPolycube(int m, vtkPoints *controlPoints, int n, vtkDoubleArray *knotPoints, int deg) {;}
  vtkGeneralizedPolycube(int m, vtkPoints *controlPoints, vtkDoubleArray *knotPoints, vtkIntArray *knotMultiplicity, int deg) {;}

  vtkTypeMacro(vtkGeneralizedPolycube,vtkUnstructuredGrid);
  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

  // Description: Surgerylines that define a more exact path in between individual starting points of the polycube. Do not have to be used
  vtkGetObjectMacro(SurgeryLines, vtkPolyData);
  vtkSetObjectMacro(SurgeryLines, vtkPolyData);

  //CUBE_TYPE
  enum CUBE_TYPE
  {
    CUBE_BRANCH = 0,
    CUBE_BIFURCATION
  };

  void Initialize() VTK_OVERRIDE;

  // Origin is front, left, bottom corner of cube when axis aligned
  int InsertGridWithOrigin(const int cellId, const double origin[3], const double dims[3],
                        const int cubetype, const int parentdirection, const int childdirection);
  int SetGridWithOrigin(const int cellId, const double origin[3], const double dims[3],
                        const int cubetype, const int parentdirection, const int childdirection);
  int SetGridWithOrigin(const int cellId, const double origin[3], const double dims[3],
                        const int cubetype, const int parentdirection, const int childdirection, const double topNormal[3], const double rightNormal[3], const int corners[4]);
  int InsertGridWithCenter(const int cellId, const double center[3], const double dims[3],
                        const int cubetype, const int parentdirection, const int childdirection);
  int SetGridWithCenter(const int cellId, const double center[3], const double dims[3],
                        const int cubetype, const int parentdirection, const int childdirection);
  int InsertGrid(const int cellId, vtkPoints *points, const int cubetype, const int parentdirection, const int childdirection);
  int SetGrid(const int cellId, vtkPoints *points, const int cubetype, const int parentdirection, const int childdirection);
  int GetFullRepresentation(vtkUnstructuredGrid *fullRepresentation) {return 0;}
  int GetGrid(const int cellId, const int spacing, vtkStructuredGrid *gridRepresentation) {return 0;}
  void SetNumberOfGrids(const int numberOfGrids);
  int GetNumberOfGrids();

protected:
  vtkGeneralizedPolycube();
  ~vtkGeneralizedPolycube();

  vtkPolyData    *SurgeryLines;

private:
  vtkGeneralizedPolycube(const vtkGeneralizedPolycube&);  // Not implemented.
  void operator=(const vtkGeneralizedPolycube&);  // Not implemented.
};

#endif
