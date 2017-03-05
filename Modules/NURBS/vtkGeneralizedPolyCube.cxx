/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkGeneralizedPolyCube.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkGeneralizedPolyCube.h"

#define vtkNew(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()


vtkStandardNewMacro(vtkGeneralizedPolyCube);

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
vtkGeneralizedPolyCube::vtkGeneralizedPolyCube()
{
  this->FullRepresentation = vtkUnstructuredGrid::New();

  this->Centerlines = vtkPolyData::New();
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
vtkGeneralizedPolyCube::~vtkGeneralizedPolyCube()
{
  if (this->FullRepresentation != NULL)
  {
    this->FullRepresentation->Delete();
  }
  if (this->Centerlines != NULL)
  {
    this->Centerlines->Delete();
  }
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkGeneralizedPolyCube::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkGeneralizedPolyCube::Initialize()
{
  this->Superclass::Initialize();
}
