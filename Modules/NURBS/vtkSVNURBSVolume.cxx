/*=========================================================================
 *
 * Copyright (c) 2014 The Regents of the University of California.
 * All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *=========================================================================*/

#include "vtkSVNURBSVolume.h"

#include "vtkCellArray.h"
#include "vtkCleanPolyData.h"
#include "vtkSVNURBSUtils.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkSparseArray.h"
#include "vtkSVGlobals.h"

// ----------------------
// StandardNewMacro
// ----------------------
vtkStandardNewMacro(vtkSVNURBSVolume);

// ----------------------
// Constructor
// ----------------------
vtkSVNURBSVolume::vtkSVNURBSVolume()
{
}

// ----------------------
// Destructor
// ----------------------
vtkSVNURBSVolume::~vtkSVNURBSVolume()
{
}

// ----------------------
// PrintSelf
// ----------------------
void vtkSVNURBSVolume::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "Number of control points in u direction: " << this->NumberOfUControlPoints << "\n";
  os << indent << "Number of knot points in u direction: " << this->NumberOfUKnotPoints << "\n";
  os << indent << "U Degree: " << this->UDegree << "\n";
  os << indent << "U Clamped: " << this->UClamped << "\n";
  os << indent << "U Closed: " << this->UClosed << "\n";
  os << "\n";
  os << indent << "Number of control points in v direction: " << this->NumberOfVControlPoints << "\n";
  os << indent << "Number of knot points in v direction: " << this->NumberOfVKnotPoints << "\n";
  os << indent << "V Degree: " << this->VDegree << "\n";
  os << indent << "V Clamped: " << this->VClamped << "\n";
  os << indent << "V Closed: " << this->VClosed << "\n";
  os << indent << "Number of control points in w direction: " << this->NumberOfWControlPoints << "\n";
  os << indent << "Number of knot points in w direction: " << this->NumberOfWKnotPoints << "\n";
  os << indent << "W Degree: " << this->WDegree << "\n";
  os << indent << "W Clamped: " << this->WClamped << "\n";
  os << indent << "W Closed: " << this->WClosed << "\n";
}


// ----------------------
// DeepCopy
// ----------------------
void vtkSVNURBSVolume::DeepCopy(vtkSVNURBSVolume *src)
{
}

// ----------------------
// Iniitialize
// ----------------------
void vtkSVNURBSVolume::Initialize()
{
  this->Superclass::Initialize();
}

// ----------------------
// GetData
// ----------------------
vtkSVNURBSVolume* vtkSVNURBSVolume::GetData(vtkInformation* info)
{
  return info? vtkSVNURBSVolume::SafeDownCast(info->Get(DATA_OBJECT())) : 0;
}

// ----------------------
// GetData
// ----------------------
vtkSVNURBSVolume* vtkSVNURBSVolume::GetData(vtkInformationVector* v, int i)
{
  return vtkSVNURBSVolume::GetData(v->GetInformationObject(i));
}

