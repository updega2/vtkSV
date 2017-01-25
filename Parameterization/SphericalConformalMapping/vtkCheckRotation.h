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


/** @file vtkCheckRotation.h
 *  @brief This is a vtk filter to map a triangulated surface to a sphere.
 *  @details This filter uses the heat flow method to map a triangulated
 *  surface to a sphere. The first step is to compute the Tutte Energy, and
 *  the second step is to perform the conformal map. For more details, see
 *  Gu et al., Genus Zero Surface Conformal Mapping and Its
 *  Application to Brain Surface Mapping, 2004.
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#ifndef vtkCheckRotation_h
#define vtkCheckRotation_h

#include "vtkPolyDataAlgorithm.h"

#include "vtkEdgeTable.h"
#include "vtkFloatArray.h"
#include "vtkPolyData.h"

class vtkCheckRotation : public vtkPolyDataAlgorithm
{
public:
  static vtkCheckRotation* New();
  vtkTypeRevisionMacro(vtkCheckRotation, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Print statements used for debugging
  vtkGetMacro(Verbose, int);
  vtkSetMacro(Verbose, int);

  vtkSetObjectMacro(OriginalPd, vtkPolyData);

  // Setup and Check Functions
  static int GetPolyDataAngles(vtkPolyData *pd, vtkFloatArray *cellAngles);
  static int ComputeMassCenter(vtkPolyData *pd, double massCenter[]);

protected:
  vtkCheckRotation();
  ~vtkCheckRotation();

  // Usual data generation method
  int RequestData(vtkInformation *vtkNotUsed(request),
		  vtkInformationVector **inputVector,
		  vtkInformationVector *outputVector);

  int MoveCenters();
  int FindAndCheckRotation();
  int GetRotationMatrix(double vec0[3], double vec1[3], double rotMatrix[16]);
  int ApplyRotationMatrix(double rotMatrix[16]);
  int CheckAnglesWithOriginal();
  int MatchPointOrder();

private:
  vtkCheckRotation(const vtkCheckRotation&);  // Not implemented.
  void operator=(const vtkCheckRotation&);  // Not implemented.

  int Verbose;
  int CellId;

  vtkPolyData *SourcePd;
  vtkPolyData *TargetPd;
  vtkPolyData *MappedPd;

  vtkPolyData *OriginalPd;
};

#endif
