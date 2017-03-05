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

// .NAME vtkPassDataArray - Get Boundary Faces from poldata and label them with integers
// .SECTION Description
// vtkPassDataArray is a filter to extract the boundary surfaces of a model, separate the surace into multiple regions and number each region.

// .SECTION Caveats
// To see the coloring of the lines you may have to set the ScalarMode
// instance variable of the mapper to SetScalarModeToUseCellData(). (This
// is only a problem if there are point data scalars.)

// .SECTION See Also
// vtkExtractEdges

/** @file vtkPassDataArray.h
 *  @brief This filter passes data information from one vtkPolyData to another.
 *  These polydatas do not need to be associated in any way. It uses
 *  vtkPointLocator and vtkCellLocators to find the closest points and pass
 *  the information. It passes the array set with PassArrayName. Will be modified
 *  in the future to pass all data arrays if specified.
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#ifndef __vtkPassDataArray_h
#define __vtkPassDataArray_h

#include "vtkPolyDataAlgorithm.h"

class vtkPassDataArray : public vtkPolyDataAlgorithm
{
public:
  static vtkPassDataArray* New();
  //vtkTypeRevisionMacro(vtkPassDataArray, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set name for data array to be used to determine the in between sections
  vtkGetStringMacro(PassArrayName);
  vtkSetStringMacro(PassArrayName);

  // Description:
  // Set/get macros for telling whether source is cell or point data
  // and how the information should be transferred
  vtkGetMacro(PassDataIsCellData, int);
  vtkSetMacro(PassDataIsCellData, int);
  vtkGetMacro(PassDataToCellData, int);
  vtkSetMacro(PassDataToCellData, int);


  // Description:
  // Set/get macros for whether to use cell centroid. Default is true
  // As false, each point of the cell will be tested and the value returned
  // the most times from the locator will be used. In the case that multiple
  // values are returned the same amount, the first is used.
  vtkGetMacro(UseCellCentroid, int);
  vtkSetMacro(UseCellCentroid, int);

protected:
  vtkPassDataArray();
  ~vtkPassDataArray();

  // Usual data generation method
  int RequestData(vtkInformation *vtkNotUsed(request),
		  vtkInformationVector **inputVector,
		  vtkInformationVector *outputVector);

  char* PassArrayName;

  vtkDataArray *PassDataArray;
  vtkDataArray *NewDataArray;

  vtkPolyData *SourcePd;
  vtkPolyData *TargetPd;

  int PassDataIsCellData;
  int PassDataToCellData;
  int UseCellCentroid;

  int GetArrays(vtkPolyData *object,int type);
  void GetMostOccuringId(vtkIdList *idList, vtkIdType &output);

  int PassDataInformation();
  int PassInformationToPoints(vtkPolyData *sourcePd, vtkPolyData *targetPd,
                              const int sourceIsCellData, vtkDataArray *sourceDataArray,
                              vtkDataArray *targetDataArray);

  int PassInformationToCells(vtkPolyData *sourcePd, vtkPolyData *targetPd,
                             const int sourceIsCellData, const int useCellCentroid,
                             vtkDataArray *sourceDataArray,
                             vtkDataArray *targetDataArray);

private:
  vtkPassDataArray(const vtkPassDataArray&);  // Not implemented.
  void operator=(const vtkPassDataArray&);  // Not implemented.
};

#endif


