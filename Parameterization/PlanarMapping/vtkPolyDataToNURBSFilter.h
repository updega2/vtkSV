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


/** @file vtkPolyDataToNURBSFilter.h
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

#ifndef vtkPolyDataToNURBSFilter_h
#define vtkPolyDataToNURBSFilter_h

#include "vtkPolyDataAlgorithm.h"

#include "vtkAppendPolyData.h"
#include "vtkEdgeTable.h"
#include "vtkFloatArray.h"
#include "vtkGeneralizedPolycube.h"
#include "vtkPolyData.h"

class vtkPolyDataToNURBSFilter : public vtkPolyDataAlgorithm
{
public:
  static vtkPolyDataToNURBSFilter* New();
  //vtkTypeRevisionMacro(vtkPolyDataToNURBSFilter, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Print statements used for debugging
  vtkGetMacro(Verbose, int);
  vtkSetMacro(Verbose, int);

  // Description:
  // Array names for the segment and group ids. Must exist on input
  // pd if also providing centerlines. Otherwise, these are the names
  // that will be used to segment the polydata with the created centerlines
  vtkGetStringMacro(SliceIdsArrayName);
  vtkSetStringMacro(SliceIdsArrayName);
  vtkGetStringMacro(GroupIdsArrayName);
  vtkSetStringMacro(GroupIdsArrayName);
  vtkGetStringMacro(SegmentIdsArrayName);
  vtkSetStringMacro(SegmentIdsArrayName);
  vtkGetStringMacro(SphereRadiusArrayName);
  vtkSetStringMacro(SphereRadiusArrayName);
  vtkGetStringMacro(BoundaryPointsArrayName);
  vtkSetStringMacro(BoundaryPointsArrayName);
  vtkGetStringMacro(InternalIdsArrayName);
  vtkSetStringMacro(InternalIdsArrayName);
  vtkGetStringMacro(DijkstraArrayName);
  vtkSetStringMacro(DijkstraArrayName);
  vtkGetStringMacro(BooleanPathArrayName);
  vtkSetStringMacro(BooleanPathArrayName);

  // Description:
  // Macro to set object centerlines
  vtkGetObjectMacro(Centerlines, vtkPolyData);
  vtkSetObjectMacro(Centerlines, vtkPolyData);
  vtkGetObjectMacro(CubeS2Pd, vtkPolyData);
  vtkSetObjectMacro(CubeS2Pd, vtkPolyData);
  vtkGetObjectMacro(OpenCubeS2Pd, vtkPolyData);
  vtkSetObjectMacro(OpenCubeS2Pd, vtkPolyData);

protected:
  vtkPolyDataToNURBSFilter();
  ~vtkPolyDataToNURBSFilter();

  // Usual data generation method
  int RequestData(vtkInformation *vtkNotUsed(request),
		  vtkInformationVector **inputVector,
		  vtkInformationVector *outputVector);

  int PerformMappings();
  int ComputeCenterlines();
  int ExtractBranches();
  int SliceAndDice();
  int GetSegment(const int segmentId, vtkPolyData *segmentPd, vtkPolyData *surgeryLinePd);
  int GetSlice(const int sliceId, vtkPolyData *segmentPd, vtkPolyData *slicePd);
  int MapBranch(const int branchId, vtkAppendPolyData *appender);
  int MapBifurcation(const int bifurcationId, vtkAppendPolyData *appender);
  int MapSliceToS2(vtkPolyData *slicePd, vtkPolyData *surgeryLinePd, vtkPolyData *sliceS2Pd,
                     vtkIntArray *firstCorners, vtkIntArray *secondCorners,
                     double xvec[3], double zvec[3]);
  int MapOpenSliceToS2(vtkPolyData *slicePd, vtkPolyData *sliceS2Pd,
                       vtkIntArray *firstCorners, vtkIntArray *secondCorners,
                       double xvec[3], double zvec[3]);
  int InterpolateMapOntoTarget(vtkPolyData *sourceS2Pd,
                               vtkPolyData *targetPd,
                               vtkPolyData *targetS2Pd,
                               vtkPolyData *mappedPd);

private:
  vtkPolyDataToNURBSFilter(const vtkPolyDataToNURBSFilter&);  // Not implemented.
  void operator=(const vtkPolyDataToNURBSFilter&);  // Not implemented.

  int Verbose;

  vtkPolyData *InputPd;
  vtkPolyData *CubeS2Pd;
  vtkPolyData *OpenCubeS2Pd;
  vtkPolyData *ParameterizedPd;
  vtkPolyData *Centerlines;
  vtkPolyData *SurgeryLines;

  char *SliceIdsArrayName;
  char *GroupIdsArrayName;
  char *SegmentIdsArrayName;
  char *SphereRadiusArrayName;
  char *InternalIdsArrayName;
  char *BoundaryPointsArrayName;
  char *DijkstraArrayName;
  char *BooleanPathArrayName;

  vtkGeneralizedPolycube *Polycube;

};

#endif
