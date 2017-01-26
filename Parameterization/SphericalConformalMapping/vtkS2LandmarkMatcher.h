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


/** @file vtkS2LandmarkMatcher.h
 *  @brief This is a vtk filter to map a triangulated surface to a sphere.
 *  @details This filter uses the heat flow method to map a triangulated
 *  surface to a sphere. The first step is to compute the Landmark Energy, and
 *  the second step is to perform the conformal map. For more details, see
 *  Gu et al., Genus Zero Surface Conformal Mapping and Its
 *  Application to Brain Surface Mapping, 2004.
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#ifndef __vtkS2LandmarkMatcher_h
#define __vtkS2LandmarkMatcher_h

#include "vtkPolyDataAlgorithm.h"

#include "vtkEdgeTable.h"
#include "vtkFloatArray.h"
#include "vtkPolyData.h"

#include <complex>
#include <vector>

class vtkS2LandmarkMatcher : public vtkPolyDataAlgorithm
{
public:
  static vtkS2LandmarkMatcher* New();
  //vtkTypeRevisionMacro(vtkS2LandmarkMatcher, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Print statements used for debugging
  vtkGetMacro(Verbose, int);
  vtkSetMacro(Verbose, int);

  // Description:
  // Get and set macros for the time step to be taken during the computation
  // of the Landmark Energy and the Harmonic Energy
  vtkGetMacro(LandmarkTimeStep, double);
  vtkSetMacro(LandmarkTimeStep, double);

  // Description:
  // Get and set macros for the Energy criterion during the computation of
  // the Landmark Energy and the Harmonic Energy
  vtkGetMacro(LandmarkEnergyCriterion, double);
  vtkSetMacro(LandmarkEnergyCriterion, double);

  // Description:
  // Get and set macros for the Landmark lagrange multiplier factor. The higher
  // the weighting, the more it attempts to exactly match landmark points.
  vtkGetMacro(LandmarkWeighting, double);
  vtkSetMacro(LandmarkWeighting, double);

  // Description:
  // Print statements used for number of iterations
  vtkGetMacro(MaxNumIterations, int);
  vtkSetMacro(MaxNumIterations, int);

  // CG Update method
  vtkGetMacro(CGUpdateMethod, int);
  vtkSetMacro(CGUpdateMethod, int);

  // Description:
  // Landmark pt ids and desired location
  vtkSetObjectMacro(SourceLandmarkPtIds, vtkIntArray);
  vtkSetObjectMacro(TargetLandmarkPtIds, vtkIntArray);

  static int ComputeMeshLandmarkLaplacian(vtkPolyData *pd, vtkEdgeTable *edgeTable,
                                          vtkFloatArray *edgeWeights, vtkIntArray *edgeNeighbors,
                                          vtkFloatArray *laplacian,
                                          double landmarkWeighting,
                                          vtkPolyData *targetS2Pd,
                                          vtkIntArray *isCurrentLandmark,
                                          vtkIntArray *landmarkPtIds);
  static int ComputePointLandmarkLaplacian(vtkIdType p0, vtkPolyData *pd,
                                           vtkEdgeTable *edgeTable, vtkFloatArray *edgeWeights,
                                           vtkIntArray *edgeNeighbors, double laplacian[],
                                           double landmarkWeighting,
                                           vtkPolyData *targetS2Pd,
                                           vtkIntArray *isCurrentLandmark,
                                           vtkIntArray *landmarkPtIds);
  static int ComputeLandmarkEnergy(vtkPolyData *pd,
                                   vtkEdgeTable *edgeTable,
                                   vtkFloatArray *edgeWeights,
                                   double &landmarkEnergy,
                                   double landmarkWeighting,
                                   vtkPolyData *targetS2Pd,
                                   vtkIntArray *isCurrentLandmark,
                                   vtkIntArray *landmarkPtIds);
  // Setup and Check Functions
protected:
  vtkS2LandmarkMatcher();
  ~vtkS2LandmarkMatcher();

  // Usual data generation method
  int RequestData(vtkInformation *vtkNotUsed(request),
		  vtkInformationVector **inputVector,
		  vtkInformationVector *outputVector);

  // Main functions in filter
  int RunLandmarkMatching();
  int InitiateCGArrays();
  int LandmarkMapping();
  int SetLandmarks();
  int FirstStep();

  // Vector functions in vtk!
  int WolfeLineSearch();
  int UpdateMap(vtkFloatArray *laplacian, int cg_update);//Sets current descent direction without cg
  int StepForward(int cg_update); // https://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
  int FRUpdateMap(); //Fletcher-Reeves
  int PRUpdateMap(); //Polak-Ribier
  int HSUpdateMap(); //Hesteness-Stiefel
  int DYUpdateMap(); //Dai-Yuan
  int CGUpdateMap(double beta[]);

private:
  vtkS2LandmarkMatcher(const vtkS2LandmarkMatcher&);  // Not implemented.
  void operator=(const vtkS2LandmarkMatcher&);  // Not implemented.

  int    Verbose;
  double LandmarkTimeStep;
  double LandmarkEnergyCriterion;
  double LandmarkWeighting;
  int    MaxNumIterations;
  int    NumLandmarks;
  int    CGUpdateMethod;

  vtkPolyData   *SourcePd;
  vtkPolyData   *TargetPd;
  vtkPolyData   *MappedPd;

  vtkEdgeTable  *SourceEdgeTable;
  vtkFloatArray *SourceEdgeWeights;
  vtkFloatArray *PrevDescent;
  vtkFloatArray *CurrDescent;
  vtkFloatArray *ConjugateDir;
  vtkIntArray   *SourceEdgeNeighbors;
  vtkIntArray   *SourceIsBoundary;

  vtkIntArray *SourceLandmarkPtIds;
  vtkIntArray *IsSourceLandmark;
  vtkIntArray *TargetLandmarkPtIds;
  vtkIntArray *IsTargetLandmark;

};

#endif


