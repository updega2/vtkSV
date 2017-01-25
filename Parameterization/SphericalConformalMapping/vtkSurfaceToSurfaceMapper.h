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


/** @file vtkSurfaceToSurfaceMapper.h
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

#ifndef __vtkSurfaceToSurfaceMapper_h
#define __vtkSurfaceToSurfaceMapper_h

#include "vtkPolyDataAlgorithm.h"

#include "vtkEdgeTable.h"
#include "vtkFloatArray.h"
#include "vtkPolyData.h"

#include <complex>
#include <vector>

class vtkSurfaceToSurfaceMapper : public vtkPolyDataAlgorithm
{
public:
  static vtkSurfaceToSurfaceMapper* New();
  vtkTypeRevisionMacro(vtkSurfaceToSurfaceMapper, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Print statements used for debugging
  vtkGetMacro(Verbose, int);
  vtkSetMacro(Verbose, int);

  // Description:
  // Get and set macros for the time step to be taken during the computation
  // of the Tutte Energy and the Harmonic Energy
  vtkGetMacro(TutteTimeStep, double);
  vtkSetMacro(TutteTimeStep, double);
  vtkGetMacro(ConformalTimeStep, double);
  vtkSetMacro(ConformalTimeStep, double);

  // Description:
  // Get and set macros for the Energy criterion during the computation of
  // the Tutte Energy and the Harmonic Energy
  vtkGetMacro(TutteEnergyCriterion, double);
  vtkSetMacro(TutteEnergyCriterion, double);
  vtkGetMacro(HarmonicEnergyCriterion, double);
  vtkSetMacro(HarmonicEnergyCriterion, double);

  // Description:
  // Print statements used for number of iterations
  vtkGetMacro(MaxNumIterations, int);
  vtkSetMacro(MaxNumIterations, int);

  // Description:
  // Print statements used for dnumber of subdivisions
  vtkGetMacro(NumSourceSubdivisions, int);
  vtkSetMacro(NumSourceSubdivisions, int);

  // Description:
  // Landmark pt ids and desired location
  vtkSetObjectMacro(SourceLandmarkPtIds, vtkIntArray);
  vtkSetObjectMacro(TargetLandmarkPtIds, vtkIntArray);

  // Functions to set up complex linear system for landmark constraint
  static int ComputeConformalFactor(std::complex<double> z, std::complex<double> &g);
  static int ConvertValueToPolyData(vtkPolyData *inPd, std::string fieldName, vtkPolyData *outPd);


  // Setup and Check Functions
protected:
  vtkSurfaceToSurfaceMapper();
  ~vtkSurfaceToSurfaceMapper();

  // Usual data generation method
  int RequestData(vtkInformation *vtkNotUsed(request),
		  vtkInformationVector **inputVector,
		  vtkInformationVector *outputVector);

  // Main functions in filter
  int RunSphericalConformalMappings();
  int RunLandmarkMatching();
  int SubdivideAndInterpolate();

  int GetMap(vtkPolyData *pd, vtkPolyData *spherePd,
             vtkFloatArray *map);

  // Functions to set landmarks
  int StereographicProjection(vtkPolyData *pd);
  int InverseStereographicProjection(vtkPolyData *pd);
  int RunMobiusTransformation();
  int SetLandmarks();

  // Functions to set up complex linear system for landmark constraint
  int ComputeComplexMobius(std::complex<double> &a, std::complex<double> &b);
  int TransformComplexValue(std::complex<double> a, std::complex<double> b);
  int GetSystem(std::vector<std::complex<double> > &lhs,
                std::vector<std::complex<double> > &rhs);

private:
  vtkSurfaceToSurfaceMapper(const vtkSurfaceToSurfaceMapper&);  // Not implemented.
  void operator=(const vtkSurfaceToSurfaceMapper&);  // Not implemented.

  int    Verbose;
  double TutteTimeStep;
  double ConformalTimeStep;
  double TutteEnergyCriterion;
  double HarmonicEnergyCriterion;
  int    MaxNumIterations;
  int    NumSourceSubdivisions;
  int    NumLandmarks;

  vtkPolyData *SourcePd;
  vtkPolyData *SourceS2Pd;
  vtkPolyData *TargetPd;
  vtkPolyData *TargetS2Pd;
  vtkPolyData *MappedPd;
  vtkPolyData *MappedS2Pd;

  vtkIntArray   *SourceLandmarkPtIds;
  vtkFloatArray *SourceLandmarkProj;
  vtkIntArray   *TargetLandmarkPtIds;
  vtkFloatArray *TargetLandmarkProj;

};

#endif


