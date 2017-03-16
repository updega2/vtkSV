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

/**
 *  \file TestHausdorffDistance.cxx
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */
#include <vtkSVHausdorffDistance.h>

#include <vtkPlaneSource.h>
#include <vtkSmartPointer.h>
#include <vtkSVGlobals.h>
#include <vtkSVHausdorffDistance.h>

int TestHausdorffDistance(int, char *[])
{
  vtkNew(vtkPlaneSource, plane0);
  plane0->SetOrigin(0.0, 0.0, 0.0);
  plane0->SetPoint1(1.0, 0.0, 0.0);
  plane0->SetPoint2(0.0, 1.0, 0.0);
  plane0->SetXResolution(4);
  plane0->SetYResolution(4);
  plane0->Update();
  vtkNew(vtkPolyData, pd0);
  pd0->ShallowCopy(plane0->GetOutput());

  vtkNew(vtkPlaneSource, plane1);
  plane1->SetOrigin(0.0, 0.0, 1.0);
  plane1->SetPoint1(1.0, 0.0, 1.0);
  plane1->SetPoint2(0.0, 1.0, 1.0);
  plane1->SetXResolution(4);
  plane1->SetYResolution(4);
  plane1->Update();
  vtkNew(vtkPolyData, pd1);
  pd1->ShallowCopy(plane1->GetOutput());

  vtkNew(vtkSVHausdorffDistance, Distancer);
  Distancer->SetInputData(0, pd0);
  Distancer->SetInputData(1, pd1);
  Distancer->SetDistanceArrayName("Distance");
  Distancer->Update();

  std::cout<<"Hausdorff Distance: "<<Distancer->GetHausdorffDistance()<<endl;
  std::cout<<"Average Distance:   "<<Distancer->GetAverageDistance()<<endl;
  if (0.9999 > Distancer->GetHausdorffDistance() ||
      1.0001 < Distancer->GetHausdorffDistance())
  {
    std::cout<<"Incorrect distance calculation"<<endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
