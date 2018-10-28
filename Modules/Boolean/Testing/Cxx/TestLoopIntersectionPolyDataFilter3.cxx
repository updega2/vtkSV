/* Copyright (c) Stanford University, The Regents of the University of
 *               California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
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
 */

#include <vtkSVLoopIntersectionPolyDataFilter.h>

#include <vtkActor.h>
#include <vtkConeSource.h>
#include <vtkCubeSource.h>
#include <vtkLinearSubdivisionFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTriangleFilter.h>

int TestLoopIntersectionPolyDataFilter3(int argc, char *argv[])
{
  vtkSmartPointer<vtkCubeSource> cubeSource =
    vtkSmartPointer<vtkCubeSource>::New();
  cubeSource->SetCenter(0.0, 0.0, 0.0);
  cubeSource->SetXLength(1.0);
  cubeSource->SetYLength(1.0);
  cubeSource->SetZLength(1.0);
  cubeSource->Update();
  vtkSmartPointer<vtkTriangleFilter> cubetriangulator =
    vtkSmartPointer<vtkTriangleFilter>::New();
  cubetriangulator->SetInputConnection(cubeSource->GetOutputPort());
  vtkSmartPointer<vtkLinearSubdivisionFilter> cubesubdivider =
    vtkSmartPointer<vtkLinearSubdivisionFilter>::New();
  cubesubdivider->SetInputConnection(cubetriangulator->GetOutputPort());
  cubesubdivider->SetNumberOfSubdivisions(3);
  vtkSmartPointer<vtkPolyDataMapper> cubeMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  cubeMapper->SetInputConnection( cubesubdivider->GetOutputPort() );
  cubeMapper->ScalarVisibilityOff();
  vtkSmartPointer<vtkActor> cubeActor =
    vtkSmartPointer<vtkActor>::New();
  cubeActor->SetMapper( cubeMapper );
  cubeActor->GetProperty()->SetOpacity(.3);
  cubeActor->GetProperty()->SetColor(1,0,0);
  cubeActor->GetProperty()->SetInterpolationToFlat();

  vtkSmartPointer<vtkConeSource> coneSource =
    vtkSmartPointer<vtkConeSource>::New();
  coneSource->SetCenter(0.0, 0.0, 0.0);
  coneSource->SetRadius(0.5);
  coneSource->SetHeight(2.0);
  coneSource->SetResolution(10.0);
  coneSource->SetDirection(1,0,0);
  vtkSmartPointer<vtkTriangleFilter> conetriangulator =
    vtkSmartPointer<vtkTriangleFilter>::New();
  conetriangulator->SetInputConnection(coneSource->GetOutputPort());
  vtkSmartPointer<vtkLinearSubdivisionFilter> conesubdivider =
    vtkSmartPointer<vtkLinearSubdivisionFilter>::New();
  conesubdivider->SetInputConnection(conetriangulator->GetOutputPort());
  conesubdivider->SetNumberOfSubdivisions(3);
  vtkSmartPointer<vtkPolyDataMapper> coneMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  coneMapper->SetInputConnection( conesubdivider->GetOutputPort() );
  coneMapper->ScalarVisibilityOff();
  vtkSmartPointer<vtkActor> coneActor =
    vtkSmartPointer<vtkActor>::New();
  coneActor->SetMapper( coneMapper );
  coneActor->GetProperty()->SetOpacity(.3);
  coneActor->GetProperty()->SetColor(0,1,0);
  coneActor->GetProperty()->SetInterpolationToFlat();

  vtkSmartPointer<vtkSVLoopIntersectionPolyDataFilter> intersectionPolyDataFilter =
    vtkSmartPointer<vtkSVLoopIntersectionPolyDataFilter>::New();
  intersectionPolyDataFilter->SetInputConnection( 0, cubesubdivider->GetOutputPort() );
  intersectionPolyDataFilter->SetInputConnection( 1, conesubdivider->GetOutputPort() );
  intersectionPolyDataFilter->SetSplitFirstOutput(0);
  intersectionPolyDataFilter->SetSplitSecondOutput(0);
  intersectionPolyDataFilter->Update();

  vtkSmartPointer<vtkPolyDataMapper> intersectionMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  intersectionMapper->SetInputConnection( intersectionPolyDataFilter->GetOutputPort() );
  intersectionMapper->ScalarVisibilityOff();

  vtkSmartPointer<vtkActor> intersectionActor =
    vtkSmartPointer<vtkActor>::New();
  intersectionActor->SetMapper( intersectionMapper );

  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  renderer->AddViewProp(cubeActor);
  renderer->AddViewProp(coneActor);
  renderer->AddViewProp(intersectionActor);
  renderer->SetBackground(.1, .2, .3);

  vtkSmartPointer<vtkRenderWindow> renderWindow
    = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer( renderer );

  vtkSmartPointer<vtkRenderWindowInteractor> renWinInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renWinInteractor->SetRenderWindow( renderWindow );

  intersectionPolyDataFilter->Print(std::cout);

  renderWindow->Render();
  renWinInteractor->Start();

  return EXIT_SUCCESS;
}
