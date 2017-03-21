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

#include <vtkSVLoopIntersectionPolyDataFilter.h>

#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>

int TestLoopIntersectionPolyDataFilter(int argc, char *argv[])
{
  vtkSmartPointer<vtkSphereSource> sphereSource1 =
    vtkSmartPointer<vtkSphereSource>::New();
  sphereSource1->SetCenter(0.0, 0.0, 0.0);
  sphereSource1->SetRadius(2.0);
  sphereSource1->SetPhiResolution(11);
  sphereSource1->SetThetaResolution(21);
  sphereSource1->Update();
  vtkSmartPointer<vtkPolyDataMapper> sphere1Mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  sphere1Mapper->SetInputConnection( sphereSource1->GetOutputPort() );
  sphere1Mapper->ScalarVisibilityOff();
  vtkSmartPointer<vtkActor> sphere1Actor =
    vtkSmartPointer<vtkActor>::New();
  sphere1Actor->SetMapper( sphere1Mapper );
  sphere1Actor->GetProperty()->SetOpacity(.3);
  sphere1Actor->GetProperty()->SetColor(1,0,0);
  sphere1Actor->GetProperty()->SetInterpolationToFlat();

  vtkSmartPointer<vtkSphereSource> sphereSource2 =
    vtkSmartPointer<vtkSphereSource>::New();
  sphereSource2->SetCenter(1.0, 0.0, 0.0);
  sphereSource2->SetRadius(2.0);
  vtkSmartPointer<vtkPolyDataMapper> sphere2Mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  sphere2Mapper->SetInputConnection( sphereSource2->GetOutputPort() );
  sphere2Mapper->ScalarVisibilityOff();
  vtkSmartPointer<vtkActor> sphere2Actor =
    vtkSmartPointer<vtkActor>::New();
  sphere2Actor->SetMapper( sphere2Mapper );
  sphere2Actor->GetProperty()->SetOpacity(.3);
  sphere2Actor->GetProperty()->SetColor(0,1,0);
  sphere2Actor->GetProperty()->SetInterpolationToFlat();

  vtkSmartPointer<vtkSVLoopIntersectionPolyDataFilter> intersectionPolyDataFilter =
    vtkSmartPointer<vtkSVLoopIntersectionPolyDataFilter>::New();
  intersectionPolyDataFilter->SetInputConnection( 0, sphereSource1->GetOutputPort() );
  intersectionPolyDataFilter->SetInputConnection( 1, sphereSource2->GetOutputPort() );
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
  renderer->AddViewProp(sphere1Actor);
  renderer->AddViewProp(sphere2Actor);
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
