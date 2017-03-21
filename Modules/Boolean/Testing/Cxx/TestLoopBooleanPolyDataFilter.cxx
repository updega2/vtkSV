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

#include "vtkActor.h"
#include "vtkSVLoopBooleanPolyDataFilter.h"
#include "vtkIntersectionPolyDataFilter.h"
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkCubeSource.h>
#include <vtkCylinderSource.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkLinearSubdivisionFilter.h>
#include <vtkTransform.h>
#include <vtkTriangleFilter.h>
#include <vtkCellData.h>
#include <vtkMath.h>

static vtkActor* GetCubeBooleanOperationActor( double x, int operation )
{
  vtkSmartPointer<vtkCubeSource> cube1 =
    vtkSmartPointer<vtkCubeSource>::New();
  cube1->SetCenter(x, 4.0, 0.0);
  cube1->SetXLength(1.0);
  cube1->SetYLength(1.0);
  cube1->SetZLength(1.0);
  cube1->Update();

  vtkSmartPointer<vtkTriangleFilter> triangulator1 =
    vtkSmartPointer<vtkTriangleFilter>::New();
  triangulator1->SetInputData(cube1->GetOutput());
  triangulator1->Update();

  vtkSmartPointer<vtkLinearSubdivisionFilter> subdivider1 =
    vtkSmartPointer<vtkLinearSubdivisionFilter>::New();
  subdivider1->SetInputData(triangulator1->GetOutput());
  subdivider1->Update();

  vtkSmartPointer<vtkCubeSource> cube2 =
    vtkSmartPointer<vtkCubeSource>::New();
  cube2->SetCenter(x + 0.3, 4.3 , 0.3);
  cube2->SetXLength(1.0);
  cube2->SetYLength(1.0);
  cube2->SetZLength(1.0);
  cube2->Update();

  vtkSmartPointer<vtkTriangleFilter> triangulator2 =
    vtkSmartPointer<vtkTriangleFilter>::New();
  triangulator2->SetInputData(cube2->GetOutput());
  triangulator2->Update();

  vtkSmartPointer<vtkLinearSubdivisionFilter> subdivider2 =
    vtkSmartPointer<vtkLinearSubdivisionFilter>::New();
  subdivider2->SetInputData(triangulator2->GetOutput());
  subdivider2->Update();

  vtkSmartPointer<vtkSVLoopBooleanPolyDataFilter> boolFilter =
    vtkSmartPointer<vtkSVLoopBooleanPolyDataFilter>::New();
  boolFilter->SetOperation( operation );
  boolFilter->SetInputConnection( 0, subdivider1->GetOutputPort() );
  boolFilter->SetInputConnection( 1, subdivider2->GetOutputPort() );
  boolFilter->Update();

  vtkSmartPointer<vtkPolyData> output =
    vtkSmartPointer<vtkPolyData>::New();
  output = boolFilter->GetOutput();
  output->GetCellData()->SetActiveScalars("FreeEdge");
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData( output );
  mapper->SetScalarRange(0,1);
  mapper->SetScalarModeToUseCellData();
  mapper->ScalarVisibilityOn();

  vtkActor *actor = vtkActor::New();
  actor->SetMapper( mapper );

  return actor;
}

static vtkActor* GetSphereBooleanOperationActor( double x, int operation )
{
  double centerSeparation = 0.15;
  vtkSmartPointer<vtkSphereSource> sphere1 =
    vtkSmartPointer<vtkSphereSource>::New();
  sphere1->SetCenter(-centerSeparation + x, 0.0, 0.0);

  vtkSmartPointer<vtkSphereSource> sphere2 =
    vtkSmartPointer<vtkSphereSource>::New();
  sphere2->SetCenter(  centerSeparation + x, 0.0, 0.0);

  vtkSmartPointer<vtkSVLoopBooleanPolyDataFilter> boolFilter =
    vtkSmartPointer<vtkSVLoopBooleanPolyDataFilter>::New();
  boolFilter->SetOperation( operation );
  boolFilter->SetInputConnection( 0, sphere1->GetOutputPort() );
  boolFilter->SetInputConnection( 1, sphere2->GetOutputPort() );
  boolFilter->Update();

  vtkSmartPointer<vtkPolyData> output =
    vtkSmartPointer<vtkPolyData>::New();
  output = boolFilter->GetOutput();
  output->GetCellData()->SetActiveScalars("FreeEdge");
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData( output );
  mapper->SetScalarRange(0,1);
  mapper->SetScalarModeToUseCellData();
  mapper->ScalarVisibilityOn();

  vtkActor *actor = vtkActor::New();
  actor->SetMapper( mapper );

  return actor;
}

static vtkActor* GetCylinderBooleanOperationActor( double x, int operation )
{
  double axis[3]; axis[0] = 0.0; axis[1] = 1.0; axis[2] = 0.0;
  double vec[3]; vec[0] = 0.0; vec[1] = 1.0; vec[2] = 0.0;
  double rotateaxis[3]; vtkMath::Cross(axis,vec,rotateaxis);
  vtkSmartPointer<vtkCylinderSource> cylinder1 =
    vtkSmartPointer<vtkCylinderSource>::New();
  cylinder1->SetCenter(0.0,0.0,0.0);
  cylinder1->SetHeight(2.0);
  cylinder1->SetRadius(0.5);
  cylinder1->SetResolution(15);
  cylinder1->Update();
  double radangle = vtkMath::AngleBetweenVectors(axis,vec);
  double degangle = vtkMath::DegreesFromRadians(radangle);
  vtkSmartPointer<vtkTransform> rotator1 =
    vtkSmartPointer<vtkTransform>::New();
  rotator1->RotateWXYZ(degangle,rotateaxis);

  vtkSmartPointer<vtkTransformPolyDataFilter> polyDataRotator1 =
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  polyDataRotator1->SetInputData(cylinder1->GetOutput());
  polyDataRotator1->SetTransform(rotator1);
  polyDataRotator1->Update();

  vtkSmartPointer<vtkTransform> mover1 =
    vtkSmartPointer<vtkTransform>::New();
  mover1->Translate(x,-4.0,0.0);

  vtkSmartPointer<vtkTransformPolyDataFilter> polyDataMover1 =
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  polyDataMover1->SetInputData(polyDataRotator1->GetOutput());
  polyDataMover1->SetTransform(mover1);
  polyDataMover1->Update();

  vtkSmartPointer<vtkTriangleFilter> triangulator1 =
    vtkSmartPointer<vtkTriangleFilter>::New();
  triangulator1->SetInputData(polyDataMover1->GetOutput());
  triangulator1->Update();

  axis[0] = 1.0; axis[1] = 0.0; axis[2] = 0.0;
  vtkMath::Cross(axis,vec,rotateaxis);
  vtkSmartPointer<vtkCylinderSource> cylinder2 =
    vtkSmartPointer<vtkCylinderSource>::New();
  cylinder2->SetCenter(0.0, 0.0, 0.0);
  cylinder2->SetHeight(2.0);
  cylinder2->SetRadius(0.5);
  cylinder2->SetResolution(15);
  cylinder2->Update();
  radangle = vtkMath::AngleBetweenVectors(axis,vec);
  degangle = vtkMath::DegreesFromRadians(radangle);
  vtkSmartPointer<vtkTransform> rotator2 =
    vtkSmartPointer<vtkTransform>::New();
  rotator2->RotateWXYZ(degangle,rotateaxis);

  vtkSmartPointer<vtkTransformPolyDataFilter> polyDataRotator2 =
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  polyDataRotator2->SetInputData(cylinder2->GetOutput());
  polyDataRotator2->SetTransform(rotator2);
  polyDataRotator2->Update();

  vtkSmartPointer<vtkTransform> mover2 =
    vtkSmartPointer<vtkTransform>::New();
  mover2->Translate(x,-4.0,0.0);

  vtkSmartPointer<vtkTransformPolyDataFilter> polyDataMover2 =
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  polyDataMover2->SetInputData(polyDataRotator2->GetOutput());
  polyDataMover2->SetTransform(mover2);
  polyDataMover2->Update();

  vtkSmartPointer<vtkTriangleFilter> triangulator2 =
    vtkSmartPointer<vtkTriangleFilter>::New();
  triangulator2->SetInputData(polyDataMover2->GetOutput());
  triangulator2->Update();

  vtkSmartPointer<vtkSVLoopBooleanPolyDataFilter> boolFilter =
    vtkSmartPointer<vtkSVLoopBooleanPolyDataFilter>::New();
  boolFilter->SetOperation( operation );
  boolFilter->SetInputConnection( 0, triangulator1->GetOutputPort() );
  boolFilter->SetInputConnection( 1, triangulator2->GetOutputPort() );
  boolFilter->Update();

  vtkSmartPointer<vtkPolyData> output =
    vtkSmartPointer<vtkPolyData>::New();
  output = boolFilter->GetOutput();
  output->GetCellData()->SetActiveScalars("FreeEdge");
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData( output );
  mapper->SetScalarRange(0,1);
  mapper->SetScalarModeToUseCellData();
  mapper->ScalarVisibilityOn();

  vtkActor *actor = vtkActor::New();
  actor->SetMapper( mapper );

  return actor;
}

int TestLoopBooleanPolyDataFilter(int argc, char *argv[])
{
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renWin =
    vtkSmartPointer<vtkRenderWindow>::New();
  renWin->AddRenderer( renderer );
  renderer->SetBackground(.1, .2, .3);

  vtkSmartPointer<vtkRenderWindowInteractor> renWinInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renWinInteractor->SetRenderWindow( renWin );

  //Sphere
  vtkActor *unionActor =
    GetSphereBooleanOperationActor( -2.0, vtkSVLoopBooleanPolyDataFilter::VTK_UNION );
  renderer->AddActor( unionActor );
  unionActor->Delete();

  vtkActor *intersectionActor =
    GetSphereBooleanOperationActor(  0.0, vtkSVLoopBooleanPolyDataFilter::VTK_INTERSECTION );
  renderer->AddActor( intersectionActor );
  intersectionActor->Delete();

  vtkActor *differenceActor =
    GetSphereBooleanOperationActor(  2.0, vtkSVLoopBooleanPolyDataFilter::VTK_DIFFERENCE );
  renderer->AddActor( differenceActor );
  differenceActor->Delete();

  //Cube
  //vtkActor *unionCubeActor =
  //  GetCubeBooleanOperationActor( -2.0, vtkSVLoopBooleanPolyDataFilter::VTK_UNION );
  //renderer->AddActor( unionCubeActor );
  //unionCubeActor->Delete();

  //vtkActor *intersectionCubeActor =
  //  GetCubeBooleanOperationActor(  0.0, vtkSVLoopBooleanPolyDataFilter::VTK_INTERSECTION );
  //renderer->AddActor( intersectionCubeActor );
  //intersectionCubeActor->Delete();

  //vtkActor *differenceCubeActor =
  //  GetCubeBooleanOperationActor(  2.0, vtkSVLoopBooleanPolyDataFilter::VTK_DIFFERENCE );
  //renderer->AddActor( differenceCubeActor );
  //differenceCubeActor->Delete();

  ////Cylinder
  //vtkActor *unionCylinderActor =
  //  GetCylinderBooleanOperationActor( -2.0, vtkSVLoopBooleanPolyDataFilter::VTK_UNION );
  //renderer->AddActor( unionCylinderActor );
  //unionCylinderActor->Delete();

  //vtkActor *intersectionCylinderActor =
  //  GetCylinderBooleanOperationActor(  0.0, vtkSVLoopBooleanPolyDataFilter::VTK_INTERSECTION );
  //renderer->AddActor( intersectionCylinderActor );
  //intersectionCylinderActor->Delete();

  //vtkActor *differenceCylinderActor =
  //  GetCylinderBooleanOperationActor(  2.0, vtkSVLoopBooleanPolyDataFilter::VTK_DIFFERENCE );
  //renderer->AddActor( differenceCylinderActor );
  //differenceCylinderActor->Delete();

  renWin->Render();
  renWinInteractor->Start();

  return EXIT_SUCCESS;
}
