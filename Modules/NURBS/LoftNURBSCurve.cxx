//
//  TestControlGrid.cxx
//
//
//  Created by Adam Updegrove on 9/14/16.
//
//

/*=========================================================================

 Program:   SimVascular

 =========================================================================*/

#include "vtkDataObject.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"
#include "vtkSVIOUtils.h"
#include "vtkSVLoftNURBSCurve.h"
#include "vtkSVNURBSUtils.h"

#include <string>
#include <sstream>
#include <iostream>

#include <unistd.h>

int main(int argc, char *argv[])
{
  if (argc != 1)
  {
      std::cout << "Incorrect Usage! Should be:" <<endl;
      std::cout << "./TestNURBSCurve" <<endl;
      return EXIT_FAILURE;
  }

  int nc = 11;

  vtkSmartPointer<vtkDoubleArray> x_data =
    vtkSmartPointer<vtkDoubleArray>::New();
  x_data->SetNumberOfTuples(nc);
  vtkSVNURBSUtils::LinSpace(0, 2*M_PI, nc, x_data);

  vtkSmartPointer<vtkPoints> inputPoints =
    vtkSmartPointer<vtkPoints>::New();
  inputPoints->Reset();
  for (int i=0; i<nc; i++)
  {
    double xval = x_data->GetTuple1(i);
    double yval = sin(xval - M_PI/2.0);
    double zval = 0.0;

    inputPoints->InsertNextPoint(xval, yval, zval);
  }

  vtkSmartPointer<vtkPolyData> inputPoly =
    vtkSmartPointer<vtkPolyData>::New();
  inputPoly->SetPoints(inputPoints);
  vtkSmartPointer<vtkSVLoftNURBSCurve> lofter =
    vtkSmartPointer<vtkSVLoftNURBSCurve>::New();
  lofter->SetInputData(inputPoly);
  lofter->SetDegree(3);
  lofter->SetPolyDataSpacing(0.01);
  lofter->SetKnotSpanType("derivative");
  lofter->SetParametricSpanType("chord");
  lofter->Update();

  std::string newDirName = getcwd(NULL, 0);
  system(("mkdir -p "+newDirName+"/../Tests").c_str());
  //Write Files
  std::cout<<"Done...Writing Files..."<<endl;
  vtkSVIOUtils::WriteVTPFile(newDirName+"/../Tests/Loft_Curve.vtp",lofter->GetOutput(0),"");
  vtkSVIOUtils::WriteVTSFile(newDirName+"/../Tests/Loft_CurveControlPoints.vts",lofter->GetCurve()->GetControlPointGrid(),"");
  std::cout<<"Done"<<endl;

  return EXIT_SUCCESS;
}
