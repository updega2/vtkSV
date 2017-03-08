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

#include "vtkSVControlGrid.h"
#include "vtkDataObject.h"
#include "vtkSVLoftNURBSSurface.h"
#include "vtkSVNURBSSurface.h"
#include "vtkSVNURBSUtils.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkSTLReader.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLStructuredGridWriter.h"

#include <string>
#include <sstream>
#include <iostream>

#include <unistd.h>

int main(int argc, char *argv[])
{
  if (argc != SV_OK)
  {
      std::cout << "Incorrect Usage! Should be:" <<endl;
      std::cout << "./TestNURBSCurve" <<endl;
      return EXIT_FAILURE;
  }

  int nU = 100;
  int nV = 243;
  vtkPolyData **inputPolys = new vtkPolyData*[nU];
  for (int i=0; i<nU; i++)
  {
    inputPolys[i] = vtkPolyData::New();
    vtkNew(vtkPoints, newPoints);
    inputPolys[i]->SetPoints(newPoints);
    inputPolys[i]->GetPoints()->SetNumberOfPoints(nV);
  }

  vtkNew(vtkDoubleArray, xdata);
  xdata->SetNumberOfTuples(nU);
  vtkNew(vtkDoubleArray, ydata);
  ydata->SetNumberOfTuples(nV);
  vtkSVNURBSUtils::LinSpace(0, 2*M_PI, nU, xdata);
  vtkSVNURBSUtils::LinSpace(0, 2*M_PI, nV, ydata);

  vtkNew(vtkDoubleArray, uders);
  uders->SetNumberOfComponents(3);
  uders->SetNumberOfTuples(nU);
  double uder[3]; uder[0] = 2.0; uder[1] = 0.0; uder[2] = 0.0;
  for (int i=0; i<nU; i++)
  {
    uders->SetTuple(i, uder);
  }

  vtkNew(vtkDoubleArray, vders);
  vders->SetNumberOfComponents(3);
  vders->SetNumberOfTuples(nV);
  double vder[3]; vder[0] = 0.0; vder[1] = 2.0; vder[2] = 0.0;
  for (int i=0; i<nV; i++)
  {
    vders->SetTuple(i, vder);
  }

  for (int i=0; i<nU; i++)
  {
    for (int j=0; j<nV; j++)
    {
      double xval = xdata->GetTuple1(i);
      double yval = ydata->GetTuple1(j);;
      double zval = sin(xval) * sin(yval);

      inputPolys[i]->GetPoints()->SetPoint(j, xval, yval, zval);
    }
  }

  vtkNew(vtkSVLoftNURBSSurface, lofter);
  for (int i=0; i<nU; i++)
  {
    lofter->AddInputData(inputPolys[i]);
  }
  lofter->SetUDegree(3);
  lofter->SetVDegree(3);
  lofter->SetPolyDataUSpacing(0.05);
  lofter->SetPolyDataVSpacing(0.05);
  lofter->SetUKnotSpanType("average");
  lofter->SetUParametricSpanType("chord");
  lofter->SetStartUDerivatives(uders);
  lofter->SetEndUDerivatives(uders);
  lofter->SetVKnotSpanType("derivative");
  lofter->SetVParametricSpanType("centripetal");
  lofter->SetStartVDerivatives(vders);
  lofter->SetEndVDerivatives(vders);
  lofter->Update();

  for (int i=0; i<nU; i++)
  {
    inputPolys[i]->Delete();
  }
  delete [] inputPolys;

  std::string newDirName = getcwd(NULL, 0);
  system(("mkdir -p "+newDirName+"/../Tests").c_str());
  //Write Files
  std::cout<<"Done...Writing Files..."<<endl;
  vtkSVIOUtils::WriteVTPFile(newDirName+"/../Tests/Loft_Surface.vtp",lofter->GetOutput(0),"");
  vtkSVIOUtils::WriteVTSFile(newDirName+"/../Tests/Loft_SurfaceControlPoints.vts",lofter->GetSurface()->GetControlPointGrid(),"");
  std::cout<<"Done"<<endl;

  return EXIT_SUCCESS;
}
