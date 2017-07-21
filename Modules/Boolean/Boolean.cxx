//
//  Boolean.cxx
//
//
//  Created by Adam Updegrove on 12/17/14.
//
//

/*=========================================================================

 Program:   SimVascular

 =========================================================================*/

#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkSVLoopBooleanPolyDataFilter.h"
#include "vtkSVLoopIntersectionPolyDataFilter.h"

#include <string>
#include <sstream>
#include <iostream>

int main(int argc, char *argv[])
{
  if (argc != 4)
  {
      std::cout << "Input Filenames Required, Need Boolean!" <<endl;
      return EXIT_FAILURE;
  }

  //Create string from input File Name on command line
  std::string inputFilename1 = argv[1];
  std::string inputFilename2 = argv[2];
  int op = (int) (argv[3][0] - '0');

  //Create pointers for reading the STL, creating the full unstructured grid,
  //creating the full poly data, and create the region poly data sets
  vtkNew(vtkPolyData, pd1);
  vtkNew(vtkPolyData, pd2);
  vtkNew(vtkSVLoopBooleanPolyDataFilter, myBoolean);

  //Call Function to Read Fie
  double testTolerance = 1.0e-6; // Doesn't matter anymore!!!
  std::cout<<"Reading Files..."<<endl;
  vtkSVIOUtils::ReadInputFile(inputFilename1,pd1);
  vtkSVIOUtils::ReadInputFile(inputFilename2,pd2);

  //BOOLEAN OPERATION EMBEDDED INTERSECTION
  vtkNew(vtkPolyData, fullpd);
  myBoolean->SetInputData(0,pd1);
  myBoolean->SetInputData(1,pd2);
  myBoolean->SetTolerance(testTolerance);

  if (op == 0)
    myBoolean->SetOperationToUnion();
  else if (op == 1)
    myBoolean->SetOperationToIntersection();
  else if (op == 2)
    myBoolean->SetOperationToDifference();

  fprintf(stdout,"Running Boolean...\n");
  myBoolean->Update();
  fullpd->DeepCopy(myBoolean->GetOutput());
  fprintf(stdout,"Done with Boolean\n");

  fullpd->GetCellData()->RemoveArray("BadTriangle");
  fullpd->GetCellData()->RemoveArray("FreeEdge");

  double dummy[2];
  vtkSVLoopIntersectionPolyDataFilter::CleanAndCheckSurface(fullpd,dummy,testTolerance);
  double fullbadtri[2], fullfreeedge[2];
  fullpd->GetCellData()->GetArray("BadTriangle")->GetRange(fullbadtri);
  fullpd->GetCellData()->GetArray("FreeEdge")->GetRange(fullfreeedge);

  std::cout<<"FULL SURFACE BAD TRI MIN: "<<fullbadtri[0]<<" MAX: "<<fullbadtri[1]<<endl;
  std::cout<<"FULL SURFACE FREE EDGE MIN: "<<fullfreeedge[0]<<" MAX: "<<fullfreeedge[1]<<endl;

  if (op == 0)
    vtkSVIOUtils::WriteVTPFile(inputFilename1,fullpd,"_"+vtkSVIOUtils::GetRawName(inputFilename2)+"_Union");
  else if (op == 1)
    vtkSVIOUtils::WriteVTPFile(inputFilename1,fullpd,"_"+vtkSVIOUtils::GetRawName(inputFilename2)+"_Intersection");
  else if (op == 2)
    vtkSVIOUtils::WriteVTPFile(inputFilename1,fullpd,"_"+vtkSVIOUtils::GetRawName(inputFilename2)+"_Difference");

  vtkSVIOUtils::WriteVTPFile(inputFilename1,myBoolean->GetOutput(1),"_"+vtkSVIOUtils::GetRawName(inputFilename2)+"_IntersectionLines");
  //Exit the program without errors
  return EXIT_SUCCESS;
}




