//
//  SphericalConformalChecker.cxx
//
//
//  Created by Adam Updegrove on 10/4/14.
//
//

/*=========================================================================

 Program:   SimVascular

 =========================================================================*/

#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkDataArray.h"
#include "vtkDataWriter.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkIntArray.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkSVCheckRotation.h"
#include "vtkSTLReader.h"
#include "vtkTriangle.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLPolyDataReader.h"

#include <string>
#include <sstream>
#include <iostream>

int main(int argc, char *argv[])
{
  if (argc < 3 || argc > 4)
  {
      std::cout << "Incorrect Usage! Should be:" <<endl;
      std::cout << "./SphericalConformalChecker [filename 1] [filename 2] [optional original]" <<endl;
      return EXIT_FAILURE;
  }

  //Create string from input File Name on command line
  std::string inputFilename1 = argv[1];
  std::string inputFilename2 = argv[2];
  std::string inputFilename3;
  if (argc == 4)
  {
    inputFilename3 = argv[3];
  }

  //creating the full poly data to read in from file and the operation filter
  vtkNew(vtkPolyData, pd1);
  vtkNew(vtkPolyData, pd2);
  vtkNew(vtkPolyData, pd3);
  vtkNew(vtkSVCheckRotation, Checker);

  //Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkSVIOUtils::ReadInputFile(inputFilename1,pd1);
  vtkSVIOUtils::ReadInputFile(inputFilename2,pd2);
  if (!inputFilename3.empty())
  {
    vtkSVIOUtils::ReadInputFile(inputFilename3,pd3);
  }

  std::string newDirName = vtkSVIOUtils::GetPath(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1);
  std::string newOutName = vtkSVIOUtils::GetPath(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1);
  system(("mkdir -p "+newDirName).c_str());
  //OPERATION
  std::cout<<"Performing Operation..."<<endl;
  Checker->SetInputData(0, pd1);
  Checker->SetInputData(1, pd2);
  if (!inputFilename3.empty())
  {
    Checker->SetOriginalPd(pd3);
  }
  Checker->SetVerbose(1);
  Checker->Update();

  //Write Files
  std::cout<<"Done"<<endl;
  vtkSVIOUtils::WriteVTPFile(newOutName+".vtp",Checker->GetOutput(0),"_CheckRotate");

  //Exit the program without errors
  return EXIT_SUCCESS;
}
