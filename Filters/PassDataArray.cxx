//
//  PassDataArray.cxx
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
#include "vtkSVPassDataArray.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"

#include <string>
#include <sstream>
#include <iostream>

int main(int argc, char *argv[])
{
  if (argc != 4)
  {
      std::cout << "Need two surfaces and array name: ./PassDataArray [Source Surface] [Target Surface] [Array Name]!" <<endl;
      return EXIT_FAILURE;
  }

  //Create string from input File Name on command line
  std::string inputFilename1 = argv[1];
  std::string inputFilename2 = argv[2];
  std::string arrayName      = argv[3];

  //creating the full poly data to read in from file and the operation filter
  vtkNew(vtkPolyData, pd1);
  vtkNew(vtkPolyData, pd2);
  vtkNew(vtkSVPassDataArray, Passer);

  //Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkSVIOUtils::ReadInputFile(inputFilename1,pd1);
  vtkSVIOUtils::ReadInputFile(inputFilename2,pd2);

  //OPERATION
  std::string newDirName = vtkSVIOUtils::GetPath(inputFilename2)+"/"+vtkSVIOUtils::GetRawName(inputFilename2);
  std::string newOutName = vtkSVIOUtils::GetPath(inputFilename2)+"/"+vtkSVIOUtils::GetRawName(inputFilename2)+"/"+vtkSVIOUtils::GetRawName(inputFilename2);
  system(("mkdir -p "+newDirName).c_str());
  std::cout<<"Performing Operation..."<<endl;
  Passer->SetInputData(0, pd1);
  Passer->SetInputData(1, pd2);
  Passer->SetPassArrayName(arrayName.c_str());
  Passer->SetPassDataIsCellData(1);
  Passer->SetPassDataToCellData(1);
  Passer->Update();

  //Write Files
  std::cout<<"Done...Writing Files..."<<endl;
  std::string attachName = "_With_"+arrayName;
  vtkSVIOUtils::WriteVTPFile(newOutName+".vtp",Passer->GetOutput(0),attachName);
  std::cout<<"Done"<<endl;

  //Exit the program without errors
  return EXIT_SUCCESS;
}
