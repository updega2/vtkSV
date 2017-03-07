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
#include "vtkSVPullApartPolyData.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"

#include <string>
#include <sstream>
#include <iostream>

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
      std::cout << "Incorrect Usage! Should be:" <<endl;
      std::cout << "./FindClosestGeodesicPoint [filename]" <<endl;
      return EXIT_FAILURE;
  }

  //Create string from input File Name on command line
  std::string inputFilename1 = argv[1];

  //creating the full poly data to read in from file and the operation filter
  vtkNew(vtkPolyData, pd1);;
  vtkNew(vtkSVPullApartPolyData, Ripper);

  //Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkSVIOUtils::ReadInputFile(inputFilename1,pd1);

  std::string newDirName = vtkSVIOUtils::GetPath(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1);
  std::string newOutName = vtkSVIOUtils::GetPath(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1);
  system(("mkdir -p "+newDirName).c_str());
  //OPERATION
  std::cout<<"Performing Operation..."<<endl;
  Ripper->SetInputData(0, pd1);
  Ripper->SetCutPointsArrayName("IsPath");
  Ripper->Update();

  //Write Files
  std::cout<<"Done"<<endl;
  vtkSVIOUtils::WriteVTPFile(newOutName+".vtp", Ripper->GetOutput(0),"_Ripped");

  //Exit the program without errors
  return EXIT_SUCCESS;
}
