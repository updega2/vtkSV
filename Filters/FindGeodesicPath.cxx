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
#include "vtkSVFindGeodesicPath.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkTriangle.h"

#include <string>
#include <sstream>
#include <iostream>

int main(int argc, char *argv[])
{
  if (argc < 3 || argc > 4)
  {
      std::cout << "Incorrect Usage! Should be:" <<endl;
      std::cout << "./FindClosestGeodesicPoint [filename] [start point id] [optional: end point id]" <<endl;
      return EXIT_FAILURE;
  }

  //Create string from input File Name on command line
  std::string inputFilename1 = argv[1];
  int startPointId = atoi(argv[2]);
  int endPointId = -1;
  if (argc == 4)
  {
    endPointId = atoi(argv[3]);
  }

  //creating the full poly data to read in from file and the operation filter
  vtkNew(vtkPolyData, pd1);
  vtkNew(vtkSVFindGeodesicPath, Finder);

  //Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkSVIOUtils::ReadInputFile(inputFilename1,pd1);

  std::string newDirName = vtkSVIOUtils::GetPath(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1);
  std::string newOutName = vtkSVIOUtils::GetPath(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1);
  system(("mkdir -p "+newDirName).c_str());
  //OPERATION
  std::cout<<"Performing Operation..."<<endl;
  Finder->SetInputData(0, pd1);
  Finder->SetStartPtId(startPointId);
  Finder->SetEndPtId(endPointId);
  Finder->SetAddPathBooleanArray(1);
  Finder->SetDijkstraArrayName("DijkstraDistance");
  Finder->SetInternalIdsArrayName("InternalIds");
  Finder->SetPathBooleanArrayName("IsPath");
  //Finder->SetRepelCloseBoundaryPoints(1);
  Finder->SetVerbose(1);
  Finder->Update();

  //Write Files
  std::cout<<"Done"<<endl;
  vtkSVIOUtils::WriteVTPFile(newOutName+".vtp", Finder->GetOutput(0),"_Geodescized");

  //Exit the program without errors
  return EXIT_SUCCESS;
}
