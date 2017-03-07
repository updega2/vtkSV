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
#include "vtkSVHausdorffDistance.h"
#include "vtkSVIOUtils.h"

#include <string>
#include <sstream>
#include <iostream>

int main(int argc, char *argv[])
{
  if (argc != 3)
  {
      std::cout << "Need two surfaces: ./HausdorffDistance [Source Surface] [Target Surface]!" <<endl;
      return EXIT_FAILURE;
  }

  //Create string from input File Name on command line
  std::string inputFilename1 = argv[1];
  std::string inputFilename2 = argv[2];

  //creating the full poly data to read in from file and the operation filter
  vtkSmartPointer<vtkPolyData> pd1 = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyData> pd2 = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkSVHausdorffDistance> Distancer =
	  vtkSmartPointer<vtkSVHausdorffDistance>::New();

  //Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkSVIOUtils::ReadInputFile(inputFilename1,pd1);
  vtkSVIOUtils::ReadInputFile(inputFilename2,pd2);

  //OPERATION
  std::string newDirName = vtkSVIOUtils::GetPath(inputFilename2)+"/"+vtkSVIOUtils::GetRawName(inputFilename2);
  std::string newOutName = vtkSVIOUtils::GetPath(inputFilename2)+"/"+vtkSVIOUtils::GetRawName(inputFilename2)+"/"+vtkSVIOUtils::GetRawName(inputFilename2);
  system(("mkdir -p "+newDirName).c_str());
  std::cout<<"Performing Operation..."<<endl;
  Distancer->SetInputData(0, pd1);
  Distancer->SetInputData(1, pd2);
  Distancer->SetDistanceArrayName("Distance");
  Distancer->Update();

  std::cout<<"Hausdorff Distance: "<<Distancer->GetHausdorffDistance()<<endl;
  std::cout<<"Average Distance:   "<<Distancer->GetAverageDistance()<<endl;
  //Write Files
  std::cout<<"Done...Writing Files..."<<endl;
  std::string attachName = "_Distanced";
  vtkSVIOUtils::WriteVTPFile(newOutName+".vtp",Distancer->GetOutput(0),attachName);
  std::cout<<"Done"<<endl;

  //Exit the program without errors
  return EXIT_SUCCESS;
}
