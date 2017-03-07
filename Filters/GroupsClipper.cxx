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
#include "vtkInformation.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkSVIOUtils.h"
#include "vtkSVPolyDataCenterlineGroupsClipper.h"

int main(int argc, char *argv[])
{
  if (argc != 3)
  {
      std::cout << "Need two surfaces: ./HausdorffDistance [Surface] [Centerlines]!" <<endl;
      return EXIT_FAILURE;
  }

  //Create string from input File Name on command line
  std::string inputFilename1 = argv[1];
  std::string inputFilename2 = argv[2];

  //creating the full poly data to read in from file and the operation filter
  vtkSmartPointer<vtkPolyData> pd1 = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyData> pd2 = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkSVPolyDataCenterlineGroupsClipper> Grouper =
	  vtkSmartPointer<vtkSVPolyDataCenterlineGroupsClipper>::New();

  //Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkSVIOUtils::ReadInputFile(inputFilename1,pd1);
  vtkSVIOUtils::ReadInputFile(inputFilename2,pd2);

  //OPERATION
  std::string newDirName = vtkSVIOUtils::GetPath(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1);
  std::string newOutName = vtkSVIOUtils::GetPath(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1);
  system(("mkdir -p "+newDirName).c_str());
  std::cout<<"Performing Operation..."<<endl;
  Grouper->SetInputData(pd1);
  Grouper->SetCenterlines(pd2);
  Grouper->SetCenterlineGroupIdsArrayName("GroupIds");
  Grouper->SetGroupIdsArrayName("GroupIds");
  Grouper->SetCenterlineRadiusArrayName("MaximumInscribedSphereRadius");
  Grouper->SetBlankingArrayName("Blanking");
  Grouper->SetCutoffRadiusFactor(1.0e16);
  Grouper->SetClipValue(0.0);
  Grouper->ClipAllCenterlineGroupIdsOn();
  Grouper->SetUseRadiusInformation(1);
  Grouper->Update();

  //Write Files
  std::cout<<"Done...Writing Files..."<<endl;
  std::string attachName = "_Grouped";
  vtkSVIOUtils::WriteVTPFile(newOutName+".vtp",Grouper->GetOutput(0),attachName);
  std::cout<<"Done"<<endl;

  //Exit the program without errors
  return EXIT_SUCCESS;
}
