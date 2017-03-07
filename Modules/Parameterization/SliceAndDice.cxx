//
//  SliceAndDice.cxx
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
#include "vtkDoubleArray.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkIntArray.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkSVIOUtils.h"
#include "vtkSVPolyDataSliceAndDiceFilter.h"

#include <string>
#include <sstream>
#include <iostream>

int main(int argc, char *argv[])
{
  if (argc != 3)
  {
      std::cout << "Need two objects: [PolyData with Ids] [Centerlines]!" <<endl;
      return EXIT_FAILURE;
  }

  //Create string from input File Name on command line
  std::string inputFilename1 = argv[1];
  std::string inputFilename2 = argv[2];

  //creating the full poly data to read in from file and the operation filter
  vtkSmartPointer<vtkPolyData> pd1 = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyData> pd2 = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkSVPolyDataSliceAndDiceFilter> Slicer =
	  vtkSmartPointer<vtkSVPolyDataSliceAndDiceFilter>::New();

  //Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkSVIOUtils::ReadInputFile(inputFilename1,pd1);
  vtkSVIOUtils::ReadInputFile(inputFilename2,pd2);

  //OPERATION
  std::string newDirName = vtkSVIOUtils::GetPath(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1);
  std::string newOutName = vtkSVIOUtils::GetPath(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1);
  system(("mkdir -p "+newDirName).c_str());
  std::cout<<"Performing Operation..."<<endl;
  Slicer->SetInputData(pd1);
  Slicer->SetCenterlines(pd2);
  Slicer->SetSliceLength(1.0);
  Slicer->SetConstructPolycube(1);
  Slicer->SetBoundaryPointsArrayName("BoundaryPoints");
  Slicer->SetGroupIdsArrayName("GroupIds");
  Slicer->SetSegmentIdsArrayName("SegmentIds");
  Slicer->SetSliceIdsArrayName("SliceIds");
  Slicer->SetSphereRadiusArrayName("MaximumInscribedSphereRadius");
  Slicer->SetInternalIdsArrayName("TmpInternalIds");
  Slicer->SetDijkstraArrayName("DijkstraDistance");
  Slicer->Update();

  vtkSmartPointer<vtkPolyData> output =
    vtkSmartPointer<vtkPolyData>::New();
  output = Slicer->GetOutput();

  //Write Files
  std::cout<<"Done...Writing Files..."<<endl;
  std::string attachName = "_Segmented";
  vtkSVIOUtils::WriteVTPFile(newOutName+".vtp",output,attachName);
  std::cout<<"Done"<<endl;

  //Exit the program without errors
  return EXIT_SUCCESS;
}
