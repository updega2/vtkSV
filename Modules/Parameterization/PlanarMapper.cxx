//
//  SphericalConformalMapper.cxx
//
//
//  Created by Adam Updegrove on 6/4/16.
//
//

/*=========================================================================

 Program:   SimVascular

 =========================================================================*/

#include "vtkSVSquareBoundaryMapper.h"
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
#include "vtkSVIOUtils.h"
#include "vtkSVPlanarMapper.h"

#include <string>
#include <sstream>
#include <iostream>

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
      std::cout << "Incorrect Usage! Should be:" <<endl;
      std::cout << "./PlanarMapper [filename]" <<endl;
      return EXIT_FAILURE;
  }

  //Create string from input File Name on command line
  std::string inputFilename1 = argv[1];

  //creating the full poly data to read in from file and the operation filter
  vtkSmartPointer<vtkPolyData> pd1 = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkSVPlanarMapper> Mapper =
	  vtkSmartPointer<vtkSVPlanarMapper>::New();

  //Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkSVIOUtils::ReadInputFile(inputFilename1,pd1);

  vtkSmartPointer<vtkIntArray> boundaryCorners =
    vtkSmartPointer<vtkIntArray>::New();
  boundaryCorners->SetNumberOfComponents(1);
  boundaryCorners->SetNumberOfTuples(4);
  // 0103_0001
  boundaryCorners->SetValue(0,4081);
  boundaryCorners->SetValue(1,370);
  boundaryCorners->SetValue(2,10);
  boundaryCorners->SetValue(3,4016);
  // 0110_0001
  //boundaryCorners->SetValue(0,5642);
  //boundaryCorners->SetValue(1,3624);
  //boundaryCorners->SetValue(2,4742);
  //boundaryCorners->SetValue(3,5610);
  // HalfSphere
  //boundaryCorners->SetValue(0,15);
  //boundaryCorners->SetValue(1,24);
  //boundaryCorners->SetValue(2,29);
  //boundaryCorners->SetValue(3,9);
  // HalfSphere 2
  //boundaryCorners->SetValue(0,32);
  //boundaryCorners->SetValue(1,104);
  //boundaryCorners->SetValue(2,67);
  //boundaryCorners->SetValue(3,72);
  // Aorta
  //boundaryCorners->SetValue(0,22320);
  //boundaryCorners->SetValue(1,22691);
  //boundaryCorners->SetValue(2,22299);
  //boundaryCorners->SetValue(3,23);
  // IliacBranchSegme
  //boundaryCorners->SetValue(0,185);
  //boundaryCorners->SetValue(1,220);
  //boundaryCorners->SetValue(2,213);
  //boundaryCorners->SetValue(3,164);
  vtkSmartPointer<vtkSVSquareBoundaryMapper> boundaryMapper =
    vtkSmartPointer<vtkSVSquareBoundaryMapper>::New();
  boundaryMapper->SetBoundaryIds(boundaryCorners);

  std::string newDirName = vtkSVIOUtils::GetPath(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1);
  std::string newOutName = vtkSVIOUtils::GetPath(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1);
  system(("mkdir -p "+newDirName).c_str());
  //OPERATION
  std::cout<<"Performing Operation..."<<endl;
  Mapper->SetInputData(pd1);
  Mapper->SetBoundaryMapper(boundaryMapper);
  Mapper->Update();

  //Write Files
  std::cout<<"Done...Writing Files..."<<endl;
  vtkSVIOUtils::WriteVTPFile(newOutName+".vtp",Mapper->GetOutput(0),"_Mapped");
  std::cout<<"Done"<<endl;

  //Exit the program without errors
  return EXIT_SUCCESS;
}
