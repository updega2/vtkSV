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

#include "vtkSVSuperSquareBoundaryMapper.h"
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
#include "vtkSVPlanarMapper.h"
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
      std::cout << "./PlanarMapper [filename]" <<endl;
      return EXIT_FAILURE;
  }

  //Create string from input File Name on command line
  std::string inputFilename1 = argv[1];

  //creating the full poly data to read in from file and the operation filter
  vtkNew(vtkPolyData, pd1);
  vtkNew(vtkSVPlanarMapper, Mapper);

  //Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkSVIOUtils::ReadInputFile(inputFilename1,pd1);

  vtkNew(vtkIntArray, boundaryCorners);
  boundaryCorners->SetNumberOfComponents(1);
  boundaryCorners->SetNumberOfTuples(10);
  int boundaryDivisions[4];
  // 0103_0001
  boundaryCorners->SetValue(0,4081);  //Boundary 1
  boundaryCorners->SetValue(1,3127);
  boundaryCorners->SetValue(2,225);
  boundaryCorners->SetValue(3,292);
  boundaryCorners->SetValue(4,370); //Boundary 2
  boundaryCorners->SetValue(5,10);  //Boundary 3
  boundaryCorners->SetValue(6,131);
  boundaryCorners->SetValue(7,145);
  boundaryCorners->SetValue(8,80);
  boundaryCorners->SetValue(9,4016); //Boundary 4
  boundaryDivisions[0] = 3; boundaryDivisions[1] = 0; boundaryDivisions[2] = 3; boundaryDivisions[3] = 0;

  vtkNew(vtkSVSuperSquareBoundaryMapper, boundaryMapper);
  boundaryMapper->SetBoundaryIds(boundaryCorners);
  boundaryMapper->SetSuperBoundaryDivisions(boundaryDivisions);

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
  vtkSVIOUtils::WriteVTPFile(newOutName+".vtp",Mapper->GetOutput(0),"_SuperSquareBoundaryMapped");
  std::cout<<"Done"<<endl;

  //Exit the program without errors
  return EXIT_SUCCESS;
}
