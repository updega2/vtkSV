//
//  PolyDataToNURBS.cxx
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
#include "vtkSVPolyDataToNURBSFilter.h"

#include <string>
#include <sstream>
#include <iostream>

int main(int argc, char *argv[])
{
  if (argc != 5)
  {
      std::cout << "Need four objects: [PolyData with Ids] [Centerlines] [S2 Matching Parameterization] [S2 Matching Open Parameterization]!" <<endl;
      return EXIT_FAILURE;
  }

  //Create string from input File Name on command line
  std::string inputFilename1 = argv[1];
  std::string inputFilename2 = argv[2];
  std::string inputFilename3 = argv[3];
  std::string inputFilename4 = argv[4];

  //creating the full poly data to read in from file and the operation filter
  vtkSmartPointer<vtkPolyData> pd1 = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyData> pd2 = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyData> pd3 = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyData> pd4 = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkSVPolyDataToNURBSFilter> Converter =
	  vtkSmartPointer<vtkSVPolyDataToNURBSFilter>::New();

  //Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkSVIOUtils::ReadInputFile(inputFilename1,pd1);
  vtkSVIOUtils::ReadInputFile(inputFilename2,pd2);
  vtkSVIOUtils::ReadInputFile(inputFilename3,pd3);
  vtkSVIOUtils::ReadInputFile(inputFilename4,pd4);

  //OPERATION
  std::string newDirName = vtkSVIOUtils::GetPath(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1);
  std::string newOutName = vtkSVIOUtils::GetPath(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1);
  system(("mkdir -p "+newDirName).c_str());
  std::cout<<"Performing Operation..."<<endl;
  Converter->SetInputData(pd1);
  Converter->SetCenterlines(pd2);
  Converter->SetCubeS2Pd(pd3);
  Converter->SetOpenCubeS2Pd(pd4);
  Converter->SetAddTextureCoordinates(1);
  Converter->SetBoundaryPointsArrayName("BoundaryPoints");
  Converter->SetGroupIdsArrayName("GroupIds");
  Converter->SetSegmentIdsArrayName("SegmentIds");
  Converter->SetSliceIdsArrayName("SliceIds");
  Converter->SetSphereRadiusArrayName("MaximumInscribedSphereRadius");
  Converter->SetInternalIdsArrayName("TmpInternalIds");
  Converter->SetDijkstraArrayName("DijkstraDistance");
  Converter->SetBooleanPathArrayName("IsPath");
  Converter->Update();

  //Write Files
  std::cout<<"Done...Writing Files..."<<endl;
  std::string attachName = "_Converted";
  vtkSVIOUtils::WriteVTPFile(newOutName+".vtp",Converter->GetOutput(), attachName);
  std::string attachName2 = "_Textured";
  vtkSVIOUtils::WriteVTPFile(newOutName+".vtp",Converter->GetTexturedPd(), attachName2);
  std::string attachName3 = "_Lofted";
  vtkSVIOUtils::WriteVTPFile(newOutName+".vtp",Converter->GetLoftedPd(), attachName3);
  std::cout<<"Done"<<endl;

  //Exit the program without errors
  return EXIT_SUCCESS;
}
