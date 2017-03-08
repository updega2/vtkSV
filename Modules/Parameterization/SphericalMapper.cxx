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
#include "vtkSVSphericalMapper.h"
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
  if (argc != 3)
  {
      std::cout << "Incorrect Usage! Should be:" <<endl;
      std::cout << "./SphericalConformalMapper [filename] [cg_update_method]" <<endl;
      std::cout << "...[cg_update_method] -> 0 = CG_NONE (steepest descent)" <<endl;
      std::cout << "...[cg_update_method] -> 1 = CG_FLETCHER_REEVES" <<endl;
      std::cout << "...[cg_update_method] -> 2 = CG_POLAK_RIBIER" <<endl;
      std::cout << "...[cg_update_method] -> 3 = CG_HESTENESS_STIEFEL" <<endl;
      std::cout << "...[cg_update_method] -> 4 = CG_DAI_YUAN" <<endl;
      return EXIT_FAILURE;
  }

  //Create string from input File Name on command line
  std::string inputFilename1 = argv[1];
  int conj_method = (int) (argv[2][0] - '0');

  //creating the full poly data to read in from file and the operation filter
  vtkNew(vtkPolyData, pd1);
  vtkNew(vtkSVSphericalMapper, Mapper);

  //Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkSVIOUtils::ReadInputFile(inputFilename1,pd1);

  //SUPPLY BOUNDARY START POINTS!!!
  vtkNew(vtkIntArray, firstLoopPts);
  vtkNew(vtkIntArray, firstLoopHelper);
  vtkNew(vtkIntArray, secondLoopPts);
  vtkNew(vtkIntArray, secondLoopHelper);
  int starts[2]; starts[0] = 687; starts[1] = -1;
  firstLoopPts->InsertNextValue(242);  firstLoopHelper->InsertNextValue(1);
  firstLoopPts->InsertNextValue(634);  firstLoopHelper->InsertNextValue(0);
  firstLoopPts->InsertNextValue(676);  firstLoopHelper->InsertNextValue(4);
  firstLoopPts->InsertNextValue(684);  firstLoopHelper->InsertNextValue(2);
  firstLoopPts->InsertNextValue(113);  firstLoopHelper->InsertNextValue(1);
  firstLoopPts->InsertNextValue(34);  firstLoopHelper->InsertNextValue(3);
  firstLoopPts->InsertNextValue(608);  firstLoopHelper->InsertNextValue(4);
  firstLoopPts->InsertNextValue(687);  firstLoopHelper->InsertNextValue(5);
  //secondLoopPts->InsertNextValue(420); secondLoopHelper->InsertNextValue(0);
  //secondLoopPts->InsertNextValue(362);  secondLoopHelper->InsertNextValue(1);
  //secondLoopPts->InsertNextValue(326);  secondLoopHelper->InsertNextValue(3);
  //secondLoopPts->InsertNextValue(310); secondLoopHelper->InsertNextValue(4);
  int cubeStart[2]; cubeStart[0] = 0; cubeStart[1] = 4;

  std::string newDirName = vtkSVIOUtils::GetPath(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1);
  std::string newOutName = vtkSVIOUtils::GetPath(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1);
  system(("mkdir -p "+newDirName).c_str());
  //OPERATION
  std::cout<<"Performing Operation..."<<endl;
  Mapper->SetInputData(pd1);
  Mapper->SetVerbose(1);
  Mapper->SetInitialTimeStep(0.1);
  //Mapper->SetTutteEnergyCriterion(0.0001);
  //Mapper->SetHarmonicEnergyCriterion(0.00000001);
  Mapper->SetTutteEnergyCriterion(1.0e-6);
  Mapper->SetHarmonicEnergyCriterion(1.0e-7);
  Mapper->SetMaxNumIterations(1e4);
  Mapper->SetIterOutputFilename(newOutName.c_str());
  Mapper->SetNumSaveIterations(1000);
  Mapper->SetCGUpdateMethod(conj_method);
  Mapper->SetBoundaryStart(starts);
  Mapper->SetObjectXAxis(1.0, 0.0, 0.0);
  Mapper->SetObjectZAxis(0.0, 0.0, 1.0);
  Mapper->SetFirstLoopPts(firstLoopPts);
  //Mapper->SetSecondLoopPts(secondLoopPts);
  Mapper->SetFirstLoopHelper(firstLoopHelper);
  //Mapper->SetSecondLoopHelper(secondLoopHelper);
  Mapper->SetCubeStart(cubeStart);
  Mapper->Update();

  //Write Files
  std::cout<<"Done...Writing Files..."<<endl;
  vtkSVIOUtils::WriteVTPFile(newOutName+".vtp",Mapper->GetOutput(0),"_Mapped_"+vtkSVIOUtils::IntToString(conj_method));
  std::cout<<"Done"<<endl;

  //Exit the program without errors
  return EXIT_SUCCESS;
}
