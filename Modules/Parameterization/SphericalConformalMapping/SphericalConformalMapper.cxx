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
#include "vtkSphericalConformalMapper.h"
#include "vtkSTLReader.h"
#include "vtkTriangle.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLPolyDataReader.h"

#include <string>
#include <sstream>
#include <iostream>

//Function to turn an integer into a string
std::string intToString(int i)
{
  std::stringstream out;
  out << i;
  return out.str();
}

//Function to get the directory from the input File Name
//For example, /User/Adam.stl returns /User
std::string getPath(std::string fullName)
{
  std::string pathName;
  unsigned split = fullName.find_last_of("/\\");
  pathName = fullName.substr(0,split);
  return pathName;
}

//Function to get the raw file name from the input File name
//For example, Adam.stl returns Adam
std::string getRawName(std::string fullName)
{
  std::string rawName;
  unsigned split = fullName.find_last_of("/\\");
  rawName = fullName.substr(split+1);
  rawName.erase(rawName.find_last_of("."),std::string::npos);
  return rawName;
}

std::string getExt(std::string fullName)
{
  std::string extName;
  unsigned split = fullName.find_last_of(".");
  extName = fullName.substr(split+1);
  return extName;
}

//Function to read in the STL file, extract the boundaries and pass the input
//Poly Data information
void ReadSTLFile(std::string inputFilename, vtkPolyData *polydata)
{
  //Create an STL reader for reading the file
  vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
  reader->SetFileName(inputFilename.c_str());
  reader->Update();

  //Save the output information from the boundary filter to a Poly Data
  //structure
  polydata->DeepCopy(reader->GetOutput());
  polydata->BuildLinks();
}

void ReadVTPFile(std::string inputFilename, vtkPolyData *polydata)
{
  //Create an STL reader for reading the file
  vtkSmartPointer<vtkXMLPolyDataReader> reader =
	  vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(inputFilename.c_str());
  reader->Update();

  //Save the output information from the boundary filter to a Poly Data
  //structure
  polydata->DeepCopy(reader->GetOutput());
  polydata->BuildLinks();
}

void ReadInputFile(std::string inputFilename, vtkPolyData *polydata)
{
  std::string ext = getExt(inputFilename);
  std:cout<<"Extension... "<<ext<<endl;
  if(!strncmp(ext.c_str(),"stl",3))
  {
    ReadSTLFile(inputFilename, polydata);
  }
  else if(!strncmp(ext.c_str(),"vtp",3))
  {
    ReadVTPFile(inputFilename, polydata);
  }
  else
  {
    std::cout<<"Unrecognized file extension, stl and vtp accepted"<<endl;
  }
}

//Function to write the polydata to a vtp
void WriteVTPFile(std::string inputFilename,vtkPolyData *writePolyData,std::string attachName)
{
  std::string rawName, pathName, outputFilename;

  vtkSmartPointer<vtkXMLPolyDataWriter> writer  = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

  pathName = getPath(inputFilename);
  rawName = getRawName(inputFilename);

  outputFilename = pathName+"/"+rawName+attachName+".vtp";

  writer->SetFileName(outputFilename.c_str());
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(writePolyData);
#else
  writer->SetInputData(writePolyData);
#endif
  //writer->SetDataModeToAscii();

  writer->Write();
}

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
  vtkSmartPointer<vtkPolyData> pd1 = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkSphericalConformalMapper> Mapper =
	  vtkSmartPointer<vtkSphericalConformalMapper>::New();

  //Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  ReadInputFile(inputFilename1,pd1);

  //SUPPLY BOUNDARY START POINTS!!!
  vtkSmartPointer<vtkIntArray> firstLoopPts =
    vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> firstLoopHelper =
    vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> secondLoopPts =
    vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> secondLoopHelper =
    vtkSmartPointer<vtkIntArray>::New();
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

  std::string newDirName = getPath(inputFilename1)+"/"+getRawName(inputFilename1);
  std::string newOutName = getPath(inputFilename1)+"/"+getRawName(inputFilename1)+"/"+getRawName(inputFilename1);
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
  WriteVTPFile(newOutName+".vtp",Mapper->GetOutput(0),"_Mapped_"+intToString(conj_method));
  std::cout<<"Done"<<endl;

  //Exit the program without errors
  return EXIT_SUCCESS;
}
