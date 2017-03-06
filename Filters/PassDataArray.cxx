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
#include "vtkSVPassDataArray.h"
#include "vtkSTLReader.h"
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
  if (argc != 4)
  {
      std::cout << "Need two surfaces and array name: ./PassDataArray [Source Surface] [Target Surface] [Array Name]!" <<endl;
      return EXIT_FAILURE;
  }

  //Create string from input File Name on command line
  std::string inputFilename1 = argv[1];
  std::string inputFilename2 = argv[2];
  std::string arrayName      = argv[3];

  //creating the full poly data to read in from file and the operation filter
  vtkSmartPointer<vtkPolyData> pd1 = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyData> pd2 = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkSVPassDataArray> Passer =
	  vtkSmartPointer<vtkSVPassDataArray>::New();

  //Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  ReadInputFile(inputFilename1,pd1);
  ReadInputFile(inputFilename2,pd2);

  //OPERATION
  std::string newDirName = getPath(inputFilename2)+"/"+getRawName(inputFilename2);
  std::string newOutName = getPath(inputFilename2)+"/"+getRawName(inputFilename2)+"/"+getRawName(inputFilename2);
  system(("mkdir -p "+newDirName).c_str());
  std::cout<<"Performing Operation..."<<endl;
  Passer->SetInputData(0, pd1);
  Passer->SetInputData(1, pd2);
  Passer->SetPassArrayName(arrayName.c_str());
  Passer->SetPassDataIsCellData(1);
  Passer->SetPassDataToCellData(1);
  Passer->Update();

  //Write Files
  std::cout<<"Done...Writing Files..."<<endl;
  std::string attachName = "_With_"+arrayName;
  WriteVTPFile(newOutName+".vtp",Passer->GetOutput(0),attachName);
  std::cout<<"Done"<<endl;

  //Exit the program without errors
  return EXIT_SUCCESS;
}
