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
#include "vtkSVPlanarMapper.h"
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
  ReadInputFile(inputFilename1,pd1);

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

  std::string newDirName = getPath(inputFilename1)+"/"+getRawName(inputFilename1);
  std::string newOutName = getPath(inputFilename1)+"/"+getRawName(inputFilename1)+"/"+getRawName(inputFilename1);
  system(("mkdir -p "+newDirName).c_str());
  //OPERATION
  std::cout<<"Performing Operation..."<<endl;
  Mapper->SetInputData(pd1);
  Mapper->SetBoundaryMapper(boundaryMapper);
  Mapper->Update();

  //Write Files
  std::cout<<"Done...Writing Files..."<<endl;
  WriteVTPFile(newOutName+".vtp",Mapper->GetOutput(0),"_Mapped");
  std::cout<<"Done"<<endl;

  //Exit the program without errors
  return EXIT_SUCCESS;
}
