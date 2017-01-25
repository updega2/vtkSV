//
//  Intersect.cxx
//
//
//  Created by Adam Updegrove on 10/4/14.
//
//

/*=========================================================================

 Program:   SimVascular

 =========================================================================*/

#include "vtkSTLReader.h"
#include "vtkOBBTree.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkDistancePolyDataFilter.h"
#include "vtkIntersectionPolyDataFilter2.h"
#include "vtkIntersectionPolyDataFilter.h"
#include "vtkPolyData.h"
#include "vtkDataArray.h"
#include "vtkIntArray.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkCellArray.h"
#include "vtkIdList.h"
#include "vtkDataWriter.h"
#include "vtkSmartPointer.h"
#include "vtkInformation.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkFillHolesFilter.h"
#include "vtkPolyDataNormals.h"
#include "vtkCleanPolyData.h"

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
      std::cout << "Input Filenames Required!" <<endl;
      return EXIT_FAILURE;
  }

  //Create string from input File Name on command line
  std::string inputFilename1 = argv[1];
  std::string inputFilename2 = argv[2];

  //Create pointers for reading the STL, creating the full unstructured grid,
  //creating the full poly data, and create the region poly data sets
  vtkSmartPointer<vtkPolyData> pd1 = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyData> pd2 = vtkSmartPointer<vtkPolyData>::New();
#ifdef USE_MINE
  vtkSmartPointer<vtkIntersectionPolyDataFilter2> PolyDataIntersection =
	  vtkSmartPointer<vtkIntersectionPolyDataFilter2>::New();
#else
  vtkSmartPointer<vtkIntersectionPolyDataFilter> PolyDataIntersection =
	  vtkSmartPointer<vtkIntersectionPolyDataFilter>::New();
#endif

  //Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  ReadSTLFile(inputFilename1,pd1);
  ReadSTLFile(inputFilename2,pd2);

  //INTERSECTION OPERATION
  std::cout<<"Performing Intersection..."<<endl;
  PolyDataIntersection->AddInputData(0,pd1);
  PolyDataIntersection->AddInputData(1,pd2);
  PolyDataIntersection->SplitFirstOutputOn();
  PolyDataIntersection->SplitSecondOutputOn();
  PolyDataIntersection->Update();

  std::cout<<"Done...Writing Files..."<<endl;
  WriteVTPFile(inputFilename1,PolyDataIntersection->GetOutput(0),"_IntersectionLines");
  WriteVTPFile(inputFilename1,PolyDataIntersection->GetOutput(1),"_IntersectionObject1");
  WriteVTPFile(inputFilename2,PolyDataIntersection->GetOutput(2),"_IntersectionObject2");
  std::cout<<"Done"<<endl;

  //Exit the program without errors
  return EXIT_SUCCESS;
}




