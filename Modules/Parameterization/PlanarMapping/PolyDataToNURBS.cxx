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
#include "vtkPolyDataToNURBSFilter.h"
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
  vtkSmartPointer<vtkPolyDataToNURBSFilter> Converter =
	  vtkSmartPointer<vtkPolyDataToNURBSFilter>::New();

  //Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  ReadInputFile(inputFilename1,pd1);
  ReadInputFile(inputFilename2,pd2);
  ReadInputFile(inputFilename3,pd3);
  ReadInputFile(inputFilename4,pd4);

  //OPERATION
  std::string newDirName = getPath(inputFilename1)+"/"+getRawName(inputFilename1);
  std::string newOutName = getPath(inputFilename1)+"/"+getRawName(inputFilename1)+"/"+getRawName(inputFilename1);
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
  WriteVTPFile(newOutName+".vtp",Converter->GetOutput(), attachName);
  std::string attachName2 = "_Textured";
  WriteVTPFile(newOutName+".vtp",Converter->GetTexturedPd(), attachName2);
  std::string attachName3 = "_Lofted";
  WriteVTPFile(newOutName+".vtp",Converter->GetLoftedPd(), attachName3);
  std::cout<<"Done"<<endl;

  //Exit the program without errors
  return EXIT_SUCCESS;
}
