//
//  Boolean.cxx
//
//
//  Created by Adam Updegrove on 12/17/14.
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
#include "vtkSVLoopIntersectionPolyDataFilter.h"
#include "vtkSVLoopBooleanPolyDataFilter.h"
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
      std::cout << "Input Filenames Required, Need Boolean!" <<endl;
      return EXIT_FAILURE;
  }

  //Create string from input File Name on command line
  std::string inputFilename1 = argv[1];
  std::string inputFilename2 = argv[2];
  int op = (int) (argv[3][0] - '0');

  //Create pointers for reading the STL, creating the full unstructured grid,
  //creating the full poly data, and create the region poly data sets
  vtkSmartPointer<vtkPolyData> pd1 = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyData> pd2 = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkSVLoopBooleanPolyDataFilter> myBoolean =
	  vtkSmartPointer<vtkSVLoopBooleanPolyDataFilter>::New();

  //Call Function to Read File
  double testTolerance = 1e-6;
  std::cout<<"Reading Files..."<<endl;
  ReadInputFile(inputFilename1,pd1);
  ReadInputFile(inputFilename2,pd2);

  //BOOLEAN OPERATION EMBEDDED INTERSECTION
  vtkSmartPointer<vtkPolyData> fullpd =
	  vtkSmartPointer<vtkPolyData>::New();
  myBoolean->SetInputData(0,pd1);
  myBoolean->SetInputData(1,pd2);
  myBoolean->SetTolerance(testTolerance);

  if (op == 0)
    myBoolean->SetOperationToUnion();
  else if (op == 1)
    myBoolean->SetOperationToIntersection();
  else if (op == 2)
    myBoolean->SetOperationToDifference();

  fprintf(stdout,"Running Boolean...\n");
  myBoolean->Update();
  fullpd->DeepCopy(myBoolean->GetOutput());
  fprintf(stdout,"Done with Boolean\n");

  fullpd->GetCellData()->RemoveArray("BadTriangle");
  fullpd->GetCellData()->RemoveArray("FreeEdge");

  double dummy[2];
  vtkSVLoopIntersectionPolyDataFilter::CleanAndCheckSurface(fullpd,dummy,testTolerance);
  double fullbadtri[2], fullfreeedge[2];
  fullpd->GetCellData()->GetArray("BadTriangle")->GetRange(fullbadtri);
  fullpd->GetCellData()->GetArray("FreeEdge")->GetRange(fullfreeedge);

  std::cout<<"FULL SURFACE BAD TRI MIN: "<<fullbadtri[0]<<" MAX: "<<fullbadtri[1]<<endl;
  std::cout<<"FULL SURFACE FREE EDGE MIN: "<<fullfreeedge[0]<<" MAX: "<<fullfreeedge[1]<<endl;

  if (op == 0)
    WriteVTPFile(inputFilename1,fullpd,"_"+getRawName(inputFilename2)+"_Union");
  else if (op == 1)
    WriteVTPFile(inputFilename1,fullpd,"_"+getRawName(inputFilename2)+"_Intersection");
  else if (op == 2)
    WriteVTPFile(inputFilename1,fullpd,"_"+getRawName(inputFilename2)+"_Difference");

  WriteVTPFile(inputFilename1,myBoolean->GetOutput(1),"_"+getRawName(inputFilename2)+"_IntersectionLines");
  //Exit the program without errors
  return EXIT_SUCCESS;
}




