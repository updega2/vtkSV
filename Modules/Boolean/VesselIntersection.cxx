//
//  VesselIntersection.cxx
//
//
//  Created by Adam Updegrove on 1/3/15.
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
#include "vtkBooleanOperationPolyDataFilter2.h"
#include "vtkIntersectionPolyDataFilter.h"
#include "vtkBooleanOperationPolyDataFilter.h"
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
#include "vtkMultiplePolyDataIntersectionFilter.h"

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

//Function to read in the STL file, extract the boundaries and pass the input
//Poly Data information
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
  std::string *inputFilenames;
  inputFilenames = new std::string[argc-1];
  vtkPolyData **inputPDs;
  inputPDs = new vtkPolyData*[argc-1];
  //Create string from input File Name on command line
  for (int i =0;i<argc-1;i++)
  {
    inputFilenames[i] = argv[i+1];
    inputPDs[i] = vtkPolyData::New();
    ReadVTPFile(inputFilenames[i],inputPDs[i]);
  }

  vtkSmartPointer<vtkPolyData> fullpd =
	  vtkSmartPointer<vtkPolyData>::New();
  //FULL COMPLETE BOOLEAN EMBEDDED BOOLEAN FILTERS
  vtkSmartPointer<vtkMultiplePolyDataIntersectionFilter> vesselInter =
    vtkSmartPointer<vtkMultiplePolyDataIntersectionFilter>::New();
  for (int i = 0; i< argc-1; i++)
  {
    vesselInter->AddInputData(inputPDs[i]);
  }
  vesselInter->SetPassInfoAsGlobal(1);
  vesselInter->SetAssignSurfaceIds(1);
  vesselInter->SetNoIntersectionOutput(1);
  vesselInter->Update();

  fullpd->DeepCopy(vesselInter->GetOutput());
  double dummy[2];
  vtkIntersectionPolyDataFilter2::CleanAndCheckSurface(fullpd,dummy,1e-6);
  double fullbadtri[2], fullfreeedge[2];
  fullpd->GetCellData()->GetArray("BadTri")->GetRange(fullbadtri,0);
  fullpd->GetCellData()->GetArray("FreeEdge")->GetRange(fullfreeedge,0);

  std::cout<<"FULL SURFACE BAD TRI MIN: "<<fullbadtri[0]<<" MAX: "<<fullbadtri[1]<<endl;
  std::cout<<"FULL SURFACE FREE EDGE MIN: "<<fullfreeedge[0]<<" MAX: "<<fullfreeedge[1]<<endl;


  WriteVTPFile(inputFilenames[0],fullpd,"_FullBoolean");

  //Delete memory used
  delete [] inputFilenames;
  for (int i =0;i<argc-1;i++)
  {
    inputPDs[i]->Delete();
  }
  delete [] inputPDs;
  //Exit the program without errors
  return EXIT_SUCCESS;
}




