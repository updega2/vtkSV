//
//  TestControlGrid.cxx
//
//
//  Created by Adam Updegrove on 9/14/16.
//
//

/*=========================================================================

 Program:   SimVascular

 =========================================================================*/

#include "vtkControlGrid.h"
#include "vtkDataObject.h"
#include "vtkLoftNURBSCurve.h"
#include "vtkNURBSCurve.h"
#include "vtkNURBSUtils.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"
#include "vtkSTLReader.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLStructuredGridWriter.h"

#include <string>
#include <sstream>
#include <iostream>

#include <unistd.h>

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

void WriteVTSFile(std::string inputFilename,vtkStructuredGrid *writeStructuredGrid,std::string attachName)
{
  std::string rawName, pathName, outputFilename;

  vtkSmartPointer<vtkXMLStructuredGridWriter> writer  = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();

  pathName = getPath(inputFilename);
  rawName = getRawName(inputFilename);

  outputFilename = pathName+"/"+rawName+attachName+".vts";

  writer->SetFileName(outputFilename.c_str());
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(writeStructuredGrid);
#else
  writer->SetInputData(writeStructuredGrid);
#endif
  //writer->SetDataModeToAscii();

  writer->Write();
}

int main(int argc, char *argv[])
{
  if (argc != 1)
  {
      std::cout << "Incorrect Usage! Should be:" <<endl;
      std::cout << "./TestNURBSCurve" <<endl;
      return EXIT_FAILURE;
  }

  int nc = 11;

  vtkSmartPointer<vtkDoubleArray> x_data =
    vtkSmartPointer<vtkDoubleArray>::New();
  x_data->SetNumberOfTuples(nc);
  vtkNURBSUtils::LinSpace(0, 2*M_PI, nc, x_data);

  vtkSmartPointer<vtkPoints> inputPoints =
    vtkSmartPointer<vtkPoints>::New();
  inputPoints->Reset();
  for (int i=0; i<nc; i++)
  {
    double xval = x_data->GetTuple1(i);
    double yval = sin(xval - M_PI/2.0);
    double zval = 0.0;

    inputPoints->InsertNextPoint(xval, yval, zval);
  }

  vtkSmartPointer<vtkPolyData> inputPoly =
    vtkSmartPointer<vtkPolyData>::New();
  inputPoly->SetPoints(inputPoints);
  vtkSmartPointer<vtkLoftNURBSCurve> lofter =
    vtkSmartPointer<vtkLoftNURBSCurve>::New();
  lofter->SetInputData(inputPoly);
  lofter->SetDegree(2);
  lofter->SetPolyDataSpacing(0.01);
  lofter->SetKnotSpanType("average");
  lofter->Update();

  std::string newDirName = getcwd(NULL, 0);
  system(("mkdir -p "+newDirName+"/../Tests").c_str());
  //Write Files
  std::cout<<"Done...Writing Files..."<<endl;
  WriteVTPFile(newDirName+"/../Tests/Loft_Curve.vtp",lofter->GetOutput(0),"");
  WriteVTSFile(newDirName+"/../Tests/Loft_CurveControlPoints.vts",lofter->GetCurve()->GetControlPointGrid(),"");
  std::cout<<"Done"<<endl;

  return EXIT_SUCCESS;
}
