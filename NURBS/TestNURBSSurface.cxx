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
#include "vtkLoftNURBSSurface.h"
#include "vtkNURBSSurface.h"
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

  int nU = 7;
  int nV = 5;
  vtkPolyData **inputPolys = new vtkPolyData*[nU];
  for (int i=0; i<nU; i++)
  {
    inputPolys[i] = vtkPolyData::New();
    vtkSmartPointer<vtkPoints> newPoints =
      vtkSmartPointer<vtkPoints>::New();
    inputPolys[i]->SetPoints(newPoints);
    inputPolys[i]->GetPoints()->SetNumberOfPoints(nV);
  }

  vtkSmartPointer<vtkDoubleArray> xdata =
    vtkSmartPointer<vtkDoubleArray>::New();
  xdata->SetNumberOfTuples(nU);
  vtkSmartPointer<vtkDoubleArray> ydata =
    vtkSmartPointer<vtkDoubleArray>::New();
  ydata->SetNumberOfTuples(nV);
  vtkNURBSUtils::LinSpace(0, 2*M_PI, nU, xdata);
  vtkNURBSUtils::LinSpace(0, 2*M_PI, nV, ydata);

  vtkSmartPointer<vtkDoubleArray> uders =
    vtkSmartPointer<vtkDoubleArray>::New();
  uders->SetNumberOfComponents(3);
  uders->SetNumberOfTuples(nU);
  double uder[3]; uder[0] = 2.0; uder[1] = 0.0; uder[2] = 0.0;
  for (int i=0; i<nU; i++)
  {
    uders->SetTuple(i, uder);
  }

  vtkSmartPointer<vtkDoubleArray> vders =
    vtkSmartPointer<vtkDoubleArray>::New();
  vders->SetNumberOfComponents(3);
  vders->SetNumberOfTuples(nV);
  double vder[3]; vder[0] = 0.0; vder[1] = 2.0; vder[2] = 0.0;
  for (int i=0; i<nV; i++)
  {
    vders->SetTuple(i, vder);
  }

  for (int i=0; i<nU; i++)
  {
    for (int j=0; j<nV; j++)
    {
      double xval = xdata->GetTuple1(i);
      double yval = ydata->GetTuple1(j);;
      double zval = sin(xval) * sin(yval);

      inputPolys[i]->GetPoints()->SetPoint(j, xval, yval, zval);
    }
  }

  vtkSmartPointer<vtkLoftNURBSSurface> lofter =
    vtkSmartPointer<vtkLoftNURBSSurface>::New();
  for (int i=0; i<nU; i++)
  {
    lofter->AddInputData(inputPolys[i]);
  }
  lofter->SetUDegree(2);
  lofter->SetVDegree(2);
  lofter->SetPolyDataUSpacing(0.01);
  lofter->SetPolyDataVSpacing(0.01);
  //lofter->SetUKnotSpanType("derivative");
  //lofter->SetStartUDerivatives(uders);
  //lofter->SetEndUDerivatives(uders);
  //lofter->SetVKnotSpanType("derivative");
  //lofter->SetStartVDerivatives(vders);
  //lofter->SetEndVDerivatives(vders);
  lofter->Update();

  for (int i=0; i<nU; i++)
  {
    inputPolys[i]->Delete();
  }
  delete [] inputPolys;

  std::string newDirName = getcwd(NULL, 0);
  system(("mkdir -p "+newDirName+"/../Tests").c_str());
  //Write Files
  std::cout<<"Done...Writing Files..."<<endl;
  WriteVTPFile(newDirName+"/../Tests/Loft_Surface.vtp",lofter->GetOutput(0),"");
  WriteVTSFile(newDirName+"/../Tests/Loft_SurfaceControlPoints.vts",lofter->GetSurface()->GetControlPointGrid(),"");
  std::cout<<"Done"<<endl;

  return EXIT_SUCCESS;
}
