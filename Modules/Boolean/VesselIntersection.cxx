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
#include "vtkSVLoopBooleanPolyDataFilter.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkSVLoopIntersectionPolyDataFilter.h"
#include "vtkSVMultiplePolyDataIntersectionFilter.h"

#include <string>
#include <sstream>
#include <iostream>

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
    vtkSVIOUtils::ReadVTPFile(inputFilenames[i],inputPDs[i]);
  }

  vtkNew(vtkPolyData, fullpd);
  //FULL COMPLETE BOOLEAN EMBEDDED BOOLEAN FILTERS
  vtkNew(vtkSVMultiplePolyDataIntersectionFilter, vesselInter);
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
  vtkSVLoopIntersectionPolyDataFilter::CleanAndCheckSurface(fullpd,dummy,1e-6);
  double fullbadtri[2], fullfreeedge[2];
  fullpd->GetCellData()->GetArray("BadTri")->GetRange(fullbadtri,0);
  fullpd->GetCellData()->GetArray("FreeEdge")->GetRange(fullfreeedge,0);

  std::cout<<"FULL SURFACE BAD TRI MIN: "<<fullbadtri[0]<<" MAX: "<<fullbadtri[1]<<endl;
  std::cout<<"FULL SURFACE FREE EDGE MIN: "<<fullfreeedge[0]<<" MAX: "<<fullfreeedge[1]<<endl;


  vtkSVIOUtils::WriteVTPFile(inputFilenames[0],fullpd,"_FullBoolean");

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




