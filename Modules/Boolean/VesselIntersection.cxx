/*=========================================================================
 *
 * Copyright (c) 2014 The Regents of the University of California.
 * All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *=========================================================================*/

/**
 *  \file VesselIntersection.cxx
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */

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




