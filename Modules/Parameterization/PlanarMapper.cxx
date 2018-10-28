/* Copyright (c) Stanford University, The Regents of the University of
 *               California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
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
 */

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
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkSVPlanarMapper.h"

#include <string>
#include <sstream>
#include <iostream>

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
  vtkNew(vtkPolyData, pd1);
  vtkNew(vtkSVPlanarMapper, Mapper);

  //Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkSVIOUtils::ReadInputFile(inputFilename1,pd1);

  vtkNew(vtkIntArray, boundaryCorners);
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
  vtkNew(vtkSVSquareBoundaryMapper, boundaryMapper);
  boundaryMapper->SetBoundaryIds(boundaryCorners);

  std::string newDirName = vtkSVIOUtils::GetPath(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1);
  std::string newOutName = vtkSVIOUtils::GetPath(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1)+"/"+vtkSVIOUtils::GetRawName(inputFilename1);
  system(("mkdir -p "+newDirName).c_str());
  //OPERATION
  std::cout<<"Performing Operation..."<<endl;
  Mapper->SetInputData(pd1);
  Mapper->SetBoundaryMapper(boundaryMapper);
  Mapper->Update();

  //Write Files
  std::cout<<"Done...Writing Files..."<<endl;
  vtkSVIOUtils::WriteVTPFile(newOutName+".vtp",Mapper->GetOutput(0),"_Mapped");
  std::cout<<"Done"<<endl;

  //Exit the program without errors
  return EXIT_SUCCESS;
}
