/*=========================================================================
 *
 * Copyright (c) 2014-2015 The Regents of the University of California.
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
#include "vtkSVPassDataArray.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"

#include <string>
#include <sstream>
#include <iostream>

int main(int argc, char *argv[])
{
  // BEGIN PROCESSING COMMAND-LINE ARGUMENTS
  // Assume the command line is okay
  bool BogusCmdLine = false;
  // Assume no options specified at command line
  bool RequestedHelp         = false;
  bool SourceProvided        = false;
  bool TargetProvided        = false;
  bool OutputProvided        = false;
  bool PassArrayNameProvided = false;

  // Variables used in processing the commandline
  int iarg, arglength;
  std::string tmpstr;

  // Filenames
  int passDataIsCellData = 1;
  int passDataToCellData = 1;
  std::string sourceFilename;
  std::string targetFilename;
  std::string outputFilename;

  // Default values for options
  std::string passArrayName = "PassArray";

  // argc is the number of strings on the command-line
  //  starting with the program name
  for(iarg=1; iarg<argc; iarg++){
      arglength = strlen(argv[iarg]);
      // replace 0..arglength-1 with argv[iarg]
      tmpstr.replace(0,arglength,argv[iarg],0,arglength);
      if(tmpstr=="-h")                {RequestedHelp = true;}
      else if(tmpstr=="-source")      {SourceProvided = true; sourceFilename = argv[++iarg];}
      else if(tmpstr=="-target")      {TargetProvided = true; targetFilename = argv[++iarg];}
      else if(tmpstr=="-output")      {OutputProvided = true; outputFilename = argv[++iarg];}
      else if(tmpstr=="-passarray")   {PassArrayNameProvided = true; passArrayName = argv[++iarg];}
      else if(tmpstr=="-iscelldata")  {passDataIsCellData = atoi(argv[++iarg]);}
      else if(tmpstr=="-passtocells") {passDataToCellData = atoi(argv[++iarg]);}
      else {BogusCmdLine = true;}
      // reset tmpstr for next argument
      tmpstr.erase(0,arglength);
  }

  if (RequestedHelp || !SourceProvided || !TargetProvided || !PassArrayNameProvided)
  {
    cout << endl;
    cout << "usage:" <<endl;
    cout << "  FindSeparateRegions -source [Source Filename] -target [Target Filename] -output [Output Filename] -passarray [Pass Array Name] ..." << endl;
    cout << endl;
    cout << "COMMAND-LINE ARGUMENT SUMMARY" << endl;
    cout << "  -h                  : Display usage and command-line argument summary"<< endl;
    cout << "  -source             : Source file name (.vtp or .stl)"<< endl;
    cout << "  -target             : Target file name (.vtp or .stl)"<< endl;
    cout << "  -ouptut             : Output file name"<< endl;
    cout << "  -passarray          : Name of data array that must be defined on the source surface [default PassArray]"<< endl;
    cout << "  -iscelldata         : Indicates whether the provided data array is cell data or point data on the source surface [default 1]"<< endl;
    cout << "  -passtocells        : Indicates whether the data should be passed to the cells or points of the target surface [default 1]"<< endl;
    cout << "END COMMAND-LINE ARGUMENT SUMMARY" << endl;
    return EXIT_FAILURE;
  }
  if (!OutputProvided)
  {
    cout << "WARNING: Output Filename not provided, setting output name based on the input filename" <<endl;
    std::string newDirName = vtkSVIOUtils::GetPath(targetFilename)+"/"+vtkSVIOUtils::GetRawName(targetFilename);
    // Only mac and linux!!!
    system(("mkdir -p "+newDirName).c_str());
    outputFilename = vtkSVIOUtils::GetPath(targetFilename)+"/"+vtkSVIOUtils::GetRawName(targetFilename)+"/"+vtkSVIOUtils::GetRawName(targetFilename)+"_With_"+vtkSVIOUtils::GetRawName(sourceFilename)+"_Data.vtp";
  }
  // Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkNew(vtkPolyData, sourcePd);
  vtkSVIOUtils::ReadInputFile(sourceFilename,sourcePd);
  vtkNew(vtkPolyData, targetPd);
  vtkSVIOUtils::ReadInputFile(targetFilename,targetPd);

  // Filter
  vtkNew(vtkSVPassDataArray, Passer);

  //OPERATION
  std::cout<<"Performing Operation..."<<endl;
  Passer->SetInputData(0, sourcePd);
  Passer->SetInputData(1, targetPd);
  Passer->SetPassArrayName(passArrayName.c_str());
  Passer->SetPassDataIsCellData(passDataIsCellData);
  Passer->SetPassDataToCellData(passDataToCellData);
  Passer->Update();
  std::cout<<"Done"<<endl;

  //Write Files
  std::cout<<"Writing Files..."<<endl;
  vtkSVIOUtils::WriteVTPFile(outputFilename, Passer->GetOutput(0));

  //Exit the program without errors
  return EXIT_SUCCESS;
}
