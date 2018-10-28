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

/**
 *  \file PolyDataBoolean.cxx
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */

#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkSVLoopBooleanPolyDataFilter.h"
#include "vtkSVLoopIntersectionPolyDataFilter.h"

#include <string>
#include <sstream>
#include <iostream>

int main(int argc, char *argv[])
{
  // BEGIN PROCESSING COMMAND-LINE ARGUMENTS
  // Assume no options specified at command line
  bool RequestedHelp = false;
  bool Input0Provided = false;
  bool Input1Provided = false;
  bool OperationProvided = false;
  bool OutputProvided = false;

  // Variables used in processing the commandline
  int iarg, arglength;
  std::string tmpstr;

  // Filenames
  std::string input0Filename;
  std::string input1Filename;
  std::string outputFilename;

  // Default values for options
  double testTolerance = 1e-6;
  double boolOperation = 0;
  double writeIntersectionLines = 0;

  // argc is the number of strings on the command-line
  //  starting with the program name
  for(iarg=1; iarg<argc; iarg++){
      arglength = strlen(argv[iarg]);
      // replace 0..arglength-1 with argv[iarg]
      tmpstr.replace(0,arglength,argv[iarg],0,arglength);
      if(tmpstr=="-h")                      {RequestedHelp = true;}
      else if(tmpstr=="-input0")            {Input0Provided = true; input0Filename = argv[++iarg];}
      else if(tmpstr=="-input1")            {Input1Provided = true; input1Filename = argv[++iarg];}
      else if(tmpstr=="-output")            {OutputProvided = true; outputFilename = argv[++iarg];}
      else if(tmpstr=="-operation")         {OperationProvided = true; boolOperation = atoi(argv[++iarg]);}
      else if(tmpstr=="-tolerance")         {testTolerance = atof(argv[++iarg]);}
      else if(tmpstr=="-writeintersection") {writeIntersectionLines = atoi(argv[++iarg]);}
      else {cout << argv[iarg] << " is not a valid argument. Ask for help with -h." << endl; RequestedHelp = true; return EXIT_FAILURE;}
      // reset tmpstr for next argument
      tmpstr.erase(0,arglength);
  }
  if (RequestedHelp || !Input0Provided || !Input1Provided)
  {
    cout << endl;
    cout << "usage:" <<endl;
    cout << "  PolyDataBoolean -input0 [Input 0 Filename] -input1 [Input 1 Filename] -operation [Operation to perform] -output [Output Filename] ..." << endl;
    cout << endl;
    cout << "COMMAND-LINE ARGUMENT SUMMARY" << endl;
    cout << "  -h                  : Display usage and command-line argument summary"<< endl;
    cout << "  -input0             : First surface input file name (.vtp or .stl)"<< endl;
    cout << "  -input1             : Second surface input file name (.vtp or .stl)"<< endl;
    cout << "  -output             : Output file name"<< endl;
    cout << "  -operation          : Operation to perform: 0 - Union, 1 - Intersection, 2 - Difference."<< endl;
    cout << "  -tolerance          : Tolerance for the boolean operation [default 1.0e-6]"<< endl;
    cout << "  -writeintersection  : Write the intersection lines [default 0]"<< endl;
    cout << "END COMMAND-LINE ARGUMENT SUMMARY" << endl;
    return EXIT_FAILURE;
  }
  if (!OperationProvided)
  {
    cout << "ERROR: Need to provide an operation to perform: 0 - Union, 1 - Intersection, 2 - Difference." << endl;
    return EXIT_FAILURE;
  }
  std::string attachName;
  if (!OutputProvided)
  {
    cout << "WARNING: Output Filename not provided, setting output name based on the input filename" <<endl;
    std::string newDirName = vtkSVIOUtils::GetPath(input0Filename)+"/"+vtkSVIOUtils::GetRawName(input0Filename);
    // Only mac and linux!!!
    system(("mkdir -p "+newDirName).c_str());
    if (boolOperation == 0)
      attachName = "Union";
    if (boolOperation == 1)
      attachName = "Intersection";
    if (boolOperation == 2)
      attachName = "Difference";
    outputFilename = vtkSVIOUtils::GetPath(input0Filename)+"/"+vtkSVIOUtils::GetRawName(input0Filename)+"/"+vtkSVIOUtils::GetRawName(input0Filename)+"_"+vtkSVIOUtils::GetRawName(input1Filename)+attachName+".vtp";
  }

  // Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkNew(vtkPolyData, input0Pd);
  if (vtkSVIOUtils::ReadInputFile(input0Filename,input0Pd) != 1)
    return EXIT_FAILURE;
  vtkNew(vtkPolyData, input1Pd);
  if (vtkSVIOUtils::ReadInputFile(input1Filename,input1Pd) != 1)
    return EXIT_FAILURE;

  // Filter
  vtkNew(vtkSVLoopBooleanPolyDataFilter, myBoolean);

  //BOOLEAN OPERATION EMBEDDED INTERSECTION
  std::cout<<"Performing Operation..."<<endl;
  myBoolean->SetInputData(0, input0Pd);
  myBoolean->SetInputData(1, input1Pd);
  myBoolean->SetTolerance(testTolerance);
  myBoolean->SetOperation(boolOperation);
  myBoolean->Update();

  vtkNew(vtkPolyData, output);
  output->DeepCopy(myBoolean->GetOutput());
  std::cout<<"Checking Surface..."<<endl;

  output->GetCellData()->RemoveArray("BadTriangle");
  output->GetCellData()->RemoveArray("FreeEdge");

  double dummy[2];
  vtkSVLoopIntersectionPolyDataFilter::CleanAndCheckSurface(output,dummy,testTolerance);
  double fullbadtri[2], fullfreeedge[2];
  output->GetCellData()->GetArray("BadTriangle")->GetRange(fullbadtri);
  output->GetCellData()->GetArray("FreeEdge")->GetRange(fullfreeedge);

  std::cout<<"FULL SURFACE BAD TRI MIN: "<<fullbadtri[0]<<" MAX: "<<fullbadtri[1]<<endl;
  std::cout<<"FULL SURFACE FREE EDGE MIN: "<<fullfreeedge[0]<<" MAX: "<<fullfreeedge[1]<<endl;

  std::cout<<"Writing Files..."<<endl;
  vtkSVIOUtils::WriteVTPFile(outputFilename, output);
  if (writeIntersectionLines)
  {
    std::string interLinesSuffix = "_IntersectionLines";
    vtkSVIOUtils::WriteVTPFile(outputFilename, myBoolean->GetOutput(1), interLinesSuffix);
  }

  std::cout<<"Done"<<endl;
  //Exit the program without errors
  return EXIT_SUCCESS;
}




