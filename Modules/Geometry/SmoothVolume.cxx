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
#include "vtkSVSmoothVolume.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkTriangle.h"

#include <string>
#include <sstream>
#include <iostream>

int main(int argc, char *argv[])
{
  // BEGIN PROCESSING COMMAND-LINE ARGUMENTS
  // Assume no options specified at command line
  bool RequestedHelp     = false;
  bool InputProvided     = false;
  bool OutputProvided    = false;

  // Variables used in processing the commandline
  int iarg, arglength;
  std::string tmpstr;

  // Filenames
  std::string inputFilename;
  std::string outputFilename;

  // Default values for options
  int numSmoothIters = 5;

  // argc is the number of strings on the command-line
  //  starting with the program name
  for(iarg=1; iarg<argc; iarg++){
      arglength = strlen(argv[iarg]);
      // replace 0..arglength-1 with argv[iarg]
      tmpstr.replace(0,arglength,argv[iarg],0,arglength);
      if(tmpstr=="-h")                      {RequestedHelp = true;}
      else if(tmpstr=="-input")             {InputProvided = true; inputFilename = argv[++iarg];}
      else if(tmpstr=="-output")            {OutputProvided = true; outputFilename = argv[++iarg];}
      else if(tmpstr=="-numsmoothiters")    {numSmoothIters = atoi(argv[++iarg]);}
      else {cout << argv[iarg] << " is not a valid argument. Ask for help with -h." << endl; RequestedHelp = true; return EXIT_FAILURE;}
      // reset tmpstr for next argument
      tmpstr.erase(0,arglength);
  }

  if (RequestedHelp || !InputProvided)
  {
    cout << endl;
    cout << "usage:" <<endl;
    cout << "  Smoother -input [Input Filename (vtu)] -output [Output Filename] -numsmoothiters [Smooth Iterations] ..." << endl;
    cout << endl;
    cout << "COMMAND-LINE ARGUMENT SUMMARY" << endl;
    cout << "  -h              : Display usage and command-line argument summary"<< endl;
    cout << "  -input          : Input file name (.vtu)"<< endl;
    cout << "  -output         : Output file name"<< endl;
    cout << "  -numsmoothiters : Point Id to start from [default 5]"<< endl;
    cout << "END COMMAND-LINE ARGUMENT SUMMARY" << endl;
    return EXIT_FAILURE;
  }
  if (!OutputProvided)
  {
    cout << "WARNING: Output Filename not provided, setting output name based on the input filename" <<endl;
    std::string newDirName = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename);
    // Only mac and linux!!!
    system(("mkdir -p "+newDirName).c_str());
    outputFilename = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"_Smoothed.vtu";
  }

  // Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkNew(vtkUnstructuredGrid, inputUg);
  if (vtkSVIOUtils::ReadVTUFile(inputFilename,inputUg) != 1)
    return EXIT_FAILURE;

  // Filter
  vtkNew(vtkSVSmoothVolume, Smoother);

  //OPERATION
  std::cout<<"Performing Operation..."<<endl;
  Smoother->SetInputData(0, inputUg);
  Smoother->SetNumberOfSmoothIterations(numSmoothIters);
  Smoother->Update();
  std::cout<<"Done"<<endl;

  if (Smoother->GetErrorCode() != 0)
  {
    std::cout << "Error in running filter" << endl;
    return SV_ERROR;
  }

  //Write Files
  std::cout<<"Writing Files..."<<endl;
  vtkSVIOUtils::WriteVTUFile(outputFilename, Smoother->GetOutput(0));

  //Exit the program without errors
  return EXIT_SUCCESS;
}
