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

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkSVFindGeodesicPath.h"
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

  // argc is the number of strings on the command-line
  //  starting with the program name
  for(iarg=1; iarg<argc; iarg++){
      arglength = strlen(argv[iarg]);
      // replace 0..arglength-1 with argv[iarg]
      tmpstr.replace(0,arglength,argv[iarg],0,arglength);
      if(tmpstr=="-h")                      {RequestedHelp = true;}
      else if(tmpstr=="-input")             {InputProvided = true; inputFilename = argv[++iarg];}
      else if(tmpstr=="-output")            {OutputProvided = true; outputFilename = argv[++iarg];}
      else {cout << argv[iarg] << " is not a valid argument. Ask for help with -h." << endl; RequestedHelp = true; return EXIT_FAILURE;}
      // reset tmpstr for next argument
      tmpstr.erase(0,arglength);
  }

  if (RequestedHelp || !InputProvided)
  {
    cout << endl;
    cout << "usage:" <<endl;
    cout << "  RawToVtk -input [Input Filename] -output [Output Filename]..." << endl;
    cout << endl;
    cout << "COMMAND-LINE ARGUMENT SUMMARY" << endl;
    cout << "  -h                  : Display usage and command-line argument summary"<< endl;
    cout << "  -input              : Input file name (.vtp)"<< endl;
    cout << "  -output             : Output file name (.raw)"<< endl;
    cout << "END COMMAND-LINE ARGUMENT SUMMARY" << endl;
    return EXIT_FAILURE;
  }
  if (!OutputProvided)
  {
    std::cout << "WARNING: Output Filename not provided, setting output name based on the input filename" <<endl;
    outputFilename = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+".raw";
  }

  if (strncmp(vtkSVIOUtils::GetExt(inputFilename).c_str(), "vtp", 3))
  {
    std::cout << "Error, input does not have vtp extension" << endl;
    return EXIT_FAILURE;
  }

  if (strncmp(vtkSVIOUtils::GetExt(outputFilename).c_str(), "raw", 3))
  {
    std::cout << "Error, output does not have raw extension" << endl;
    return EXIT_FAILURE;
  }

  // Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkNew(vtkPolyData, inputPd);
  if (vtkSVIOUtils::ReadVTPFile(inputFilename,inputPd) != SV_OK)
    return EXIT_FAILURE;


  //Write Files
  std::cout<<"Writing Files..."<<endl;
  if (vtkSVIOUtils::WriteRawFile(outputFilename, inputPd) != SV_OK)
    return EXIT_FAILURE;

  //Exit the program without errors
  return EXIT_SUCCESS;
}
