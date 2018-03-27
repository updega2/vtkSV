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
#include "vtkInformation.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkSplineFilter.h"

#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkSVSurfaceCuboidPatcher.h"
#include "vtkSVPassDataArray.h"

int main(int argc, char *argv[])
{
  // BEGIN PROCESSING COMMAND-LINE ARGUMENTS
  // Assume no options specified at command line
  bool RequestedHelp       = false;
  bool InputProvided       = true;
  bool OutputProvided      = true;
  bool CenterlinesProvided = true;
  bool PolycubeProvided    = false;

  // Variables used in processing the commandline
  int iarg, arglength;
  std::string tmpstr;

  // Filenames
  std::string inputFilename;
  std::string outputFilename;
  std::string centerlinesFilename;
  std::string polycubeFilename;

  // Default values for options
  int useRadiusInfo = 1;

  std::string groupIdsArrayName = "GroupIds";
  std::string radiusArrayName   = "MaximumInscribedSphereRadius";
  std::string blankingArrayName = "Blanking";

  // argc is the number of strings on the command-line
  //  starting with the program name
  for(iarg=1; iarg<argc; iarg++){
      arglength = strlen(argv[iarg]);
      // replace 0..arglength-1 with argv[iarg]
      tmpstr.replace(0,arglength,argv[iarg],0,arglength);
      if(tmpstr=="-h")                  {RequestedHelp = true;}
      else if(tmpstr=="-input")         {InputProvided = true; inputFilename = argv[++iarg];}
      else if(tmpstr=="-mergedcenterlines")   {CenterlinesProvided = true; centerlinesFilename = argv[++iarg];}
      else if(tmpstr=="-output")        {OutputProvided = true; outputFilename = argv[++iarg];}
      else if(tmpstr=="-polycube")      {PolycubeProvided = true; polycubeFilename = argv[++iarg];}
      else if(tmpstr=="-groupids")      {groupIdsArrayName = argv[++iarg];}
      else if(tmpstr=="-radius")        {radiusArrayName = argv[++iarg];}
      else if(tmpstr=="-blanking")      {blankingArrayName = argv[++iarg];}
      else if(tmpstr=="-useradiusinfo") {useRadiusInfo = atoi(argv[++iarg]);}
      else {cout << argv[iarg] << " is not a valid argument. Ask for help with -h." << endl; RequestedHelp = true; return EXIT_FAILURE;}
      // reset tmpstr for next argument
      tmpstr.erase(0,arglength);
  }

  if (RequestedHelp || !InputProvided || !CenterlinesProvided)
  {
    cout << endl;
    cout << "usage:" <<endl;
    cout << "  SurfaceCuboidPatcher -input [Input Filename] -mergedcenterlines [Centerlines] -output [Output Filename] -groupids [GroupIds Array Name] ..." << endl;
    cout << endl;
    cout << "COMMAND-LINE ARGUMENT SUMMARY" << endl;
    cout << "  -h                              : Display usage and command-line argument summary"<< endl;
    cout << "  -input                          : Input file name (.vtp or .stl)"<< endl;
    cout << "  -mergedcenterlines              : Split and merged centerlines file name (.vtp)"<< endl;
    cout << "  -output                         : Output file name"<< endl;
    cout << "  -polycube                       : Polycube file name (.vtp)"<< endl;
    cout << "  -groupids                       : Name to be used for group ids [default GroupIds]"<< endl;
    cout << "  -radius                         : Name on centerlines describing maximum inscribed sphere radius [default MaximumInscribedSphereRadius]"<< endl;
    cout << "  -blanking                       : Name on centerlines describing whether line is part of bifurcation region or not [default Blanking]"<< endl;
    cout << "  -useradiusinfo                  : Use radius to help in clipping operation [default 1]"<< endl;
    cout << "END COMMAND-LINE ARGUMENT SUMMARY" << endl;
    return EXIT_FAILURE;
  }
  if (!OutputProvided)
  {
    cout << "WARNING: Output Filename not provided, setting output name based on the input filename" <<endl;
    std::string newDirName = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename);
    // Only mac and linux!!!
    system(("mkdir -p "+newDirName).c_str());
    outputFilename = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"_Cuboid_Patched.vtp";
  }

  // Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkNew(vtkPolyData, inputPd);
  if (vtkSVIOUtils::ReadInputFile(inputFilename,inputPd) != 1)
    return EXIT_FAILURE;
  vtkNew(vtkPolyData, centerlinesPd);
  if (vtkSVIOUtils::ReadInputFile(centerlinesFilename,centerlinesPd) != 1)
    return EXIT_FAILURE;

  vtkNew(vtkPolyData, polycubePd);
  if (PolycubeProvided)
  {
    if (vtkSVIOUtils::ReadInputFile(polycubeFilename, polycubePd) != 1)
      return EXIT_FAILURE;
  }

  //OPERATION
  std::cout<<"Performing Operation..."<<endl;
  // Filter
  vtkNew(vtkSVSurfaceCuboidPatcher, Patcher);

  Patcher->SetInputData(inputPd);
  Patcher->SetMergedCenterlines(centerlinesPd);
  if (PolycubeProvided)
  {
    Patcher->SetPolycubePd(polycubePd);
  }
  Patcher->SetCenterlineGroupIdsArrayName(groupIdsArrayName.c_str());
  Patcher->SetGroupIdsArrayName(groupIdsArrayName.c_str());
  Patcher->SetCenterlineRadiusArrayName(radiusArrayName.c_str());
  Patcher->SetBlankingArrayName(blankingArrayName.c_str());
  Patcher->SetCenterlineIdsArrayName("CenterlineIds");
  Patcher->SetTractIdsArrayName("TractIds");
  Patcher->DebugOn();
  Patcher->Update();
  std::cout<<"Done"<<endl;

  //Write Files
  std::cout<<"Writing Files..."<<endl;
  vtkSVIOUtils::WriteVTPFile(outputFilename, Patcher->GetOutput(0));

  if (Patcher->GetErrorCode() != 0)
  {
    std::cerr << "Surface centerline grouper failed. " <<endl;
    return EXIT_FAILURE;
  }

  //Exit the program without errors
  return EXIT_SUCCESS;
}
