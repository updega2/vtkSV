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
#include "vtkTriangleFilter.h"

#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkSVCenterlines.h"
#include "vtkSVCenterlineBranchSplitter.h"
#include "vtkSVSeedSelector.h"
#include "vtkSVPickPointSeedSelector.h"
#include "vtkSVOpenProfilesSeedSelector.h"

#include "vtkvmtkMergeCenterlines.h"
#include "vtkvmtkCapPolyData.h"
#include "vtkvmtkPolyDataCenterlines.h"
#include "vtkvmtkCenterlineBranchExtractor.h"

int main(int argc, char *argv[])
{
  // BEGIN PROCESSING COMMAND-LINE ARGUMENTS
  // Assume no options specified at command line
  bool RequestedHelp       = false;
  bool InputProvided       = false;
  bool OutputProvided      = false;

  // Variables used in processing the commandline
  int iarg, arglength;
  std::string tmpstr;

  // Filenames
  std::string inputFilename;
  std::string outputFilename;

  // Default values for options
  int useVmtk = 0;
  int useAbsoluteMergeDistance = 0;
  double mergeDistance = 0.1;
  double radiusMergeRatio = 0.5;
  std::string radiusArrayName   = "MaximumInscribedSphereRadius";

  // argc is the number of strings on the command-line
  //  starting with the program name
  for(iarg=1; iarg<argc; iarg++){
      arglength = strlen(argv[iarg]);
      // replace 0..arglength-1 with argv[iarg]
      tmpstr.replace(0,arglength,argv[iarg],0,arglength);
      if(tmpstr=="-h")                     {RequestedHelp = true;}
      else if(tmpstr=="-input")            {InputProvided = true; inputFilename = argv[++iarg];}
      else if(tmpstr=="-output")           {OutputProvided = true; outputFilename = argv[++iarg];}
      else if(tmpstr=="-radius")           {radiusArrayName = argv[++iarg];}
      else if(tmpstr=="-usevmtk")          {useVmtk = atoi(argv[++iarg]);}
      else if(tmpstr=="-radiusmergeratio") {radiusMergeRatio = atof(argv[++iarg]);}
      else if(tmpstr=="-usemergedistance") {useAbsoluteMergeDistance = atoi(argv[++iarg]);}
      else if(tmpstr=="-mergedistance")    {mergeDistance = atof(argv[++iarg]);}
      else {cout << argv[iarg] << " is not a valid argument. Ask for help with -h." << endl; RequestedHelp = true; return EXIT_FAILURE;}
      // reset tmpstr for next argument
      tmpstr.erase(0,arglength);
  }

  if (RequestedHelp || !InputProvided)
  {
    cout << endl;
    cout << "usage:" <<endl;
    cout << "  CenterlinesSplitter -input [Input Filename] -output [Output Filename] ..." << endl;
    cout << endl;
    cout << "COMMAND-LINE ARGUMENT SUMMARY" << endl;
    cout << "  -h                  : Display usage and command-line argument summary"<< endl;
    cout << "  -input              : Input file name (.vtp or .stl)"<< endl;
    cout << "  -output             : Output file name"<< endl;
    cout << "  -radius             : Name on centerlines describing maximum inscribed sphere radius [default MaximumInscribedSphereRadius]"<< endl;
    cout << "  -usevmtk            : Use the vmtk centerlines extractor rather than vtksv [default 0]"<< endl;
    cout << "  -radiusmergeratio   : When extracting centerline branches, the portion of the radius to use (radius at bifurcation location) to use as the merging distance [default 0.35]"<< endl;
    cout << "  -usemergedistance   : Instead of using a ratio to the radius, use an absolute distance for the merge distance [default 0]" << endl;
    cout << "  -mergedistance      : The merge distance; only used is usemergedistance is on [default 0.1]" << endl;
    cout << "END COMMAND-LINE ARGUMENT SUMMARY" << endl;
    return EXIT_FAILURE;
  }
  if (!OutputProvided)
  {
    cout << "WARNING: Output Filename not provided, setting output name based on the input filename" <<endl;
    std::string newDirName = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename);
    // Only mac and linux!!!
    system(("mkdir -p "+newDirName).c_str());
    outputFilename = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"_Split.vtp";
  }

  // Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkNew(vtkPolyData, inputPd);
  if (vtkSVIOUtils::ReadInputFile(inputFilename,inputPd) != 1)
    return EXIT_FAILURE;

  std::cout<<"Splitting Centerlines..."<<endl;

  if (useVmtk)
  {
    vtkNew(vtkvmtkCenterlineBranchExtractor, BranchSplitter);
    BranchSplitter->SetInputData(inputPd);
    //BranchSplitter->SetGroupingModeToFirstPoint();
    BranchSplitter->SetBlankingArrayName("Blanking");
    BranchSplitter->SetRadiusArrayName(radiusArrayName.c_str());
    BranchSplitter->SetGroupIdsArrayName("GroupIds");
    BranchSplitter->SetCenterlineIdsArrayName("CenterlineIds");
    BranchSplitter->SetTractIdsArrayName("TractIds");
    BranchSplitter->Update();

    std::cout<<"Done"<<endl;

    //Write Files
    std::cout<<"Writing Files..."<<endl;
    vtkSVIOUtils::WriteVTPFile(outputFilename, BranchSplitter->GetOutput(0));
  }
  else
  {
    vtkNew(vtkSVCenterlineBranchSplitter, BranchSplitter);
    //vtkNew(vtkvmtkCenterlineBranchExtractor, BranchSplitter);
    BranchSplitter->SetInputData(inputPd);
    BranchSplitter->SetGroupingModeToFirstPoint();
    BranchSplitter->SetBlankingArrayName("Blanking");
    BranchSplitter->SetRadiusArrayName(radiusArrayName.c_str());
    BranchSplitter->SetGroupIdsArrayName("GroupIds");
    BranchSplitter->SetCenterlineIdsArrayName("CenterlineIds");
    BranchSplitter->SetTractIdsArrayName("TractIds");
    BranchSplitter->SetRadiusMergeRatio(radiusMergeRatio);
    BranchSplitter->SetUseAbsoluteMergeDistance(useAbsoluteMergeDistance);
    BranchSplitter->SetMergeDistance(mergeDistance);
    BranchSplitter->Update();

    std::cout<<"Done"<<endl;

    //Write Files
    std::cout<<"Writing Files..."<<endl;
    vtkSVIOUtils::WriteVTPFile(outputFilename, BranchSplitter->GetOutput(0));
  }

  //Exit the program without errors
  return EXIT_SUCCESS;
}
