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
#include "vtkSVSurfaceCenterlineGrouper.h"
#include "vtkSVPassDataArray.h"

#include "vtkvmtkPolyDataCenterlineGroupsClipper.h"

int main(int argc, char *argv[])
{
  // BEGIN PROCESSING COMMAND-LINE ARGUMENTS
  // Assume no options specified at command line
  bool RequestedHelp       = false;
  bool InputProvided       = false;
  bool OutputProvided      = false;
  bool CenterlinesProvided = false;
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
  int useVmtk = 0;
  int useRadiusInfo = 1;
  int groupSurface = 1;
  int enforcePolycubeConnectivity = 0;
  int enforceCenterlinesConnectivity = 0;

  std::string groupIdsArrayName = "GroupIds";
  std::string radiusArrayName   = "MaximumInscribedSphereRadius";
  std::string blankingArrayName = "Blanking";
  std::string patchIdsArrayName = "PatchIds";
  std::string slicePointsArrayName = "SlicePoints";

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
      else if(tmpstr=="-usevmtk")       {useVmtk = atoi(argv[++iarg]);}
      else if(tmpstr=="-polycube")      {PolycubeProvided = true; polycubeFilename = argv[++iarg];}
      else if(tmpstr=="-groupids")      {groupIdsArrayName = argv[++iarg];}
      else if(tmpstr=="-radius")        {radiusArrayName = argv[++iarg];}
      else if(tmpstr=="-blanking")      {blankingArrayName = argv[++iarg];}
      else if(tmpstr=="-patchids")      {patchIdsArrayName = argv[++iarg];}
      else if(tmpstr=="-slicepoints")   {slicePointsArrayName = argv[++iarg];}
      else if(tmpstr=="-useradiusinfo") {useRadiusInfo = atoi(argv[++iarg]);}
      else if(tmpstr=="-groupsurface")  {groupSurface = atoi(argv[++iarg]);}
      else if(tmpstr=="-enforcecenterlinesconnectivity")     {enforceCenterlinesConnectivity = atoi(argv[++iarg]);}
      else if(tmpstr=="-enforcepolycubeconnectivity")     {enforcePolycubeConnectivity = atoi(argv[++iarg]);}
      else {cout << argv[iarg] << " is not a valid argument. Ask for help with -h." << endl; RequestedHelp = true; return EXIT_FAILURE;}
      // reset tmpstr for next argument
      tmpstr.erase(0,arglength);
  }

  if (RequestedHelp || !InputProvided || !CenterlinesProvided)
  {
    cout << endl;
    cout << "usage:" <<endl;
    cout << "  SurfaceCenterlineGrouper -input [Input Filename] -mergedcenterlines [Centerlines] -output [Output Filename] -groupids [GroupIds Array Name] ..." << endl;
    cout << endl;
    cout << "COMMAND-LINE ARGUMENT SUMMARY" << endl;
    cout << "  -h                              : Display usage and command-line argument summary"<< endl;
    cout << "  -input                          : Input file name (.vtp or .stl)"<< endl;
    cout << "  -mergedcenterlines              : Split and merged centerlines file name (.vtp)"<< endl;
    cout << "  -output                         : Output file name"<< endl;
    cout << "  -usevmtk                        : Use the vmtk clipper rather than vtksv [default 0]"<< endl;
    cout << "  -polycube                       : Polycube file name (.vtp)"<< endl;
    cout << "  -groupids                       : Name to be used for group ids [default GroupIds]"<< endl;
    cout << "  -radius                         : Name on centerlines describing maximum inscribed sphere radius [default MaximumInscribedSphereRadius]"<< endl;
    cout << "  -blanking                       : Name on centerlines describing whether line is part of bifurcation region or not [default Blanking]"<< endl;
    cout << "  -patchids                       : Name on polycube indicating the cuboid patch id [default PatchIds]"<< endl;
    cout << "  -slicepoints                    : Name to be placed on surface to indicate where patch divisions will occur if preparing for cuboid patching. Only used if enforcepolycubeconnectivity is turned on [default SlicePoints]"<< endl;
    cout << "  -useradiusinfo                  : Use radius to help in clipping operation [default 1]"<< endl;
    cout << "  -enforcecenterlinesconnectivity : Enforce the connectivity of the centerlines on the surface [default 0]" << endl;
    cout << "  -enforcepolycubeconnectivity    : Enforce the connectivity of the polycube on the surface [default 0]" << endl;
    cout << "  -groupsurface                   : Group the surface using centerlines. If surface already has groups and the polycube connectivity needs to be enforced, this can be turned off [default 1]" << endl;
    cout << "END COMMAND-LINE ARGUMENT SUMMARY" << endl;
    return EXIT_FAILURE;
  }
  if (!OutputProvided)
  {
    cout << "WARNING: Output Filename not provided, setting output name based on the input filename" <<endl;
    std::string newDirName = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename);
    // Only mac and linux!!!
    system(("mkdir -p "+newDirName).c_str());
    outputFilename = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"_Surface_Centerline_Grouped.vtp";
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
  if (useVmtk)
  {
    vtkNew(vtkSplineFilter, Resampler);
    Resampler->SetInputData(centerlinesPd);
    //Resampler->SetInputData(this->MergedCenterlines);
    Resampler->SetSubdivideToLength();
    Resampler->SetLength(centerlinesPd->GetLength()/100.);
    Resampler->Update();

    vtkNew(vtkvmtkPolyDataCenterlineGroupsClipper, BranchClipper);
    BranchClipper->SetInputData(inputPd);
    BranchClipper->SetCenterlines(Resampler->GetOutput());
    BranchClipper->SetGroupIdsArrayName(groupIdsArrayName.c_str());
    BranchClipper->SetCenterlineGroupIdsArrayName(groupIdsArrayName.c_str());
    BranchClipper->SetCenterlineRadiusArrayName(radiusArrayName.c_str());
    BranchClipper->SetBlankingArrayName(blankingArrayName.c_str());
    BranchClipper->SetCutoffRadiusFactor(VTK_SV_LARGE_DOUBLE);
    BranchClipper->SetClipValue(0.0);
    BranchClipper->SetUseRadiusInformation(useRadiusInfo);
    BranchClipper->SetClipAllCenterlineGroupIds(1);
    BranchClipper->Update();

    vtkNew(vtkSVPassDataArray, DataPasser);
    DataPasser->SetInputData(0, BranchClipper->GetOutput());
    DataPasser->SetInputData(1, inputPd);
    DataPasser->SetPassArrayName(groupIdsArrayName.c_str());
    DataPasser->SetPassDataIsCellData(0);
    DataPasser->SetPassDataToCellData(1);
    DataPasser->Update();

    //Write Files
    std::cout<<"Writing Files..."<<endl;
    vtkSVIOUtils::WriteVTPFile(outputFilename, DataPasser->GetOutput(0));

  }
  else
  {
    // Filter
    vtkNew(vtkSVSurfaceCenterlineGrouper, Grouper);

    Grouper->SetInputData(inputPd);
    Grouper->SetMergedCenterlines(centerlinesPd);
    if (PolycubeProvided)
    {
      Grouper->SetPolycubePd(polycubePd);
    }
    Grouper->SetCenterlineGroupIdsArrayName(groupIdsArrayName.c_str());
    Grouper->SetGroupIdsArrayName(groupIdsArrayName.c_str());
    Grouper->SetCenterlineRadiusArrayName(radiusArrayName.c_str());
    Grouper->SetBlankingArrayName(blankingArrayName.c_str());
    Grouper->SetCenterlineIdsArrayName("CenterlineIds");
    Grouper->SetTractIdsArrayName("TractIds");
    Grouper->SetPatchIdsArrayName(patchIdsArrayName.c_str());
    Grouper->SetSlicePointsArrayName(slicePointsArrayName.c_str());
    Grouper->SetUseRadiusInformation(useRadiusInfo);
    Grouper->SetEnforceCenterlinesConnectivity(enforceCenterlinesConnectivity);
    Grouper->SetEnforcePolycubeConnectivity(enforcePolycubeConnectivity);
    Grouper->SetGroupSurface(groupSurface);
    Grouper->DebugOn();
    Grouper->Update();
    std::cout<<"Done"<<endl;

    if (Grouper->GetErrorCode() != 0)
    {
      std::cerr << "Surface centerline grouper failed. " <<endl;
      return EXIT_FAILURE;
    }

    //Write Files
    std::cout<<"Writing Files..."<<endl;
    vtkSVIOUtils::WriteVTPFile(outputFilename, Grouper->GetOutput(0));

  }

  //Exit the program without errors
  return EXIT_SUCCESS;
}
