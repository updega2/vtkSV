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
#include "vtkSVSurfaceCenterlineAttributesPasser.h"
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
  int isVasculature = 1;
  int enforceBoundaryDirections = 0;
  int boundaryEnforceFactor = 1;

  double normalsWeighting = 0.8;

  std::string groupIdsArrayName = "GroupIds";
  std::string radiusArrayName   = "MaximumInscribedSphereRadius";
  std::string blankingArrayName = "Blanking";
  std::string patchIdsArrayName = "PatchIds";
  std::string slicePointsArrayName = "SlicePoints";
  std::string parallelTransportVectorArrayName = "ParallelTransportVector";
  std::string transformedNormalsArrayName = "NormalsTransformedToCenterlines";

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
      else if(tmpstr=="-patchids")      {patchIdsArrayName = argv[++iarg];}
      else if(tmpstr=="-slicepoints")   {slicePointsArrayName = argv[++iarg];}
      else if(tmpstr=="-useradiusinfo") {useRadiusInfo = atoi(argv[++iarg]);}
      else if(tmpstr=="-isvasculature")                 {isVasculature = atoi(argv[++iarg]);}
      else if(tmpstr=="-paralleltransportvector")       {parallelTransportVectorArrayName = argv[++iarg];}
      else if(tmpstr=="-transformednormals")            {transformedNormalsArrayName = argv[++iarg];}
      else if(tmpstr=="-enforceboundarydirections")     {enforceBoundaryDirections = atoi(argv[++iarg]);}
      else if(tmpstr=="-boundaryenforcefactor")         {boundaryEnforceFactor = atoi(argv[++iarg]);}
      else if(tmpstr=="-normalsweighting")              {normalsWeighting = atof(argv[++iarg]);}
      else {cout << argv[iarg] << " is not a valid argument. Ask for help with -h." << endl; RequestedHelp = true; return EXIT_FAILURE;}
      // reset tmpstr for next argument
      tmpstr.erase(0,arglength);
  }

  if (RequestedHelp || !InputProvided || !CenterlinesProvided)
  {
    cout << endl;
    cout << "usage:" <<endl;
    cout << "  SurfaceCenterlineAttributesPasser -input [Input Filename] -mergedcenterlines [Centerlines] -output [Output Filename] -groupids [GroupIds Array Name] ..." << endl;
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
    cout << "  -patchids                       : Name on polycube indicating the cuboid patch id [default PatchIds]"<< endl;
    cout << "  -slicepoints                    : Name to be placed on surface to indicate where patch divisions will occur if preparing for cuboid patching. Only used if enforcepolycubeconnectivity is turned on [default SlicePoints]"<< endl;
    cout << "  -paralleltransportvector        : Name on centerlines for centerlines local coordinate frame. If not provided, the coordinate frame will be computed and given this name. [default ParallelTransportVector]"<< endl;
    cout << "  -transformednormals            : Name to be placed on surface after the normals are transformed to the local coordinate system on the centerlines. [default NormalsTransformedToCenterlines]"<< endl;
    cout << "  -useradiusinfo                  : Use radius to help in clipping operation [default 1]"<< endl;
    cout << "  -normalsweighting              : For the individual branch clustering, the weighting to put on normals. Should vary between 0 and 1. 1.0 will cluster only based on surface normals. 0.0 will cluster only based on position around the centerline [default 0.8]" << endl;
    cout << "  -isvasculature                 : Flag to indicate whether model is a vascular model with truncated boundaries. If model is not vasculature, the ends of the centerlines must be removed and the ends of the vessels need to be clustered based on position [default 1]" << endl;
    cout << "  -enforceboundarydirections     : At separating patches, enforce the boundary directions by modifying the clustering vectors [default 0]" << endl;
    cout << "  -boundaryenforcefactor         : Approximately represents the number of centerline points to enforce per branch. Typically a fairly low integer works well. The larger the value, the larger the portion of the vessel is set explicitly, and sometimes this can cause large problems. [default 1]" << endl;
    cout << "END COMMAND-LINE ARGUMENT SUMMARY" << endl;
    return EXIT_FAILURE;
  }
  if (!OutputProvided)
  {
    cout << "WARNING: Output Filename not provided, setting output name based on the input filename" <<endl;
    std::string newDirName = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename);
    // Only mac and linux!!!
    system(("mkdir -p "+newDirName).c_str());
    outputFilename = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"_Attributes_Passer.vtp";
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
  vtkNew(vtkSVSurfaceCenterlineAttributesPasser, Passer);

  Passer->SetInputData(inputPd);
  Passer->SetMergedCenterlines(centerlinesPd);
  if (PolycubeProvided)
  {
    Passer->SetPolycubePd(polycubePd);
  }
  Passer->SetCenterlineGroupIdsArrayName(groupIdsArrayName.c_str());
  Passer->SetGroupIdsArrayName(groupIdsArrayName.c_str());
  Passer->SetCenterlineRadiusArrayName(radiusArrayName.c_str());
  Passer->SetBlankingArrayName(blankingArrayName.c_str());
  Passer->SetCenterlineIdsArrayName("CenterlineIds");
  Passer->SetTractIdsArrayName("TractIds");
  Passer->SetPatchIdsArrayName(patchIdsArrayName.c_str());
  Passer->SetSlicePointsArrayName(slicePointsArrayName.c_str());
  Passer->SetParallelTransportVectorArrayName(parallelTransportVectorArrayName.c_str());
  Passer->SetTransformedNormalsArrayName(transformedNormalsArrayName.c_str());
  Passer->SetUseRadiusInformation(useRadiusInfo);
  Passer->SetIsVasculature(isVasculature);
  Passer->SetEnforceBoundaryDirections(enforceBoundaryDirections);
  Passer->SetNormalsWeighting(normalsWeighting);
  Passer->SetBoundaryEnforceFactor(boundaryEnforceFactor);
  Passer->DebugOn();
  Passer->Update();
  std::cout<<"Done"<<endl;

  //Write Files
  std::cout<<"Writing Files..."<<endl;
  vtkSVIOUtils::WriteVTPFile(outputFilename, Passer->GetOutput(0));

  if (Passer->GetErrorCode() != 0)
  {
    std::cerr << "Surface centerline grouper failed. " <<endl;
    return EXIT_FAILURE;
  }

  //Exit the program without errors
  return EXIT_SUCCESS;
}
