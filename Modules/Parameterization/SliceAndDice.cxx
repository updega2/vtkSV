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

/** @file SliceAndDice.cxx
 *  @brief This implements the vtkSVPolyDataSliceAndDiceFilter filter as a class
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkSVPolyDataSliceAndDiceFilter.h"

#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"

#include <string>
#include <sstream>
#include <iostream>

/**
 * \brief This creates an executable to process, segment, and create a
 * polycube of an arbitrary vascular model using its centerlines.
 */

int main(int argc, char *argv[])
{
  // BEGIN PROCESSING COMMAND-LINE ARGUMENTS
  // Assume the command line is okay
  bool BogusCmdLine = false;
  // Assume no options specified at command line
  bool RequestedHelp = false;
  bool InputProvided = false;
  bool CenterlinesProvided = false;
  bool OutputProvided = false;

  // Variables used in processing the commandline
  int iarg, arglength;
  std::string tmpstr;

  // Filenames
  std::string inputFileName;
  std::string centerlinesFileName;
  std::string outputFileName;

  // Default values for options
  double sliceLength = 1.0;
  int constructPolycube = 0;
  std::string boundaryPointsArrayName = "BoundaryPoints";
  std::string groupIdsArrayName = "GroupIds";
  std::string segmentIdsArrayName = "SegmentIds";
  std::string sliceIdsArrayName = "SliceIds";
  std::string radiusArrayName = "MaximumInscribedSphereRadius";
  std::string internalIdsArrayName = "TmpInternalIds";
  std::string dijkstraDistanceArrayName = "DijkstraDistance";

  // argc is the number of strings on the command-line
  //  starting with the program name
  for(iarg=1; iarg<argc; iarg++){
      arglength = strlen(argv[iarg]);
      // replace 0..arglength-1 with argv[iarg]
      tmpstr.replace(0,arglength,argv[iarg],0,arglength);
      if(tmpstr=="-h") {RequestedHelp = true;}
      else if(tmpstr=="-input") {InputProvided = true; inputFileName = argv[++iarg];}
      else if(tmpstr=="-centerlines") {CenterlinesProvided = true; centerlinesFileName = argv[++iarg];}
      else if(tmpstr=="-output") {OutputProvided = true; outputFileName = argv[++iarg];}
      else if(tmpstr=="-boundarypoints") {boundaryPointsArrayName = argv[++iarg];}
      else if(tmpstr=="-groupids") {groupIdsArrayName = argv[++iarg];}
      else if(tmpstr=="-segmentids") {segmentIdsArrayName = argv[++iarg];}
      else if(tmpstr=="-sliceids") {sliceIdsArrayName = argv[++iarg];}
      else if(tmpstr=="-radius") {radiusArrayName = argv[++iarg];}
      else if(tmpstr=="-internalids") {internalIdsArrayName = argv[++iarg];}
      else if(tmpstr=="-dijkstra") {dijkstraDistanceArrayName = argv[++iarg];}
      else if(tmpstr=="-polycube") {constructPolycube = atoi(argv[++iarg]);}
      else {BogusCmdLine = true;}
      // reset tmpstr for next argument
      tmpstr.erase(0,arglength);
  }

  if (RequestedHelp || !InputProvided || !CenterlinesProvided)
  {
    cout << endl;
    cout << "usage:" <<endl;
    cout << "  SliceAndDice -input [Input Filename] -centerlines [Centerlines Filename] -output [Output Filename] ..." << endl;
    cout << endl;
    cout << "COMMAND-LINE ARGUMENT SUMMARY" << endl;
    cout << "  -h                  : Display usage and command-line argument summary"<< endl;
    cout << "  -input              : Input file name (.vtp or .stl)"<< endl;
    cout << "  -centerlines        : Centerlines file name (.vtp)"<< endl;
    cout << "  -ouptut             : Output file name"<< endl;
    cout << "  -boundarypoints     : Name to be used for boundary points [default BoundaryPoints]"<< endl;
    cout << "  -groupids           : Name to be used for group ids [default GroupIds]"<< endl;
    cout << "  -segmentids         : Name to be used for segments ids [default SegmentIds]"<< endl;
    cout << "  -sliceids           : Name to be used for slice ids [default SliceIds]"<< endl;
    cout << "  -radius             : Name on centerlines describing maximum inscribed sphere radius [default MaximumInscribedSphereRadius]"<< endl;
    cout << "  -internalids        : Name to be used for ids used internal to filter [default TmpInternalIds]"<< endl;
    cout << "  -dijkstra           : Name to be used for distance calculated from dijkstra filter [default DijkstraDistance]"<< endl;
    cout << "  -polycube           : Construct polycube if turned on [default 0]"<< endl;
    cout << "END COMMAND-LINE ARGUMENT SUMMARY" << endl;
    return EXIT_FAILURE;
  }
  if (!OutputProvided)
  {
    cout << "WARNING: Output Filename not provided, setting output name based on the input filename" <<endl;
    std::string newDirName = vtkSVIOUtils::GetPath(inputFileName)+"/"+vtkSVIOUtils::GetRawName(inputFileName);
    // Only mac and linux!!!
    system(("mkdir -p "+newDirName).c_str());
    outputFileName = vtkSVIOUtils::GetPath(inputFileName)+"/"+vtkSVIOUtils::GetRawName(inputFileName)+"/"+vtkSVIOUtils::GetRawName(inputFileName)+"_Segmented.vtp";
  }

  // Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkNew(vtkPolyData, inputPd);
  vtkSVIOUtils::ReadInputFile(inputFileName,inputPd);
  vtkNew(vtkPolyData, centerlinesPd);
  vtkSVIOUtils::ReadInputFile(centerlinesFileName,centerlinesPd);

  // Filter
  vtkNew(vtkSVPolyDataSliceAndDiceFilter, Slicer);

  // OPERATION
  std::cout<<"Performing Operation..."<<endl;
  Slicer->SetInputData(inputPd);
  Slicer->SetCenterlinesPd(centerlinesPd);
  Slicer->SetSliceLength(sliceLength);
  Slicer->SetConstructPolycube(constructPolycube);
  Slicer->SetBoundaryPointsArrayName(boundaryPointsArrayName.c_str());
  Slicer->SetGroupIdsArrayName(groupIdsArrayName.c_str());
  Slicer->SetSegmentIdsArrayName(segmentIdsArrayName.c_str());
  Slicer->SetSliceIdsArrayName(sliceIdsArrayName.c_str());
  Slicer->SetSphereRadiusArrayName(radiusArrayName.c_str());
  Slicer->SetInternalIdsArrayName(internalIdsArrayName.c_str());
  Slicer->SetDijkstraArrayName(dijkstraDistanceArrayName.c_str());
  Slicer->Update();

  // Get output
  vtkNew(vtkPolyData, output);
  output = Slicer->GetOutput();

  // Write Files
  std::cout<<"Writing Files..."<<endl;
  vtkSVIOUtils::WriteVTPFile(outputFileName, output);
  std::cout<<"Done"<<endl;

  // Exit the program without errors
  return EXIT_SUCCESS;
}
