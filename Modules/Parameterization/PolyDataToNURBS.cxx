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

/**
 *  \file PolyDataToNURBS.cxx
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */

#include "vtkSVPolyDataToNURBSFilter.h"

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"

#include <string>
#include <sstream>
#include <iostream>

/**
 * \brief This creates an executable to take an input vasuclar surface and
 * its centerlines and form a globally valid parameterization.
 */
int main(int argc, char *argv[])
{
  // BEGIN PROCESSING COMMAND-LINE ARGUMENTS
  // Assume no options specified at command line
  bool RequestedHelp = false;
  bool InputProvided = false;
  bool CenterlinesProvided = false;
  bool OutputProvided = false;

  // Variables used in processing the commandline
  int iarg, arglength;
  std::string tmpstr;

  // Filenames
  std::string inputFilename;
  std::string centerlinesFilename;
  std::string outputFilename;

  // Default values for options
  int BaseDomainXResolution = 36;
  int BaseDomainYResolution = 9;
  int writeTexturedPd       = 0;
  int writeLoftedPd         = 0;
  std::string boundaryPointsArrayName   = "BoundaryPoints";
  std::string groupIdsArrayName         = "GroupIds";
  std::string segmentIdsArrayName       = "SegmentIds";
  std::string sliceIdsArrayName         = "SliceIds";
  std::string radiusArrayName           = "MaximumInscribedSphereRadius";
  std::string internalIdsArrayName      = "TmpInternalIds";
  std::string dijkstraDistanceArrayName = "DijkstraDistance";
  std::string pathBooleanArrayName      = "IsPath";

  // argc is the number of strings on the command-line
  //  starting with the program name
  for(iarg=1; iarg<argc; iarg++){
      arglength = strlen(argv[iarg]);
      // replace 0..arglength-1 with argv[iarg]
      tmpstr.replace(0,arglength,argv[iarg],0,arglength);
      if(tmpstr=="-h")                      {RequestedHelp = true;}
      else if(tmpstr=="-input")             {InputProvided = true; inputFilename = argv[++iarg];}
      else if(tmpstr=="-centerlines")       {CenterlinesProvided = true; centerlinesFilename = argv[++iarg];}
      else if(tmpstr=="-output")            {OutputProvided = true; outputFilename = argv[++iarg];}
      else if(tmpstr=="-boundarypoints")    {boundaryPointsArrayName = argv[++iarg];}
      else if(tmpstr=="-groupids")          {groupIdsArrayName = argv[++iarg];}
      else if(tmpstr=="-segmentids")        {segmentIdsArrayName = argv[++iarg];}
      else if(tmpstr=="-sliceids")          {sliceIdsArrayName = argv[++iarg];}
      else if(tmpstr=="-radius")            {radiusArrayName = argv[++iarg];}
      else if(tmpstr=="-internalids")       {internalIdsArrayName = argv[++iarg];}
      else if(tmpstr=="-dijkstra")          {dijkstraDistanceArrayName = argv[++iarg];}
      else if(tmpstr=="-pathboolean")       {pathBooleanArrayName = argv[++iarg];}
      else if(tmpstr=="-writetextured")     {writeTexturedPd = atoi(argv[++iarg]);}
      else if(tmpstr=="-writelofted")       {writeLoftedPd = atoi(argv[++iarg]);}
      else {cout << argv[iarg] << " is not a valid argument. Ask for help with -h." << endl; RequestedHelp = true; return EXIT_FAILURE;}
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
    cout << "  -output             : Output file name"<< endl;
    cout << "  -boundarypoints     : Name to be used for boundary points [default BoundaryPoints]"<< endl;
    cout << "  -groupids           : Name to be used for group ids [default GroupIds]"<< endl;
    cout << "  -segmentids         : Name to be used for segments ids [default SegmentIds]"<< endl;
    cout << "  -sliceids           : Name to be used for slice ids [default SliceIds]"<< endl;
    cout << "  -radius             : Name on centerlines describing maximum inscribed sphere radius [default MaximumInscribedSphereRadius]"<< endl;
    cout << "  -internalids        : Name to be used for ids used internal to filter [default TmpInternalIds]"<< endl;
    cout << "  -dijkstra           : Name to be used for distance calculated from dijkstra filter [default DijkstraDistance]"<< endl;
    cout << "  -polycube           : Construct polycube if turned on [default 1]"<< endl;
    cout << "  -writetextured      : Write the input with the texture map coordinates [default 0]"<< endl;
    cout << "  -writelofted        : Write the lofted NURBS polydata representation [default 0]"<< endl;
    cout << "END COMMAND-LINE ARGUMENT SUMMARY" << endl;
    return EXIT_FAILURE;
  }
  if (!OutputProvided)
  {
    cout << "WARNING: Output Filename not provided, setting output name based on the input filename" <<endl;
    std::string newDirName = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename);
    // Only mac and linux!!!
    system(("mkdir -p "+newDirName).c_str());
    outputFilename = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"_Converted.vtp";
  }

  // Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkNew(vtkPolyData, inputPd);
  if (vtkSVIOUtils::ReadInputFile(inputFilename,inputPd) != 1)
    return EXIT_FAILURE;
  vtkNew(vtkPolyData, centerlinesPd);
  if (vtkSVIOUtils::ReadInputFile(centerlinesFilename,centerlinesPd) != 1)
    return EXIT_FAILURE;

  // Filter
  vtkNew(vtkSVPolyDataToNURBSFilter, Converter);

  //OPERATION
  std::cout<<"Performing Operation..."<<endl;
  Converter->SetInputData(inputPd);
  Converter->SetCenterlinesPd(centerlinesPd);
  Converter->SetAddTextureCoordinates(1);
  Converter->SetBoundaryPointsArrayName(boundaryPointsArrayName.c_str());
  Converter->SetGroupIdsArrayName(groupIdsArrayName.c_str());
  Converter->SetSegmentIdsArrayName(segmentIdsArrayName.c_str());
  Converter->SetSliceIdsArrayName(sliceIdsArrayName.c_str());
  Converter->SetSphereRadiusArrayName(radiusArrayName.c_str());
  Converter->SetInternalIdsArrayName(internalIdsArrayName.c_str());
  Converter->SetDijkstraArrayName(dijkstraDistanceArrayName.c_str());
  Converter->SetBooleanPathArrayName(pathBooleanArrayName.c_str());
  Converter->Update();

  //Write Files
  std::cout<<"Writing Files..."<<endl;
  vtkSVIOUtils::WriteVTPFile(outputFilename, Converter->GetOutput());
  if (writeTexturedPd)
    vtkSVIOUtils::WriteVTPFile(outputFilename, Converter->GetTexturedPd(), "_Textured");
  if (writeLoftedPd)
    vtkSVIOUtils::WriteVTPFile(outputFilename, Converter->GetLoftedPd(), "_Lofted");
  std::cout<<"Done"<<endl;

  //Exit the program without errors
  return EXIT_SUCCESS;
}
