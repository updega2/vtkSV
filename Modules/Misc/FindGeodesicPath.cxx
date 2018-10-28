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
  bool StartPtIdProvided = false;
  bool EndPtIdProvided   = false;
  bool ClosePtProvided   = false;
  bool OutputProvided    = false;

  // Variables used in processing the commandline
  int iarg, arglength;
  std::string tmpstr;

  // Filenames
  std::string inputFilename;
  std::string outputFilename;

  // Default values for options
  int startPtId        = -1;
  int endPtId          = -1;
  int repelBoundaryPts = 0;
  double closePt[3]; closePt[0] = 0.0; closePt[1] = 0.0; closePt[2] = 0.0;
  std::string dijkstraArrayName    = "DijkstraDistance";
  std::string internalIdsArrayName = "TmpInternalIds";
  std::string pathBooleanArrayName = "IsPath";

  // argc is the number of strings on the command-line
  //  starting with the program name
  for(iarg=1; iarg<argc; iarg++){
      arglength = strlen(argv[iarg]);
      // replace 0..arglength-1 with argv[iarg]
      tmpstr.replace(0,arglength,argv[iarg],0,arglength);
      if(tmpstr=="-h")                      {RequestedHelp = true;}
      else if(tmpstr=="-input")             {InputProvided = true; inputFilename = argv[++iarg];}
      else if(tmpstr=="-output")            {OutputProvided = true; outputFilename = argv[++iarg];}
      else if(tmpstr=="-startptid")         {StartPtIdProvided = true; startPtId = atoi(argv[++iarg]);}
      else if(tmpstr=="-endptid")           {EndPtIdProvided = true; endPtId = atoi(argv[++iarg]);}
      else if(tmpstr=="-closept")           {ClosePtProvided = true; closePt[0] = atof(argv[++iarg]); closePt[1] = atof(argv[++iarg]); closePt[2] = atof(argv[++iarg]);}
      else if(tmpstr=="-internalids")       {internalIdsArrayName = argv[++iarg];}
      else if(tmpstr=="-dijkstra")          {dijkstraArrayName = argv[++iarg];}
      else if(tmpstr=="-pathboolean")       {pathBooleanArrayName = argv[++iarg];}
      else if(tmpstr=="-repelboundarypts")  {repelBoundaryPts = atoi(argv[++iarg]);}
      else {cout << argv[iarg] << " is not a valid argument. Ask for help with -h." << endl; RequestedHelp = true; return EXIT_FAILURE;}
      // reset tmpstr for next argument
      tmpstr.erase(0,arglength);
  }

  if (RequestedHelp || !InputProvided || !StartPtIdProvided)
  {
    cout << endl;
    cout << "usage:" <<endl;
    cout << "  FindGeodesicPath -input [Input Filename] -output [Output Filename] -startptid [Start Pt Id] ..." << endl;
    cout << endl;
    cout << "COMMAND-LINE ARGUMENT SUMMARY" << endl;
    cout << "  -h                  : Display usage and command-line argument summary"<< endl;
    cout << "  -input              : Input file name (.vtp or .stl)"<< endl;
    cout << "  -output             : Output file name"<< endl;
    cout << "  -startptid          : Point Id to start from [default 0]"<< endl;
    cout << "  -endptid            : Point Id to end at [default not used]"<< endl;
    cout << "  -closept            : If end pt id, is not provided, this 3d point location can be used to get the boundary closest to the point to use as the ending location of the shortes distance calculation [default not used, usage [-closept 0.0 0.0 0.0]"<< endl;
    cout << "  -internalids        : Name to be used for ids used internal to filter [default TmpInternalIds]"<< endl;
    cout << "  -dijkstra           : Name to be used for distance calculated from dijkstra filter [default DijkstraDistance]"<< endl;
    cout << "  -pathboolean        : Name to be used for point data indicating wheterh point is on shortes distance path [default IsPath]"<< endl;
    cout << "  -repelboundarypts   : Will ask the shortest distance algorithm to stay away from all other points on a boundary near the start and end point ids [default 0]"<< endl;
    cout << "END COMMAND-LINE ARGUMENT SUMMARY" << endl;
    return EXIT_FAILURE;
  }
  if (!EndPtIdProvided && !ClosePtProvided)
  {
    cout << "WARNING: No end point id, '-endpointid',  and no point close to the ending boundary, '-closept', were provided, which means no path will be found. Output will contain array indicating distance of every point from the input point but that is it." <<endl;
  }
  if (EndPtIdProvided && ClosePtProvided)
  {
    cout << "WARNING: An end point id, '-endpointid',  and a point close to the ending boundary, '-closept', were provided, which is impossible to do. Ending point id will be used" <<endl;
  }
  if (!OutputProvided)
  {
    cout << "WARNING: Output Filename not provided, setting output name based on the input filename" <<endl;
    std::string newDirName = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename);
    // Only mac and linux!!!
    system(("mkdir -p "+newDirName).c_str());
    outputFilename = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"_Geodescized.vtp";
  }

  // Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkNew(vtkPolyData, inputPd);
  if (vtkSVIOUtils::ReadInputFile(inputFilename,inputPd) != 1)
    return EXIT_FAILURE;

  // Filter
  vtkNew(vtkSVFindGeodesicPath, Finder);

  //OPERATION
  std::cout<<"Performing Operation..."<<endl;
  Finder->SetInputData(0, inputPd);
  Finder->SetStartPtId(startPtId);
  Finder->SetEndPtId(endPtId);
  Finder->SetClosePt(closePt);
  if (ClosePtProvided || EndPtIdProvided)
    Finder->SetAddPathBooleanArray(1);
  Finder->SetDijkstraArrayName(dijkstraArrayName.c_str());
  Finder->SetInternalIdsArrayName(internalIdsArrayName.c_str());
  Finder->SetPathBooleanArrayName(pathBooleanArrayName.c_str());
  if (repelBoundaryPts)
    Finder->SetRepelCloseBoundaryPoints(1);
  Finder->Update();
  std::cout<<"Done"<<endl;

  if (Finder->GetErrorCode() != 0)
  {
    std::cerr << "Error in filter" << endl;
    return EXIT_FAILURE;
  }

  //Write Files
  std::cout<<"Writing Files..."<<endl;
  vtkSVIOUtils::WriteVTPFile(outputFilename, Finder->GetOutput(0));

  //Exit the program without errors
  return EXIT_SUCCESS;
}
