//
//  SliceAndDice.cxx
//
//
//  Created by Adam Updegrove on 10/4/14.
//
//

/*=========================================================================

 Program:   SimVascular

 =========================================================================*/

#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkDataArray.h"
#include "vtkDataWriter.h"
#include "vtkDoubleArray.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkIntArray.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkSVPolyDataSliceAndDiceFilter.h"

#include <string>
#include <sstream>
#include <iostream>

int main(int argc, char *argv[])
{
  /* BEGIN PROCESSING COMMAND-LINE ARGUMENTS */
  /* Assume the command line is okay */
  bool BogusCmdLine = false;
  /* Assume no options specified at command line */
  bool RequestedHelp = false;
  bool InputProvided = false;
  bool CenterlinesProvided = false;
  bool OutputProvided = false;

  /* variables used in processing the commandline */
  int iarg, arglength;
  std::string tmpstr;

  //Filenames
  std::string inputFileName;
  std::string centerlinesFileName;
  std::string outputFileName;

 //Default values for options
 double sliceLength = 1.0;
 int constructPolycube = 0;
 std::string boundaryPointsArrayName = "BoundaryPoints";
 std::string groupIdsArrayName = "GroupIds";
 std::string segmentIdsArrayName = "SegmentIds";
 std::string sliceIdsArrayName = "SliceIds";
 std::string radiusArrayName = "MaximumInscribedSphereRadius";
 std::string internalIdsArrayName = "TmpInternalIds";
 std::string dijkstraDistanceArrayName = "DijkstraDistance";

  /* argc is the number of strings on the command-line */
  /*  starting with the program name */
  for(iarg=1; iarg<argc; iarg++){
      arglength = strlen(argv[iarg]);
      /* replace 0..arglength-1 with argv[iarg] */
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
      /* reset tmpstr for next argument */
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
    //Only mac and linux!!!
    system(("mkdir -p "+newDirName).c_str());
    outputFileName = vtkSVIOUtils::GetPath(inputFileName)+"/"+vtkSVIOUtils::GetRawName(inputFileName)+"/"+vtkSVIOUtils::GetRawName(inputFileName)+"_Segmented.vtp";
  }

  //Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkNew(vtkPolyData, inputPd);
  vtkSVIOUtils::ReadInputFile(inputFileName,inputPd);
  vtkNew(vtkPolyData, centerlinesPd);
  vtkSVIOUtils::ReadInputFile(centerlinesFileName,centerlinesPd);

  //Filter
  vtkNew(vtkSVPolyDataSliceAndDiceFilter, Slicer);

  //OPERATION
  std::cout<<"Performing Operation..."<<endl;
  Slicer->SetInputData(inputPd);
  Slicer->SetCenterlines(centerlinesPd);
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

  vtkNew(vtkPolyData, output);
  output = Slicer->GetOutput();

  //Write Files
  std::cout<<"Writing Files..."<<endl;
  vtkSVIOUtils::WriteVTPFile(outputFileName, output);
  std::cout<<"Done"<<endl;

  //Exit the program without errors
  return EXIT_SUCCESS;
}
