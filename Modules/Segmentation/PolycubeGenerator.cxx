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
#include "vtkDataSetSurfaceFilter.h"
#include "vtkFeatureEdges.h"
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
#include "vtkSVCleanUnstructuredGrid.h"
#include "vtkSVSeedSelector.h"
#include "vtkSVPickPointSeedSelector.h"
#include "vtkSVOpenProfilesSeedSelector.h"

#include "vtkSVPolycubeGenerator.h"
#include "vtkvmtkCapPolyData.h"
#include "vtkvmtkPolyDataCenterlines.h"

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
  std::string groupIdsArrayName = "GroupIds";
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
      else if(tmpstr=="-radiusname")       {radiusArrayName = argv[++iarg];}
      else if(tmpstr=="-groupidsname")     {radiusArrayName = argv[++iarg];}
      else {cout << argv[iarg] << " is not a valid argument. Ask for help with -h." << endl; RequestedHelp = true; return EXIT_FAILURE;}
      // reset tmpstr for next argument
      tmpstr.erase(0,arglength);
  }

  if (RequestedHelp || !InputProvided)
  {
    cout << endl;
    cout << "usage:" <<endl;
    cout << "  PolycubeGenerator -input [Input Filename] -output [Output Filename] ..." << endl;
    cout << endl;
    cout << "COMMAND-LINE ARGUMENT SUMMARY" << endl;
    cout << "  -h                  : Display usage and command-line argument summary"<< endl;
    cout << "  -input              : Input file name (.vtp)"<< endl;
    cout << "  -output             : Output file name"<< endl;
    cout << "  -radiusname         : Name on centerlines describing maximum inscribed sphere radius [default MaximumInscribedSphereRadius]"<< endl;
    cout << "  -groupidsname       : Name on centerlines differentiating the different centerlines [default GroupIds]"<< endl;
    cout << "END COMMAND-LINE ARGUMENT SUMMARY" << endl;
    return EXIT_FAILURE;
  }
  if (!OutputProvided)
  {
    cout << "WARNING: Output Filename not provided, setting output name based on the input filename" <<endl;
    std::string newDirName = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename);
    // Only mac and linux!!!
    system(("mkdir -p "+newDirName).c_str());
    outputFilename = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"_Polycube.vtp";
  }

  // Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkNew(vtkPolyData, inputPd);
  if (vtkSVIOUtils::ReadInputFile(inputFilename,inputPd) != 1)
    return EXIT_FAILURE;

  std::cout<<"Merging Centerlines..."<<endl;
  vtkNew(vtkSVPolycubeGenerator, Polycuber);
  Polycuber->SetInputData(inputPd);
  Polycuber->SetCenterlineRadiusArrayName(radiusArrayName.c_str());
  Polycuber->SetCenterlineGroupIdsArrayName(groupIdsArrayName.c_str());
  Polycuber->Update();
  std::cout<<"Done"<<endl;

  //Write Files
  std::cout<<"Writing Files..."<<endl;
  vtkSVIOUtils::WriteVTPFile(outputFilename, Polycuber->GetOutput(0));
  vtkSVIOUtils::WriteVTUFile(outputFilename, Polycuber->GetVolumePolycubeUg(), "_Volume");
  vtkSVIOUtils::WriteVTPFile(outputFilename, Polycuber->GetGraphPd(), "_Graph");
  vtkSVIOUtils::WriteVTPFile(outputFilename, Polycuber->GetCenterlineGraph()->Lines, "_Centerlines");

  // Check to make sure correct output
  int numEnds = 0;
  vtkIdType npts, *pts;
  vtkNew(vtkIdList, point0Cells);
  vtkNew(vtkIdList, pointNCells);
  for (int i=0; i<inputPd->GetNumberOfCells(); i++)
  {
    inputPd->GetCellPoints(i, npts, pts);

    inputPd->GetPointCells(pts[0], point0Cells);
    inputPd->GetPointCells(pts[npts-1], pointNCells);

    if (point0Cells->GetNumberOfIds() == 1)
    {
      numEnds++;
    }
    if (pointNCells->GetNumberOfIds() == 1)
    {
      numEnds++;
    }
  }

  int numPolycubeCells;
  if (numEnds == 2)
  {
    numPolycubeCells = 6;
  }
  else
  {
    numPolycubeCells = numEnds*5 + (inputPd->GetNumberOfCells() - numEnds) * 4;
  }

  vtkNew(vtkFeatureEdges, featurer);
  featurer->SetInputData(Polycuber->GetOutput());
  featurer->BoundaryEdgesOn();
  featurer->FeatureEdgesOff();
  featurer->ManifoldEdgesOff();
  featurer->NonManifoldEdgesOn();
  featurer->Update();

  if (featurer->GetOutput()->GetNumberOfPoints() != 0 || featurer->GetOutput()->GetNumberOfPoints() != 0)
  {
    std::cerr << "Open edges on surface polycube. "<< featurer->GetOutput()->GetNumberOfPoints() <<" points and " << featurer->GetOutput()->GetNumberOfCells() << " cells." <<endl;
    return EXIT_FAILURE;
  }

  if (Polycuber->GetGraphPd()->GetNumberOfCells() != inputPd->GetNumberOfCells())
  {
    std::cerr << "Incorrect number of cells on graph: "<< Polycuber->GetGraphPd()->GetNumberOfCells() <<". But should be " << inputPd->GetNumberOfCells() <<endl;
    return EXIT_FAILURE;
  }

  if (Polycuber->GetOutput()->GetNumberOfCells() != numPolycubeCells)
  {
    std::cerr << "Incorrect number of cells on surface polycube: "<< Polycuber->GetOutput()->GetNumberOfCells() <<". But should be " << numPolycubeCells <<endl;
    return EXIT_FAILURE;
  }

  vtkDataArray *polycubeDivisions = Polycuber->GetVolumePolycubeUg()->GetFieldData()->GetArray("PolycubeDivisions");
  if (polycubeDivisions == NULL)
  {
    std::cerr << "Field data array with name PolycubeDivisions does not exist on volume" << endl;
    return EXIT_FAILURE;
  }
  if (polycubeDivisions->GetNumberOfTuples() != inputPd->GetNumberOfCells())
  {
    std::cerr << "Incorrect number of cells on volume polycube: "<< polycubeDivisions->GetNumberOfTuples() <<". But should be " << inputPd->GetNumberOfCells() <<endl;
    return EXIT_FAILURE;
  }

  double divs[4];
  int totalPoints = 0;
  for (int i=0; i<polycubeDivisions->GetNumberOfTuples(); i++)
  {
    polycubeDivisions->GetTuple(i, divs);
    totalPoints += divs[1] * divs[2] * divs[3];
  }

  if (totalPoints != Polycuber->GetVolumePolycubeUg()->GetNumberOfPoints())
  {
    std::cerr << "Incorrect number of points on volume polycube: "<< totalPoints <<". But should be " << Polycuber->GetVolumePolycubeUg()->GetNumberOfPoints() <<endl;
    return EXIT_FAILURE;
  }

  vtkNew(vtkSVCleanUnstructuredGrid, cleaner);
  cleaner->SetInputData(Polycuber->GetVolumePolycubeUg());
  cleaner->ToleranceIsAbsoluteOn();
  cleaner->SetAbsoluteTolerance(1.0e-6);
  cleaner->Update();

  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(cleaner->GetOutput());
  surfacer->Update();

  featurer->SetInputData(surfacer->GetOutput());
  featurer->BoundaryEdgesOn();
  featurer->FeatureEdgesOff();
  featurer->ManifoldEdgesOff();
  featurer->NonManifoldEdgesOff();
  featurer->Update();

  if (featurer->GetOutput()->GetNumberOfPoints() != 0 || featurer->GetOutput()->GetNumberOfPoints() != 0)
  {
    std::cerr << "Open edges on volume polycube. "<< featurer->GetOutput()->GetNumberOfPoints() <<" points and " << featurer->GetOutput()->GetNumberOfCells() << " cells." <<endl;
    return EXIT_FAILURE;
  }


  //Exit the program without errors
  return EXIT_SUCCESS;
}
