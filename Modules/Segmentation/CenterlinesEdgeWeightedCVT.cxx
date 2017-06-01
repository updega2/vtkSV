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
 *  \file EdgeWeightedCVT.cxx
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */

#include "vtkSVCenterlinesEdgeWeightedCVT.h"

#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkPolyDataNormals.h"
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
  // Assume no options specified at command line
  bool RequestedHelp = false;
  bool InputProvided = false;
  bool CVTDataArrayNameProvided = false;
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
  int numberOfRings             = 2;
  int maximumNumberOfIterations = 1.0e3;
  double threshold              = 2;
  double edgeWeight             = 1.0;
  int useRadiusInfo             = 1;

  std::string patchIdsArrayName = "GroupIds";
  std::string cvtDataArrayName  = "CellNormals";
  std::string groupIdsArrayName = "GroupIds";
  std::string radiusArrayName   = "MaximumInscribedSphereRadius";
  std::string blankingArrayName = "Blanking";

  // argc is the number of strings on the command-line
  //  starting with the program name
  for(iarg=1; iarg<argc; iarg++){
      arglength = strlen(argv[iarg]);
      // replace 0..arglength-1 with argv[iarg]
      tmpstr.replace(0,arglength,argv[iarg],0,arglength);
      if(tmpstr=="-h")                       {RequestedHelp = true;}
      else if(tmpstr=="-input")              {InputProvided = true; inputFilename = argv[++iarg];}
      else if(tmpstr=="-centerlines")        {CenterlinesProvided = true; centerlinesFilename = argv[++iarg];}
      else if(tmpstr=="-output")             {OutputProvided = true; outputFilename = argv[++iarg];}
      else if(tmpstr=="-patchids")           {patchIdsArrayName = argv[++iarg];}
      else if(tmpstr=="-cvtdata")            {CVTDataArrayNameProvided; cvtDataArrayName = argv[++iarg];}
      else if(tmpstr=="-numberofrings")      {numberOfRings = atoi(argv[++iarg]);}
      else if(tmpstr=="-maximumnumberofiterations")        {maximumNumberOfIterations = atoi(argv[++iarg]);}
      else if(tmpstr=="-threshold")          {threshold = atof(argv[++iarg]);}
      else if(tmpstr=="-edgeweight")         {edgeWeight = atof(argv[++iarg]);}
      else if(tmpstr=="-groupids")           {groupIdsArrayName = argv[++iarg];}
      else if(tmpstr=="-radius")             {radiusArrayName = argv[++iarg];}
      else if(tmpstr=="-blanking")           {blankingArrayName = argv[++iarg];}
      else if(tmpstr=="-useradiusinfo")      {useRadiusInfo = atoi(argv[++iarg]);}
      else {cout << argv[iarg] << " is not a  valid argument. Ask for help with -h." << endl; RequestedHelp = true; return EXIT_FAILURE;}
      // reset tmpstr for next argument
      tmpstr.erase(0,arglength);
  }

  if (RequestedHelp || !InputProvided || !CenterlinesProvided)
  {
    cout << endl;
    cout << "usage:" <<endl;
    cout << "  EdgeWeightedCVT -input [Input Filename] -centerlines [Centerlines Filename] -output [Output Filename] ..." << endl;
    cout << endl;
    cout << "COMMAND-LINE ARGUMENT SUMMARY" << endl;
    cout << "  -h                  : Display usage and command-line argument summary"<< endl;
    cout << "  -input              : Input file name (.vtp or .stl)"<< endl;
    cout << "  -centerlines        : Centerlines file name (.vtp)"<< endl;
    cout << "  -output             : Output file name"<< endl;
    cout << "  -numberofrings      : Number of rings to consider neighborhood of an element [default 2]"<< endl;
    cout << "  -patchids           : Name to be used for patches found using cvt [default GroupIds]"<< endl;
    cout << "  -cvtdata            : Name on input with data to be used for cvt [no default]"<< endl;
    cout << "  -threshold          : Threshold criteria for when to stop cvt iterations [default 2]"<< endl;
    cout << "  -edgeweight         : How much influence the edge weighting term should have [default 1.0]"<< endl;
    cout << "  -maximumnumberofiterations        : Set a maximum number of iterations [default 1000]"<< endl;
    cout << "  -groupids           : Name to be used for group ids [default GroupIds]"<< endl;
    cout << "  -radius             : Name on centerlines describing maximum inscribed sphere radius [default MaximumInscribedSphereRadius]"<< endl;
    cout << "  -blanking           : Name on centerlines describing whether line is part of bifurcation region or not [default Blanking]"<< endl;
    cout << "  -useradiusinfo      : Use radius to help in clipping operation [default 1]"<< endl;
    cout << "END COMMAND-LINE ARGUMENT SUMMARY" << endl;
    return EXIT_FAILURE;
  }
  if (!OutputProvided)
  {
    cout << "WARNING: Output Filename not provided, setting output name based on the input filename" <<endl;
    std::string newDirName = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename);
    // Only mac and linux!!!
    system(("mkdir -p "+newDirName).c_str());
    outputFilename = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"_CVT.vtp";
  }

  std::cout<<"Reading Files..."<<endl;
  vtkNew(vtkPolyData, inputPd);
  if (vtkSVIOUtils::ReadInputFile(inputFilename,inputPd) != 1)
    return EXIT_FAILURE;
  vtkNew(vtkPolyData, centerlinesPd);
  if (vtkSVIOUtils::ReadInputFile(centerlinesFilename,centerlinesPd) != 1)
    return EXIT_FAILURE;

  // Using normals as cvtdata
  if (!CVTDataArrayNameProvided)
  {
    vtkNew(vtkPolyDataNormals, normaler);
    normaler->SetInputData(inputPd);
    normaler->SplittingOff();
    normaler->AutoOrientNormalsOn();
    normaler->ComputePointNormalsOff();
    normaler->ComputeCellNormalsOn();
    normaler->Update();

    vtkFloatArray *floatArray = vtkFloatArray::SafeDownCast(normaler->GetOutput()->GetCellData()->GetArray("Normals"));
    vtkNew(vtkDoubleArray, doubleArray);
    doubleArray->SetNumberOfComponents(3);
    doubleArray->SetNumberOfTuples(floatArray->GetNumberOfTuples());
    for (int i=0; i<floatArray->GetNumberOfTuples(); i++)
    {
      for (int j=0; j<3; j++)
      {
        doubleArray->SetComponent(i, j, floatArray->GetComponent(i, j));
      }
    }

    doubleArray->SetName("CellNormals");
    inputPd->GetCellData()->AddArray(doubleArray);
  }

  // Call Function to Read File
  // Filter
  vtkNew(vtkSVCenterlinesEdgeWeightedCVT, CVT);

  // OPERATION
  std::cout<<"Performing Operation..."<<endl;
  CVT->SetInputData(inputPd);
  CVT->SetGenerators(centerlinesPd);
  CVT->SetNumberOfRings(numberOfRings);
  CVT->SetThreshold(threshold);
  CVT->SetEdgeWeight(edgeWeight);
  CVT->SetMaximumNumberOfIterations(maximumNumberOfIterations);
  CVT->SetPatchIdsArrayName(patchIdsArrayName.c_str());
  CVT->SetCVTDataArrayName(cvtDataArrayName.c_str());
  CVT->SetGroupIdsArrayName(groupIdsArrayName.c_str());
  CVT->SetCenterlineRadiusArrayName(radiusArrayName.c_str());
  CVT->SetBlankingArrayName(blankingArrayName.c_str());
  CVT->SetUseRadiusInformation(useRadiusInfo);
  CVT->Update();

  // Get output
  vtkNew(vtkPolyData, output);
  output = CVT->GetOutput();

  // Write Files
  std::cout<<"Writing Files..."<<endl;
  vtkSVIOUtils::WriteVTPFile(outputFilename, output);
  std::cout<<"Done"<<endl;

  // Exit the program without errors
  return EXIT_SUCCESS;
}
