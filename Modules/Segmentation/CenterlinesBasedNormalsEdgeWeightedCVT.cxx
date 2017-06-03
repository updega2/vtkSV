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

#include "vtkSVCenterlinesBasedNormals.h"
#include "vtkSVEdgeWeightedCVT.h"

#include "vtkAppendPolyData.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmartPointer.h"
#include "vtkSVGeneralUtils.h"
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
  bool GeneratorsProvided = false;
  bool CenterlinesProvided = false;
  bool OutputProvided = false;

  // Variables used in processing the commandline
  int iarg, arglength;
  std::string tmpstr;

  // Filenames
  std::string inputFilename;
  std::string generatorsFilename;
  std::string centerlinesFilename;
  std::string outputFilename;

  // Default values for options
  int usePointArray                     = 0;
  int useCellArray                      = 1;
  int useGeneratorsArray                 = 0;
  int numberOfRings                     = 2;
  int useTransferredPatchesAsThreshold  = 1;
  int maximumNumberOfIterations         = 1.0e3;
  double threshold                      = 2;
  double edgeWeight                     = 1.0;
  std::string patchIdsArrayName         = "PatchIds";
  std::string groupIdsArrayName         = "GroupIds";
  std::string cvtDataArrayName          = "CellNormals";
  std::string newCellArrayName          = "CenterlinesBasedCellNormals";
  std::string generatorsArrayName       = "SliceIds";

  // argc is the number of strings on the command-line
  //  starting with the program name
  for(iarg=1; iarg<argc; iarg++){
      arglength = strlen(argv[iarg]);
      // replace 0..arglength-1 with argv[iarg]
      tmpstr.replace(0,arglength,argv[iarg],0,arglength);
      if(tmpstr=="-h")                       {RequestedHelp = true;}
      else if(tmpstr=="-input")              {InputProvided = true; inputFilename = argv[++iarg];}
      else if(tmpstr=="-generators")         {GeneratorsProvided = true; generatorsFilename = argv[++iarg];}
      else if(tmpstr=="-centerlines")        {CenterlinesProvided = true; centerlinesFilename = argv[++iarg];}
      else if(tmpstr=="-output")             {OutputProvided = true; outputFilename = argv[++iarg];}
      else if(tmpstr=="-patchids")           {patchIdsArrayName = argv[++iarg];}
      else if(tmpstr=="-groupids")           {groupIdsArrayName = argv[++iarg];}
      else if(tmpstr=="-newcellarray")       {newCellArrayName = argv[++iarg];}
      else if(tmpstr=="-cvtdata")            {CVTDataArrayNameProvided; cvtDataArrayName = argv[++iarg];}
      else if(tmpstr=="-generatorsarray")    {generatorsArrayName = argv[++iarg];}
      else if(tmpstr=="-numberofrings")      {numberOfRings = atoi(argv[++iarg]);}
      else if(tmpstr=="-usetransferredpatchesasthreshold") {useTransferredPatchesAsThreshold = atoi(argv[++iarg]);}
      else if(tmpstr=="-maximumnumberofiterations")        {maximumNumberOfIterations = atoi(argv[++iarg]);}
      else if(tmpstr=="-threshold")          {threshold = atof(argv[++iarg]);}
      else if(tmpstr=="-edgeweight")         {edgeWeight = atof(argv[++iarg]);}
      else if(tmpstr=="-usepointarray")      {usePointArray = atoi(argv[++iarg]);}
      else if(tmpstr=="-usecellarray")       {useCellArray = atoi(argv[++iarg]);}
      else if(tmpstr=="-usegeneratorsarray") {useGeneratorsArray = atoi(argv[++iarg]);}
      else {cout << argv[iarg] << " is not a  valid argument. Ask for help with -h." << endl; RequestedHelp = true; return EXIT_FAILURE;}
      // reset tmpstr for next argument
      tmpstr.erase(0,arglength);
  }

  if (RequestedHelp || !InputProvided || !CenterlinesProvided)
  {
    cout << endl;
    cout << "usage:" <<endl;
    cout << "  EdgeWeightedCVT -input [Input Filename] -centerlines [Centerlines Filename]-generators [Generators Filename] -output [Output Filename] ..." << endl;
    cout << endl;
    cout << "COMMAND-LINE ARGUMENT SUMMARY" << endl;
    cout << "  -h                  : Display usage and command-line argument summary"<< endl;
    cout << "  -input              : Input file name (.vtp or .stl)"<< endl;
    cout << "  -generators         : Generators file name (.vtp)"<< endl;
    cout << "  -centerlines        : Centerlines file name (.vtp)"<< endl;
    cout << "  -output             : Output file name"<< endl;
    cout << "  -usepointarray      : Use data on points for cvt [default 0]"<< endl;
    cout << "  -usecellarray       : Use data on cellss for cvt [default 1]"<< endl;
    cout << "  -usegeneratorsarray : Instead of using generators points, use array on generator points as generator location. Usefull if more than three generator components [default 0]"<< endl;
    cout << "  -numberofrings      : Number of rings to consider neighborhood of an element [default 2]"<< endl;
    cout << "  -patchids           : Name to be used for patches found using cvt [default PatchIds]"<< endl;
    cout << "  -groupids           : Name to be used for group [default GroupIds]"<< endl;
    cout << "  -newcellarray       : Name to be given to output of centerlines normaler [default CenterlinesBasedCellNormals]"<< endl;
    cout << "  -cvtdata            : Name on input with data to be used for cvt [no default]"<< endl;
    cout << "  -generatorsarray    : Name on generators with data. Requied is usegeneratorsarray turned on [no default]"<< endl;
    cout << "  -threshold          : Threshold criteria for when to stop cvt iterations [default 2]"<< endl;
    cout << "  -edgeweight         : How much influence the edge weighting term should have [default 1.0]"<< endl;
    cout << "  -usetransferredpatchesasthreshold : Instead of using an energy criteria, the number of cells changing corresponding generators in an iteration determines convergence [default 1]"<< endl;
    cout << "  -maximumnumberofiterations        : Set a maximum number of iterations [default 1000]"<< endl;
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

  vtkNew(vtkPolyData, generatorsPd);
  if (!GeneratorsProvided)
  {
    vtkNew(vtkPoints, generatorsPts);
    generatorsPts->SetNumberOfPoints(6);
    generatorsPts->SetPoint(0, 1.0, 0.0, 0.0);
    generatorsPts->SetPoint(1, -1.0, 0.0, 0.0);
    generatorsPts->SetPoint(2, 0.0, 1.0, 0.0);
    generatorsPts->SetPoint(3, 0.0, -1.0, 0.0);
    generatorsPts->SetPoint(4, 0.0, 0.0, 1.0);
    generatorsPts->SetPoint(5, 0.0, 0.0, -1.0);

    generatorsPd->SetPoints(generatorsPts);
  }
  else
  {
    if (vtkSVIOUtils::ReadInputFile(generatorsFilename,generatorsPd) != 1)
      return EXIT_FAILURE;
  }

  // Using normals as cvtdata
  if (!CVTDataArrayNameProvided)
  {
    useCellArray  = 1;
    usePointArray = 0;

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

  vtkNew(vtkSVCenterlinesBasedNormals, newNormaler);
  newNormaler->SetInputData(inputPd);
  newNormaler->SetCenterlinesPd(centerlinesPd);
  newNormaler->SetCellArrayName(cvtDataArrayName.c_str());
  newNormaler->SetGroupIdsArrayName(groupIdsArrayName.c_str());
  newNormaler->SetNewCellArrayName(newCellArrayName.c_str());
  newNormaler->SetUsePointArray(usePointArray);
  newNormaler->SetUseCellArray(useCellArray);
  newNormaler->Update();

  //vtkNew(vtkIdList, centerlineGroupIds);
  //for (int i=0; i<inputPd->GetCellData()->GetArray(groupIdsArrayName.c_str())->GetNumberOfTuples(); i++)
  //{
  //  centerlineGroupIds->InsertUniqueId(static_cast<vtkIdType>(vtkMath::Round(inputPd->GetCellData()->GetArray(groupIdsArrayName.c_str())->GetComponent(i,0))));
  //}
  //int numGroups = centerlineGroupIds->GetNumberOfIds();

  //vtkNew(vtkAppendPolyData, appender);

  //for (int i=0; i<numGroups; i++)
  //{
  //  int groupId = centerlineGroupIds->GetId(i);

  //  vtkNew(vtkPolyData, better);
  //  vtkSVGeneralUtils::ThresholdPd(newNormaler->GetOutput(), groupId, groupId, 1, groupIdsArrayName, better);
    // Filter
    vtkNew(vtkSVEdgeWeightedCVT, CVT);

    // OPERATION
    std::cout<<"Performing Operation..."<<endl;
    CVT->SetInputData(newNormaler->GetOutput());
    CVT->SetGenerators(generatorsPd);
    CVT->SetUseCellArray(useCellArray);
    CVT->SetUsePointArray(usePointArray);
    CVT->SetUseGeneratorsArray(useGeneratorsArray);
    CVT->SetNumberOfRings(numberOfRings);
    CVT->SetThreshold(threshold);
    CVT->SetEdgeWeight(edgeWeight);
    CVT->SetUseTransferredPatchesAsThreshold(useTransferredPatchesAsThreshold);
    CVT->SetMaximumNumberOfIterations(maximumNumberOfIterations);
    CVT->SetPatchIdsArrayName(patchIdsArrayName.c_str());
    CVT->SetCVTDataArrayName(newCellArrayName.c_str());
    CVT->SetGeneratorsArrayName(generatorsArrayName.c_str());
    CVT->Update();

  //  appender->AddInputData(CVT->GetOutput());
  //}

  //appender->Update();

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
