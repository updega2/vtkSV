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
#include "vtkDataSetSurfaceFilter.h"
#include "vtkDataWriter.h"
#include "vtkInformation.h"
#include "vtkMassProperties.h"
#include "vtkPointData.h"
#include "vtkPointLocator.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkSplineFilter.h"
#include "vtkThreshold.h"
#include "vtkTriangleFilter.h"

#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkSVGeneralUtils.h"
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
  int appendEndPoints = 0;
  int pickSeedPoints = 0;
  int useAbsoluteMergeDistance = 0;
  int medialEdgeThreshold = 5;
  int absoluteThreshold = 3;
  double mergeDistance = 0.1;
  double radiusMergeRatio = 0.35;
  double relativeThreshold = 0.5;

  std::string radiusArrayName   = "MaximumInscribedSphereRadius";
  std::string seedSelector = "none";

  // argc is the number of strings on the command-line
  //  starting with the program name
  for(iarg=1; iarg<argc; iarg++){
      arglength = strlen(argv[iarg]);
      // replace 0..arglength-1 with argv[iarg]
      tmpstr.replace(0,arglength,argv[iarg],0,arglength);
      if(tmpstr=="-h")                        {RequestedHelp = true;}
      else if(tmpstr=="-input")               {InputProvided = true; inputFilename = argv[++iarg];}
      else if(tmpstr=="-output")              {OutputProvided = true; outputFilename = argv[++iarg];}
      else if(tmpstr=="-radius")              {radiusArrayName = argv[++iarg];}
      else if(tmpstr=="-usevmtk")             {useVmtk = atoi(argv[++iarg]);}
      else if(tmpstr=="-appendendpoints")     {appendEndPoints = atoi(argv[++iarg]);}
      else if(tmpstr=="-seedselector")        {seedSelector = argv[++iarg];}
      else if(tmpstr=="-radiusmergeratio")    {radiusMergeRatio = atof(argv[++iarg]);}
      else if(tmpstr=="-usemergedistance")    {useAbsoluteMergeDistance = atoi(argv[++iarg]);}
      else if(tmpstr=="-mergedistance")       {mergeDistance = atof(argv[++iarg]);}
      else if(tmpstr=="-absolutethreshold")   {absoluteThreshold = atoi(argv[++iarg]);}
      else if(tmpstr=="-relativethreshold")   {relativeThreshold = atof(argv[++iarg]);}
      else if(tmpstr=="-medialedgethreshold") {medialEdgeThreshold = atoi(argv[++iarg]);}
      else {cout << argv[iarg] << " is not a valid argument. Ask for help with -h." << endl; RequestedHelp = true; return EXIT_FAILURE;}
      // reset tmpstr for next argument
      tmpstr.erase(0,arglength);
  }

  if (RequestedHelp || !InputProvided)
  {
    cout << endl;
    cout << "usage:" <<endl;
    cout << "  CenterlinesExtractor -input [Input Filename] -output [Output Filename] ..." << endl;
    cout << endl;
    cout << "COMMAND-LINE ARGUMENT SUMMARY" << endl;
    cout << "  -h                   : Display usage and command-line argument summary"<< endl;
    cout << "  -input               : Input file name (.vtp or .stl)"<< endl;
    cout << "  -output              : Output file name"<< endl;
    cout << "  -radius              : Name on centerlines describing maximum inscribed sphere radius [default MaximumInscribedSphereRadius]"<< endl;
    cout << "  -usevmtk             : Use the vmtk centerlines extractor rather than vtksv [default 0]"<< endl;
    cout << "  -appendendpoints     : Append end points to the end of the centerline paths to touch surface [default 0]"<< endl;
    cout << "  -seedselector        : How to choose source and target seeds if desired, options are none, pickpoints, and openprofiles [default none]"<< endl;
    cout << "  -radiusmergeratio    : When extracting centerline branches, the portion of the radius to use (radius at bifurcation location) to use as the merging distance [default 0.35]"<< endl;
    cout << "  -usemergedistance    : Instead of using a ratio to the radius, use an absolute distance for the merge distance [default 0]" << endl;
    cout << "  -mergedistance       : The merge distance; only used is usemergedistance is on [default 0.1]" << endl;
    cout << "  -absolutethreshold   : The threshold for the absolute persistance metric in the cell thinning process. [default 3]" << endl;
    cout << "  -relativethreshold   : The threshold for the relative persistance metric in the cell thinning process. [default 0.5]" << endl;
    cout << "  -medialedgethreshold : During the cell thinning process, medial edges are extracted. Connected components of more than this threshold value are defined as the medial axis. A higher value will remove small spurious lines. [default 5]" << endl;
    cout << "END COMMAND-LINE ARGUMENT SUMMARY" << endl;
    return EXIT_FAILURE;
  }
  if (!OutputProvided)
  {
    cout << "WARNING: Output Filename not provided, setting output name based on the input filename" <<endl;
    std::string newDirName = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename);
    // Only mac and linux!!!
    system(("mkdir -p "+newDirName).c_str());
    outputFilename = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"_Centerlines.vtp";
  }

  // Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkNew(vtkPolyData, inputPd);
  if (vtkSVIOUtils::ReadInputFile(inputFilename,inputPd) != 1)
    return EXIT_FAILURE;

  std::cout << "Cleaning Surface..." << endl;
  vtkNew(vtkCleanPolyData, cleaner);
  cleaner->SetInputData(inputPd);
  cleaner->Update();

  std::cout << "Triangulating Surface..." << endl;
  vtkNew(vtkTriangleFilter, triangulator);
  triangulator->SetInputData(cleaner->GetOutput());
  triangulator->PassLinesOff();
  triangulator->PassVertsOff();
  triangulator->Update();

  inputPd->DeepCopy(triangulator->GetOutput());

  vtkSVSeedSelector *seedPointPicker;
  vtkNew(vtkIdList, capCenterIds);
  if (seedSelector == "pickpoints")
  {
    int numNonTriangleCells = 0;
    int numNonManifoldEdges = 0;
    int numOpenEdges = 0;
    int surfaceGenus = 0;

    vtkSVGeneralUtils::CheckSurface(inputPd, numNonTriangleCells,
                                    numNonManifoldEdges, numOpenEdges,
                                    surfaceGenus);

    int allGood = 1;
    if (numOpenEdges > 0)
    {
      std::cerr << "The surface has free edges, which means the seedselector needs to be either openprofiles or maxradiusprofile" << endl;
      allGood = 0;
    }
    if (numNonTriangleCells > 0)
    {
      std::cerr << "Surface contains non-triangle cells. Number of non-triangle cells: " << numNonTriangleCells << endl;
      allGood = 0;
    }
    if (numNonManifoldEdges > 0)
    {
      std::cerr << "Surface contains non-manifold edges. Number of non-manifold edges: " << numNonManifoldEdges << endl;
      allGood = 0;
    }
    if (surfaceGenus > 0)
    {
      std::cerr <<  "Surface genus is greater than 0. Surface genus is: " << surfaceGenus << endl;
      allGood = 0;
    }

    if (!allGood)
      return EXIT_FAILURE;

    vtkSVPickPointSeedSelector *pointPicker =
      vtkSVPickPointSeedSelector::New();
    pointPicker->SetInputData(inputPd);
    pointPicker->Update();

    int numSourceSeeds = pointPicker->GetSourceSeedIds()->GetNumberOfIds();
    int numTargetSeeds = pointPicker->GetTargetSeedIds()->GetNumberOfIds();

    std::cout <<"Number of Source Seeds: " << numSourceSeeds << endl;
    std::cout <<"  Source Seeds Ids: " << endl;
    for (int i=0; i<numSourceSeeds; i++)
      std::cout <<" " << pointPicker->GetSourceSeedIds()->GetId(i) << endl;
    std::cout << endl;

    std::cout <<"Number of Target Seeds: " << numTargetSeeds << endl;
    std::cout <<"  Target Seeds Ids: " << endl;
    for (int i=0; i<numTargetSeeds; i++)
      std::cout <<" " << pointPicker->GetTargetSeedIds()->GetId(i) << endl;
    std::cout << endl;

    seedPointPicker = pointPicker;
  }
  else if (seedSelector == "openprofiles" || seedSelector == "maxareaprofile")
  {
    std::cout << "Capping Surface..." << endl;

    vtkNew(vtkvmtkCapPolyData, surfaceCapper);
    surfaceCapper->SetInputData(inputPd);
    surfaceCapper->SetDisplacement(0.0);
    surfaceCapper->SetInPlaneDisplacement(0.0);
    surfaceCapper->SetCellEntityIdsArrayName("CapIds");
    surfaceCapper->Update();

    capCenterIds->DeepCopy(surfaceCapper->GetCapCenterIds());

    std::cout << "Cleaning Capped Surface..." << endl;

    // Clean
    cleaner->SetInputData(surfaceCapper->GetOutput());
    cleaner->Update();

    inputPd->DeepCopy(cleaner->GetOutput());

    // Remove non-triangle cells
    for (int i=0; i<inputPd->GetNumberOfCells(); i++)
    {
      if (inputPd->GetCellType(i) != VTK_TRIANGLE)
      {
        inputPd->DeleteCell(i);
      }
    }

    inputPd->RemoveDeletedCells();
    inputPd->BuildLinks();

    // Reset cap ids if the polydata happened to change
    vtkNew(vtkPointLocator, pointLocator);
    pointLocator->SetDataSet(inputPd);
    pointLocator->BuildLocator();

    int currCapCenterId, newCapCenterId;
    double pt[3];
    for (int i=0; i<capCenterIds->GetNumberOfIds(); i++)
    {
      currCapCenterId = capCenterIds->GetId(i);
      surfaceCapper->GetOutput()->GetPoint(currCapCenterId, pt);

      newCapCenterId = pointLocator->FindClosestPoint(pt);
      capCenterIds->SetId(i, newCapCenterId);
    }

    int numNonTriangleCells = 0;
    int numNonManifoldEdges = 0;
    int numOpenEdges = 0;
    int surfaceGenus = 0;

    vtkSVGeneralUtils::CheckSurface(inputPd, numNonTriangleCells,
                                    numNonManifoldEdges, numOpenEdges,
                                    surfaceGenus);

    int allGood = 1;
    if (numOpenEdges > 0)
    {
      std::cerr << "The surface has free edges, which means the capper did not work correctly, something must be wrong with the surface" << endl;
      allGood = 0;
    }
    if (numNonTriangleCells > 0)
    {
      std::cerr << "Surface contains non-triangle cells. Number of non-triangle cells: " << numNonTriangleCells << endl;
      allGood = 0;
    }
    if (numNonManifoldEdges > 0)
    {
      std::cerr << "Surface contains non-manifold edges. Number of non-manifold edges: " << numNonManifoldEdges << endl;
      allGood = 0;
    }
    if (surfaceGenus > 0)
    {
      std::cerr <<  "Surface genus is greater than 0. Surface genus is: " << surfaceGenus << endl;
      allGood = 0;
    }

    if (!allGood)
      return EXIT_FAILURE;

    if (seedSelector == "openprofiles")
    {
      vtkSVOpenProfilesSeedSelector *openProfilesPicker =
        vtkSVOpenProfilesSeedSelector::New();
      openProfilesPicker->SetInputData(inputPd);
      openProfilesPicker->SetSeedIds(surfaceCapper->GetCapCenterIds());
      openProfilesPicker->Update();

      int numSourceSeeds = openProfilesPicker->GetSourceSeedIds()->GetNumberOfIds();
      int numTargetSeeds = openProfilesPicker->GetTargetSeedIds()->GetNumberOfIds();

      std::cout << "Number of Source Seeds: " << numSourceSeeds << endl;
      std::cout << "  Source Seeds Ids: " << endl;
      for (int i=0; i<numSourceSeeds; i++)
        std::cout << " " << openProfilesPicker->GetSourceSeedIds()->GetId(i) << endl;
      std::cout << endl;

      std::cout << "Number of Target Seeds: " << numTargetSeeds << endl;
      std::cout << "  Target Seeds Ids: " << endl;
      for (int i=0; i<numTargetSeeds; i++)
        std::cout << " " << openProfilesPicker->GetTargetSeedIds()->GetId(i) << endl;
      std::cout << endl;

      seedPointPicker = openProfilesPicker;
    }
  }
  else if (seedSelector == "none")
  {
    std::cout << "No seed points given" << endl;

    int numNonTriangleCells = 0;
    int numNonManifoldEdges = 0;
    int numOpenEdges = 0;
    int surfaceGenus = 0;

    vtkSVGeneralUtils::CheckSurface(inputPd, numNonTriangleCells,
                                    numNonManifoldEdges, numOpenEdges,
                                    surfaceGenus);

    int allGood = 1;

    if (numOpenEdges > 0)
    {
      std::cerr << "The surface has free edges, which means the seedselector needs to be either openprofiles or maxradiusprofile " << endl;
      allGood = 0;
    }
    if (numNonTriangleCells > 0)
    {
      std::cerr << "Surface contains non-triangle cells. Number of non-triangle cells: " << numNonTriangleCells << endl;
      allGood = 0;
    }
    if (numNonManifoldEdges > 0)
    {
      std::cerr << "Surface contains non-manifold edges. Number of non-manifold edges: " << numNonManifoldEdges << endl;
      allGood = 0;
    }
    if (surfaceGenus > 0)
    {
      std::cerr <<  "Surface genus is greater than 0. Surface genus is: " << surfaceGenus << endl;
      allGood = 0;
    }

    if (!allGood)
      return EXIT_FAILURE;
  }
  else
  {
    std::cerr << "Incorrect seedselector given, must be none, pickpoints or openprofiles" << endl;
  }

  // Filter
  if (useVmtk)
  {
    vtkNew(vtkvmtkPolyDataCenterlines, CenterlineFilter);

    //OPERATION
    std::cout<<"Getting Centerlines..."<<endl;
    if (seedSelector == "pickpoints" || seedSelector == "openprofiles")
    {
      CenterlineFilter->SetSourceSeedIds(seedPointPicker->GetSourceSeedIds());
      CenterlineFilter->SetTargetSeedIds(seedPointPicker->GetTargetSeedIds());
      if (seedSelector == "openprofiles")
        CenterlineFilter->SetCapCenterIds(capCenterIds);
    }
    else if (seedSelector == "maxareaprofile")
    {
      vtkNew(vtkThreshold, capThresholder);
      capThresholder->SetInputData(inputPd);
      capThresholder->SetInputArrayToProcess(0, 0, 0, 1, "CapIds");

      vtkNew(vtkDataSetSurfaceFilter, surfacer);

      vtkNew(vtkMassProperties, measurer);

      double maxArea = -1.0;
      int maxAreaCapId = -1;

      int numCaps = capCenterIds->GetNumberOfIds();
      for (int i=0; i<numCaps; i++)
      {
        capThresholder->ThresholdBetween(i+2, i+2);
        capThresholder->Update();

        surfacer->SetInputData(capThresholder->GetOutput());
        surfacer->Update();

        measurer->SetInputData(surfacer->GetOutput());
        measurer->Update();

        if (measurer->GetSurfaceArea() > maxArea)
        {
          maxArea = measurer->GetSurfaceArea();
          maxAreaCapId = i;
        }
      }

      vtkNew(vtkIdList, newSourceId);
      newSourceId->InsertNextId(maxAreaCapId);

      vtkNew(vtkIdList, newTargetIds);
      int testId;
      for (int i=0; i<numCaps; i++)
      {
        if (i != maxAreaCapId)
        {
          newTargetIds->InsertNextId(i);
        }
      }

      CenterlineFilter->SetSourceSeedIds(newSourceId);
      CenterlineFilter->SetTargetSeedIds(newTargetIds);
      CenterlineFilter->SetCapCenterIds(capCenterIds);
    }
    CenterlineFilter->SetInputData(inputPd);
    CenterlineFilter->SetRadiusArrayName(radiusArrayName.c_str());
    CenterlineFilter->SetCostFunction("1/R");
    CenterlineFilter->SetSimplifyVoronoi(0);
    CenterlineFilter->SetAppendEndPointsToCenterlines(appendEndPoints);
    CenterlineFilter->SetCenterlineResampling(0);
    CenterlineFilter->SetResamplingStepLength(0.0);
    CenterlineFilter->DebugOn();
    CenterlineFilter->Update();

    std::cout<<"Done"<<endl;

    if (seedSelector == "pickpoints" || seedSelector == "openprofiles")
      seedPointPicker->Delete();

    //Write Files
    std::cout<<"Writing Files..."<<endl;
    vtkSVIOUtils::WriteVTPFile(outputFilename, CenterlineFilter->GetOutput(0));

    if (seedSelector == "openprofiles" || seedSelector == "maxareaprofile")
    {
      if (CenterlineFilter->GetOutput()->GetNumberOfLines() != capCenterIds->GetNumberOfIds() - 1)
      {
        std::cerr << "Incorrect number of lines found: "<< CenterlineFilter->GetOutput()->GetNumberOfLines() << ". But should be " << capCenterIds->GetNumberOfIds() -1 << endl;
        return EXIT_FAILURE;
      }
    }

  }
  else
  {
    vtkNew(vtkSVCenterlines, CenterlineFilter);

    //OPERATION
    std::cout<<"Getting Centerlines..."<<endl;
    if (seedSelector == "pickpoints" || seedSelector == "openprofiles")
    {
      CenterlineFilter->SetSourceSeedIds(seedPointPicker->GetSourceSeedIds());
      CenterlineFilter->SetTargetSeedIds(seedPointPicker->GetTargetSeedIds());
      if (seedSelector == "openprofiles")
        CenterlineFilter->SetCapCenterIds(capCenterIds);
    }
    else if (seedSelector == "maxareaprofile")
    {
      vtkNew(vtkThreshold, capThresholder);
      capThresholder->SetInputData(inputPd);
      capThresholder->SetInputArrayToProcess(0, 0, 0, 1, "CapIds");

      vtkNew(vtkDataSetSurfaceFilter, surfacer);

      vtkNew(vtkMassProperties, measurer);

      double maxArea = -1.0;
      int maxAreaCapId = -1;

      int numCaps = capCenterIds->GetNumberOfIds();
      for (int i=0; i<numCaps; i++)
      {
        capThresholder->ThresholdBetween(i+2, i+2);
        capThresholder->Update();

        surfacer->SetInputData(capThresholder->GetOutput());
        surfacer->Update();

        measurer->SetInputData(surfacer->GetOutput());
        measurer->Update();

        if (measurer->GetSurfaceArea() > maxArea)
        {
          maxArea = measurer->GetSurfaceArea();
          maxAreaCapId = i;
        }
      }

      vtkNew(vtkIdList, newSourceId);
      newSourceId->InsertNextId(maxAreaCapId);

      vtkNew(vtkIdList, newTargetIds);
      int testId;
      for (int i=0; i<numCaps; i++)
      {
        if (i != maxAreaCapId)
        {
          newTargetIds->InsertNextId(i);
        }
      }

      CenterlineFilter->SetSourceSeedIds(newSourceId);
      CenterlineFilter->SetTargetSeedIds(newTargetIds);
      CenterlineFilter->SetCapCenterIds(capCenterIds);
    }
    CenterlineFilter->SetInputData(inputPd);
    CenterlineFilter->SetRadiusArrayName(radiusArrayName.c_str());
    CenterlineFilter->SetCostFunction("1/R");
    CenterlineFilter->SetSimplifyVoronoi(0);
    CenterlineFilter->SetAppendEndPointsToCenterlines(appendEndPoints);
    CenterlineFilter->SetCenterlineResampling(0);
    CenterlineFilter->SetResamplingStepLength(0.0);
    CenterlineFilter->SetAbsoluteThreshold(absoluteThreshold);
    CenterlineFilter->SetRelativeThreshold(relativeThreshold);
    CenterlineFilter->SetMedialEdgeThreshold(medialEdgeThreshold);
    //CenterlineFilter->SetProcessCenterlinesIntoTree(0);
    CenterlineFilter->DebugOn();
    CenterlineFilter->Update();

    std::cout<<"Done"<<endl;

    if (seedSelector == "pickpoints" || seedSelector == "openprofiles")
      seedPointPicker->Delete();

    //Write Files
    std::cout<<"Writing Files..."<<endl;
    vtkSVIOUtils::WriteVTPFile(outputFilename, CenterlineFilter->GetOutput(0));

    if (seedSelector == "openprofiles" || seedSelector == "maxareaprofile")
    {
      if (CenterlineFilter->GetOutput()->GetNumberOfLines() != capCenterIds->GetNumberOfIds() - 1)
      {
        std::cerr << "Incorrect number of lines found: "<< CenterlineFilter->GetOutput()->GetNumberOfLines() << ". But should be " << capCenterIds->GetNumberOfIds() -1 << endl;
        return EXIT_FAILURE;
      }
    }

  }


  //Exit the program without errors
  return EXIT_SUCCESS;
}
