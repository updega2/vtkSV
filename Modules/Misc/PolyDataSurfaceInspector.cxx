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
#include "vtkInformation.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkSVPolyDataSurfaceInspector.h"
#include "vtkTriangleFilter.h"

#if !defined(_WIN32) || defined(__CYGWIN__)
# include <unistd.h> /* unlink */
#else
# include <io.h> /* unlink */
#endif


int main(int argc, char *argv[])
{
  // BEGIN PROCESSING COMMAND-LINE ARGUMENTS
  // Assume no options specified at command line
  bool RequestedHelp   = false;
  bool InputProvided   = false;
  bool OutputProvided  = false;

  // Variables used in processing the commandline
  int iarg, arglength;
  std::string tmpstr;

  // Filenames
  std::string inputFilename;
  std::string outputFilename;

  // Default values for options

  // argc is the number of strings on the command-line
  //  starting with the program name
  for(iarg=1; iarg<argc; iarg++){
      arglength = strlen(argv[iarg]);
      // replace 0..arglength-1 with argv[iarg]
      tmpstr.replace(0,arglength,argv[iarg],0,arglength);
      if(tmpstr=="-h")                      {RequestedHelp = true;}
      else if(tmpstr=="-input")             {InputProvided = true; inputFilename = argv[++iarg];}
      else if(tmpstr=="-output")            {OutputProvided = true; outputFilename = argv[++iarg];}
      else {cout << argv[iarg] << " is not a valid argument. Ask for help with -h." << endl; RequestedHelp = true; return EXIT_FAILURE;}
      // reset tmpstr for next argument
      tmpstr.erase(0,arglength);
  }

  if (RequestedHelp || !InputProvided)
  {
    cout << endl;
    cout << "usage:" <<endl;
    cout << "  PolyDataSurfaceInspector -input [Input Filename] -output [Output Filename]..." << endl;
    cout << endl;
    cout << "COMMAND-LINE ARGUMENT SUMMARY" << endl;
    cout << "  -h                  : Display usage and command-line argument summary"<< endl;
    cout << "  -input              : Input file name (.vtp or .stl)"<< endl;
    cout << "  -output             : Output file information (.txt)"<< endl;
    cout << "END COMMAND-LINE ARGUMENT SUMMARY" << endl;
    return EXIT_FAILURE;
  }

  if (!OutputProvided)
  {
    cout << "WARNING: Output Filename not provided, setting output name based on the input filename" <<endl;
    std::string newDirName = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename);
    // Only mac and linux!!!
    system(("mkdir -p "+newDirName).c_str());
    outputFilename = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"_ModelInfo.txt";
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

  // Filter
  vtkNew(vtkSVPolyDataSurfaceInspector, PolyDataSurfaceInspector);

  //OPERATION
  std::cout<<"Performing Operation..."<<endl;
  PolyDataSurfaceInspector->SetInputData(inputPd);
  PolyDataSurfaceInspector->CheckNumberOfConnectedRegionsOn();
  PolyDataSurfaceInspector->CheckNumberOfHolesOn();
  PolyDataSurfaceInspector->Update();

  if (PolyDataSurfaceInspector->GetErrorCode() != 0)
  {
    std::cerr << "Error in filter" << endl;
    return EXIT_FAILURE;
  }

  int numCells  = PolyDataSurfaceInspector->GetNumberOfElements();
  int numPoints = PolyDataSurfaceInspector->GetNumberOfPoints();
  int numEdges  = PolyDataSurfaceInspector->GetNumberOfEdges();
  int numNonTriangleCells = PolyDataSurfaceInspector->GetNumberOfNonTriangularElements();
  int numNonManifoldEdges = PolyDataSurfaceInspector->GetNumberOfNonManifoldEdges();
  int numOpenEdges        = PolyDataSurfaceInspector->GetNumberOfOpenEdges();
  int surfaceGenus        = PolyDataSurfaceInspector->GetSurfaceGenus();
  int numConnectedRegions = PolyDataSurfaceInspector->GetNumberOfConnectedRegions();
  int numHoles            = PolyDataSurfaceInspector->GetNumberOfHoles();

  FILE *fp;

  if ((fp = fopen(outputFilename.c_str(), "w")) == NULL)
  {
    std::cerr<<"Could not open file"<<endl;
    return EXIT_FAILURE;
  }

  fprintf(stdout, "%s MODEL STATISTICS\n", vtkSVIOUtils::GetRawName(inputFilename).c_str());
  fprintf(stdout, "\n");
  fprintf(stdout, "  num_cells:              %d\n", numCells);
  fprintf(stdout, "  num_points:             %d\n", numPoints);
  fprintf(stdout, "  num_edges:              %d\n", numEdges);
  fprintf(stdout, "  num_non_triangle_cells: %d\n", numNonTriangleCells);
  fprintf(stdout, "  num_non_manifold_edges: %d\n", numNonManifoldEdges);
  fprintf(stdout, "  num_open_edges:         %d\n", numOpenEdges);
  fprintf(stdout, "  num_connected_regions:  %d\n", numConnectedRegions);
  fprintf(stdout, "  num_holes:              %d\n", numHoles);
  fprintf(stdout, "  genus:                  %d\n", surfaceGenus);

  fprintf(fp, "%s MODEL STATISTICS\n", vtkSVIOUtils::GetRawName(inputFilename).c_str());
  fprintf(fp, "\n");
  fprintf(fp, "  num_cells:              %d\n", numCells);
  fprintf(fp, "  num_points:             %d\n", numPoints);
  fprintf(fp, "  num_edges:              %d\n", numEdges);
  fprintf(fp, "  num_non_triangle_cells: %d\n", numNonTriangleCells);
  fprintf(fp, "  num_non_manifold_edges: %d\n", numNonManifoldEdges);
  fprintf(fp, "  num_open_edges:         %d\n", numOpenEdges);
  fprintf(fp, "  num_connected_regions:  %d\n", numConnectedRegions);
  fprintf(fp, "  num_holes:              %d\n", numHoles);
  fprintf(fp, "  genus:                  %d\n", surfaceGenus);

  fclose (fp);

  //// Make sure that the input is triangulated, watertight, genus 0 surface
  //if (numNonTriangleCells > 0)
  //{
  //  return EXIT_FAILURE;
  //}
  //if (numNonManifoldEdges > 0)
  //{
  //  return EXIT_FAILURE;
  //}
  //if (numOpenEdges > 0)
  //{
  //  return EXIT_FAILURE;
  //}
  //if (surfaceGenus > 0)
  //{
  //  return EXIT_FAILURE;
  //}
  // ------------------------------------------------------------------------

  fprintf(stdout,"RETURNING HERE\n");
  //Exit the program without errors
  return EXIT_SUCCESS;
}
