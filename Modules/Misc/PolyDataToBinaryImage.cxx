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
#include "vtkImageStencil.h"
#include "vtkInformation.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataToImageStencil.h"
#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkSVGetBoundaryFaces.h"

int main(int argc, char *argv[])
{
  // BEGIN PROCESSING COMMAND-LINE ARGUMENTS
  // Assume no options specified at command line
  bool RequestedHelp          = false;
  bool InputProvided          = false;
  bool OutputProvided         = false;

  // Variables used in processing the commandline
  int iarg, arglength;
  std::string tmpstr;

  // Filenames
  std::string inputFilename;
  std::string outputFilename;

  // Default values for options
  double gridSpacing = 0.5;

  // argc is the number of strings on the command-line
  //  starting with the program name
  for(iarg=1; iarg<argc; iarg++){
      arglength = strlen(argv[iarg]);
      // replace 0..arglength-1 with argv[iarg]
      tmpstr.replace(0,arglength,argv[iarg],0,arglength);
      if(tmpstr=="-h")             {RequestedHelp = true;}
      else if(tmpstr=="-input")    {InputProvided = true; inputFilename = argv[++iarg];}
      else if(tmpstr=="-output")   {OutputProvided = true; outputFilename = argv[++iarg];}
      else if(tmpstr=="-spacing")  {gridSpacing = atof(argv[++iarg]);}
      else {cout << argv[iarg] << " is not a valid argument. Ask for help with -h." << endl; RequestedHelp = true; return EXIT_FAILURE;}
      // reset tmpstr for next argument
      tmpstr.erase(0,arglength);
  }

  if (RequestedHelp || !InputProvided)
  {
    cout << endl;
    cout << "usage:" <<endl;
    cout << "  PolyDataToBinaryImage -input [Input Filename (.vtp)] -output [Output Filename (.mhd)] -spacing [Image spacing] ..." << endl;
    cout << endl;
    cout << "COMMAND-LINE ARGUMENT SUMMARY" << endl;
    cout << "  -h        : Display usage and command-line argument summary"<< endl;
    cout << "  -input    : Input file name (.vtp or .stl)"<< endl;
    cout << "  -output   : Output file name (.mhd)"<< endl;
    cout << "  -spacing  : Spacing for output image [default 0.5]"<< endl;
    cout << "END COMMAND-LINE ARGUMENT SUMMARY" << endl;
    return EXIT_FAILURE;
  }
  if (!OutputProvided)
  {
    cout << "WARNING: Output Filename not provided, setting output name based on the input filename" <<endl;
    std::string newDirName = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename);
    // Only mac and linux!!!
    system(("mkdir -p "+newDirName).c_str());
    outputFilename = vtkSVIOUtils::GetPath(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"/"+vtkSVIOUtils::GetRawName(inputFilename)+"_Binary_Image.mhd";
  }

  // Call Function to Read File
  std::cout<<"Reading Files..."<<endl;
  vtkNew(vtkPolyData, inputPd);
  if (vtkSVIOUtils::ReadInputFile(inputFilename,inputPd) != 1)
    return EXIT_FAILURE;

  // Clean input
  vtkNew(vtkCleanPolyData, cleaner);
  cleaner->SetInputData(inputPd);
  cleaner->Update();

  inputPd->DeepCopy(cleaner->GetOutput());

  // Bounds
  double bounds[6];
  inputPd->GetBounds(bounds);

  // Spacing
  double spacing[3];
  spacing[0] = gridSpacing;
  spacing[1] = gridSpacing;
  spacing[2] = gridSpacing;

  // Image
  vtkNew(vtkImageData, whiteImage);
  whiteImage->SetSpacing(spacing);

  // Dim
  int dims[3];
  dims[0] = (int) ceil((bounds[1] - bounds[0]) / spacing[0]);
  dims[1] = (int) ceil((bounds[3] - bounds[2]) / spacing[1]);
  dims[2] = (int) ceil((bounds[5] - bounds[4]) / spacing[2]);
  whiteImage->SetDimensions(dims);
  whiteImage->SetExtent(0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1);

  // Origin
  double origin[3];
  origin[0] = (bounds[0] + spacing[0] / 2.);
  origin[1] = (bounds[2] + spacing[1] / 2.);
  origin[2] = (bounds[4] + spacing[2] / 2.);

  // Scalars
  std::cout<<"Allocating scalars..."<<endl;
  whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

  int inval = 255;
  int outval = 0;

  std::cout<<"Filling component..."<<endl;
  int count = whiteImage->GetNumberOfPoints();
  whiteImage->GetPointData()->GetScalars()->FillComponent(0, inval);
  //for (int i=0; i<count; i++)
  //{
  //  whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);
  //}

  std::cout<<"Running Filter..."<<endl;
  // PolyData to image
  vtkNew(vtkPolyDataToImageStencil, pol2Stenc);
  pol2Stenc->SetInputData(inputPd);
  pol2Stenc->SetOutputOrigin(origin);
  pol2Stenc->SetOutputSpacing(spacing);
  pol2Stenc->SetOutputWholeExtent(whiteImage->GetExtent());
  pol2Stenc->Update();

  // Cut white image and set background
  vtkNew(vtkImageStencil, imgStenc);
  imgStenc->SetInputData(whiteImage);
  imgStenc->SetStencilConnection(pol2Stenc->GetOutputPort());
  imgStenc->ReverseStencilOff();
  imgStenc->SetBackgroundValue(outval);
  imgStenc->Update();

  // Write out image
  std::cout<<"Writing Files..."<<endl;
  vtkSVIOUtils::WriteMHDFile(outputFilename, imgStenc->GetOutput());

  //Exit the program without errors
  return EXIT_SUCCESS;
}
