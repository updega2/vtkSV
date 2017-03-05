/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMultiplePolyDataIntersectionFilter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkMultiplePolyDataIntersectionFilter
// .SECTION Description
//
// vtkMultiplePolyDataIntersectionFilter 

#ifndef __vtkMultiplePolyDataIntersectionFilter_h
#define __vtkMultiplePolyDataIntersectionFilter_h

#include "vtkFiltersCoreModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"

class vtkCellArray;
class vtkDataArray;
class vtkPoints;
class vtkPolyData;

class VTKFILTERSCORE_EXPORT vtkMultiplePolyDataIntersectionFilter : public vtkPolyDataAlgorithm
{
public:
  static vtkMultiplePolyDataIntersectionFilter *New();

  vtkTypeMacro(vtkMultiplePolyDataIntersectionFilter,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // UserManagedInputs allows the user to set inputs by number instead of
  // using the AddInput/RemoveInput functions. Calls to
  // SetNumberOfInputs/SetInputConnectionByNumber should not be mixed with calls
  // to AddInput/RemoveInput. By default, UserManagedInputs is false.
  vtkSetMacro(UserManagedInputs,int);
  vtkGetMacro(UserManagedInputs,int);
  vtkBooleanMacro(UserManagedInputs,int);

  // Description:
  // Add a dataset to the list of data to append. Should not be
  // used when UserManagedInputs is true, use SetInputByNumber instead.
  void AddInputData(vtkPolyData *);

  // Description:
  // Remove a dataset from the list of data to append. Should not be
  // used when UserManagedInputs is true, use SetInputByNumber (NULL) instead.
  void RemoveInputData(vtkPolyData *);

//BTX
  // Description:
  // Get any input of this filter.
  vtkPolyData *GetInput(int idx);
  vtkPolyData *GetInput() { return this->GetInput( 0 ); };
//ETX

  // Description:
  // Directly set(allocate) number of inputs, should only be used
  // when UserManagedInputs is true.
  void SetNumberOfInputs(int num);

  // Set Nth input, should only be used when UserManagedInputs is true.
  void SetInputConnectionByNumber(int num, vtkAlgorithmOutput *input);
  void SetInputDataByNumber(int num, vtkPolyData *ds);

  // Description:
  // ParallelStreaming is for a particular application.
  // It causes this filter to ask for a different piece
  // from each of its inputs.  If all the inputs are the same,
  // then the output of this append filter is the whole dataset
  // pieced back together.  Duplicate points are create
  // along the seams.  The purpose of this feature is to get
  // data parallelism at a course scale.  Each of the inputs
  // can be generated in a different process at the same time.
  vtkSetMacro(ParallelStreaming, int);
  vtkGetMacro(ParallelStreaming, int);
  vtkBooleanMacro(ParallelStreaming, int);

  // Description:
  // Set/get the boolean determing the output when two objects don't 
  // intersect. With a value of 1, either objects output. With a value of 1,
  // both objects are output. 
  vtkSetMacro(NoIntersectionOutput,int);
  vtkGetMacro(NoIntersectionOutput,int);

  // Description:
  // Set/get boolean to determine whether info such as BoundaryPoints
  // is passed as global information. Local boundary scalars will be
  // passed to the full boolean
  vtkSetMacro(PassInfoAsGlobal,int);
  vtkGetMacro(PassInfoAsGlobal,int);

  // Description:
  // Set/get boolean to determine whether surfaces are given id information
  // the output scalary array will be defined as "ModelFaceId"
  vtkSetMacro(AssignSurfaceIds,int);
  vtkGetMacro(AssignSurfaceIds,int);

//ETX
protected:
  vtkMultiplePolyDataIntersectionFilter();
  ~vtkMultiplePolyDataIntersectionFilter();

  // Flag for selecting parallel streaming behavior
  int ParallelStreaming;

  // Usual data generation method
  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **, vtkInformationVector *);
  virtual int RequestUpdateExtent(vtkInformation *,
                                  vtkInformationVector **, vtkInformationVector *);
  virtual int FillInputPortInformation(int, vtkInformation *);

 private:
  // hide the superclass' AddInput() from the user and the compiler
  void AddInputData(vtkDataObject *)
    { vtkErrorMacro( << "AddInput() must be called with a vtkPolyData not a vtkDataObject."); };

  //User defined booleans for filter management
  int UserManagedInputs;
  int NoIntersectionOutput;
  int PassInfoAsGlobal;
  int AssignSurfaceIds;

  int **IntersectionTable;
  int *inResult;
  vtkPolyData *BooleanObject;

  //Function to build the table defining where intersections occur.
  int BuildIntersectionTable(vtkPolyData* inputs[], int numInputs);
  //Function to run the intersection on intersecting polydatas
  int ExecuteIntersection(vtkPolyData *inputs[],int numInputs,int start);
  //Function to set the boundary point information as global information
  void PreSetGlobalArrays(vtkPolyData *input);
  void PostSetGlobalArrays(int numIntersections);
  //Function to set surface id
  void SetSurfaceId(vtkPolyData *input,int surfaceid);

  //Function to print intersectiontable
  void PrintTable(int numInputs);

private:
  vtkMultiplePolyDataIntersectionFilter(const vtkMultiplePolyDataIntersectionFilter&);  // Not implemented.
  void operator=(const vtkMultiplePolyDataIntersectionFilter&);  // Not implemented.
};

#endif


