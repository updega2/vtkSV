/*=========================================================================

Program:   VMTK
Module:    $RCSfile: vtkSVPolyDataCenterlineGroupsClipper.h,v $
Language:  C++
Date:      $Date: 2006/04/06 16:46:43 $
Version:   $Revision: 1.7 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
  // .NAME vtkSVPolyDataCenterlineGroupsClipper - .
  // .SECTION Description
  // ...

#ifndef __vtkSVPolyDataCenterlineGroupsClipper_h
#define __vtkSVPolyDataCenterlineGroupsClipper_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyData.h"
#include "vtkIdList.h"

class vtkSVPolyDataCenterlineGroupsClipper : public vtkPolyDataAlgorithm
{
  public:
  vtkTypeMacro(vtkSVPolyDataCenterlineGroupsClipper,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkSVPolyDataCenterlineGroupsClipper *New();

  vtkSetObjectMacro(Centerlines,vtkPolyData);
  vtkGetObjectMacro(Centerlines,vtkPolyData);

  vtkSetObjectMacro(CenterlineGroupIds,vtkIdList);
  vtkGetObjectMacro(CenterlineGroupIds,vtkIdList);

  vtkSetStringMacro(CenterlineGroupIdsArrayName);
  vtkGetStringMacro(CenterlineGroupIdsArrayName);

  vtkSetStringMacro(CenterlineRadiusArrayName);
  vtkGetStringMacro(CenterlineRadiusArrayName);

  vtkSetStringMacro(GroupIdsArrayName);
  vtkGetStringMacro(GroupIdsArrayName);

  vtkSetStringMacro(BlankingArrayName);
  vtkGetStringMacro(BlankingArrayName);

  vtkSetMacro(ClipAllCenterlineGroupIds,int);
  vtkGetMacro(ClipAllCenterlineGroupIds,int);
  vtkBooleanMacro(ClipAllCenterlineGroupIds,int);

  vtkSetMacro(GenerateClippedOutput,int);
  vtkGetMacro(GenerateClippedOutput,int);
  vtkBooleanMacro(GenerateClippedOutput,int);

  vtkPolyData* GetClippedOutput();

  vtkSetMacro(CutoffRadiusFactor,double);
  vtkGetMacro(CutoffRadiusFactor,double);

  vtkSetMacro(ClipValue,double);
  vtkGetMacro(ClipValue,double);

  vtkSetMacro(UseRadiusInformation,int);
  vtkGetMacro(UseRadiusInformation,int);
  vtkBooleanMacro(UseRadiusInformation,int);

  protected:
  vtkSVPolyDataCenterlineGroupsClipper();
  ~vtkSVPolyDataCenterlineGroupsClipper();

  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int ReplaceDataOnCells(vtkPointSet *pointset,
                         const int replaceVal, const int currVal,
                         const std::string &arrName);
  int FindGroupSeparatingPoints(vtkPolyData *pd,
                                vtkIdList *separateIds);
  int GetPointGroups(vtkPolyData *pd, std::string arrayName,
                     const int pointId, vtkIdList *groupIds);

  int SplitGroups(vtkPolyData *pd, vtkIdList *separateIds, vtkPoints *newPoints);
  int FillGroups(vtkPolyData *pd, vtkPoints *newPoints);
  int FillRegionGroups(vtkPolyData *pd, vtkPolyData *boundary, double center[3]);
  int ThresholdPd(vtkPolyData *pd, int minVal, int maxVal, int dataType,
                  std::string arrayName, vtkPolyData *returnPd);
  int ComputeMassCenter(vtkPolyData *pd, double massCenter[3]);
  int GetClosestPointConnectedRegion(vtkPolyData *inPd,
                                     double origin[3],
                                     vtkPolyData *outPd);
  vtkPolyData* Centerlines;

  vtkIdList* CenterlineGroupIds;

  char* CenterlineGroupIdsArrayName;
  char* CenterlineRadiusArrayName;

  char* GroupIdsArrayName;
  char* BlankingArrayName;

  int ClipAllCenterlineGroupIds;
  double CutoffRadiusFactor;
  double ClipValue;

  int GenerateClippedOutput;

  int UseRadiusInformation;

  private:
  vtkSVPolyDataCenterlineGroupsClipper(const vtkSVPolyDataCenterlineGroupsClipper&);  // Not implemented.
  void operator=(const vtkSVPolyDataCenterlineGroupsClipper&);  // Not implemented.
};

#endif
