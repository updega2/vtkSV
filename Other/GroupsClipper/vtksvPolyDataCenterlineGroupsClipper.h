/*=========================================================================

Program:   VMTK
Module:    $RCSfile: vtksvPolyDataCenterlineGroupsClipper.h,v $
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
  // .NAME vtksvPolyDataCenterlineGroupsClipper - .
  // .SECTION Description
  // ...

#ifndef __vtksvPolyDataCenterlineGroupsClipper_h
#define __vtksvPolyDataCenterlineGroupsClipper_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyData.h"
#include "vtkIdList.h"

class vtksvPolyDataCenterlineGroupsClipper : public vtkPolyDataAlgorithm
{
  public:
  vtkTypeMacro(vtksvPolyDataCenterlineGroupsClipper,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtksvPolyDataCenterlineGroupsClipper *New();

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
  vtksvPolyDataCenterlineGroupsClipper();
  ~vtksvPolyDataCenterlineGroupsClipper();

  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int ReplaceDataOnCells(vtkPointSet *pointset,
                         const int replaceVal, const int currVal,
                         const std::string &arrName);
  int FindGroupSeparatingPoints(vtkPolyData *pd,
                                vtkIdList *separateIds);
  int GetPointGroups(vtkPolyData *pd, std::string arrayName,
                     const int pointId, vtkIdList *groupIds);

  int SplitGroups(vtkPolyData *pd, vtkIdList *separateIds, vtkIdList *newPointIds,
                  vtkDoubleArray *averageDistances);
  int AddNewCriticalPoint(vtkPolyData *pd,
                          vtkPoints *separatePoints,
                          vtkIdList *groupIds,
                          vtkIdList *newPointIds);
  int FillGroups(vtkPolyData *pd, vtkIdList *newPointIds, vtkDoubleArray *averageDistances);
  int FillRegionGroups(vtkPolyData *pd, const int pointId, const double rayDist);
  int FillCellGroups(vtkPolyData *pd, const int cellId, const int pointId,
                     vtkIdList *cellList,
                     vtkIdList *groupList);
  int GetNewPointSpecial(vtkPolyData *pd, double centroid[3],
                         int &newCell, double newPt[3]);
  int ThresholdPd(vtkPolyData *pd, int minVal, int maxVal, int dataType,
                  std::string arrayName, vtkPolyData *returnPd);
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
  vtksvPolyDataCenterlineGroupsClipper(const vtksvPolyDataCenterlineGroupsClipper&);  // Not implemented.
  void operator=(const vtksvPolyDataCenterlineGroupsClipper&);  // Not implemented.
};

#endif
