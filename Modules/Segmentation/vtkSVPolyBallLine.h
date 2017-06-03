/*=========================================================================

Program:   VMTK
Module:    $RCSfile: vtkSVPolyBallLine.h,v $
Language:  C++
Date:      $Date: 2006/04/06 16:46:43 $
Version:   $Revision: 1.4 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

/**
 *  \class vtkSVPolyBallLine
 */

#ifndef vtkSVPolyBallLine_h
#define vtkSVPolyBallLine_h

#include "vtkImplicitFunction.h"
#include "vtkSVSegmentationModule.h" // For export
#include "vtkPolyData.h"
#include "vtkIdList.h"

class VTKSVSEGMENTATION_EXPORT vtkSVPolyBallLine : public vtkImplicitFunction
{
public:

  static vtkSVPolyBallLine *New();
  vtkTypeMacro(vtkSVPolyBallLine,vtkImplicitFunction);
  void PrintSelf(ostream& os, vtkIndent indent);

  //@{
  /// \brief Evaluate polyball.
  double EvaluateFunction(double x[3]);
  double EvaluateFunction(double x, double y, double z)
  {return this->vtkImplicitFunction::EvaluateFunction(x, y, z); } ;
  //@}

  //@{
  /// \brief Evaluate polyball gradient.
  void EvaluateGradient(double x[3], double n[3]);
  //@}
  //
  //@{
  /// \brief Set / get input poly data.
  vtkSetVector3Macro(LastLocalCoordX, double);
  vtkGetVector3Macro(LastLocalCoordX, double);
  vtkSetVector3Macro(LastLocalCoordY, double);
  vtkGetVector3Macro(LastLocalCoordY, double);
  vtkSetVector3Macro(LastLocalCoordZ, double);
  vtkGetVector3Macro(LastLocalCoordZ, double);
  //@}

  //@{
  /// \brief Set / get input poly data.
  vtkSetVector3Macro(PointNormal, double);
  vtkGetVector3Macro(PointNormal, double);
  //@}

  //@{
  /// \brief Set / get input poly data.
  vtkSetObjectMacro(Input,vtkPolyData);
  vtkGetObjectMacro(Input,vtkPolyData);
  //@}

  //@{
  /// \brief Set / get input cell ids used for the function.
  vtkSetObjectMacro(InputCellIds,vtkIdList);
  vtkGetObjectMacro(InputCellIds,vtkIdList);
  //@}

  //@{
  /// \brief Set / get a single input cell id used for the function.
  vtkSetMacro(InputCellId,vtkIdType);
  vtkGetMacro(InputCellId,vtkIdType);
  //@}

  //@{
  /// \brief Set / get poly ball radius array name.
  vtkSetStringMacro(PolyBallRadiusArrayName);
  vtkGetStringMacro(PolyBallRadiusArrayName);
  //@}

  //@{
  /// \brief Set / get poly ball radius array name.
  vtkSetStringMacro(LocalCoordinatesArrayName);
  vtkGetStringMacro(LocalCoordinatesArrayName);
  //@}

  //@{
  /// \brief Get the id of the last nearest poly ball center.
  vtkGetMacro(LastPolyBallCellId,vtkIdType);
  vtkGetMacro(LastPolyBallCellSubId,vtkIdType);
  vtkGetMacro(LastPolyBallCellPCoord,double);
  vtkGetVectorMacro(LastPolyBallCenter,double,3);
  vtkGetMacro(LastPolyBallCenterRadius,double);
  //@}

  //@{
  /// \brief Use radius information
  vtkSetMacro(UseRadiusInformation,int);
  vtkGetMacro(UseRadiusInformation,int);
  vtkBooleanMacro(UseRadiusInformation,int);
  //@}

  //@{
  /// \brief Use given normal for directional help
  vtkSetMacro(UsePointNormal,int);
  vtkGetMacro(UsePointNormal,int);
  vtkBooleanMacro(UsePointNormal,int);
  //@}

  //@{
  /// \brief Use given normal for directional help
  vtkSetMacro(UseRadiusWeighting,int);
  vtkGetMacro(UseRadiusWeighting,int);
  vtkBooleanMacro(UseRadiusWeighting,int);
  //@}

  //@{
  /// \brief Use radius information
  vtkSetMacro(UseLocalCoordinates,int);
  vtkGetMacro(UseLocalCoordinates,int);
  vtkBooleanMacro(UseLocalCoordinates,int);
  //@}

  //@{
  /// \brief Remove end points from centerline possible points to check
  vtkSetMacro(RemoveEndPoints,int);
  vtkGetMacro(RemoveEndPoints,int);
  vtkBooleanMacro(RemoveEndPoints,int);
  //@}
  //@{
  /// \brief Remove end points from centerline possible points to check
  vtkSetMacro(ControlEndPoints,int);
  vtkGetMacro(ControlEndPoints,int);
  vtkBooleanMacro(ControlEndPoints,int);
  //@}

  static double ComplexDot(double x[4], double y[4]);

protected:
  vtkSVPolyBallLine();
  ~vtkSVPolyBallLine();

  char* PolyBallRadiusArrayName;
  char* LocalCoordinatesArrayName;

  vtkPolyData* Input;
  vtkIdList* InputCellIds;
  vtkIdType InputCellId;

  vtkIdType LastPolyBallCellId;
  vtkIdType LastPolyBallCellSubId;

  int UseRadiusInformation;
  int UsePointNormal;
  int UseRadiusWeighting;
  int UseLocalCoordinates;
  int RemoveEndPoints;
  int ControlEndPoints;

  double LastPolyBallCellPCoord;
  double LastPolyBallCenter[3];
  double LastPolyBallCenterRadius;
  double LastLocalCoordX[3];
  double LastLocalCoordY[3];
  double LastLocalCoordZ[3];
  double PointNormal[3];

private:
  vtkSVPolyBallLine(const vtkSVPolyBallLine&);  // Not implemented.
  void operator=(const vtkSVPolyBallLine&);  // Not implemented.
};

#endif


