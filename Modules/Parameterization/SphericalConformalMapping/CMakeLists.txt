cmake_minimum_required(VERSION 2.8)

PROJECT(SphericalConformalMapper)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(SphericalConformalMapper SphericalConformalMapper.cxx
  vtkPlacePointsOnS2.cxx
  vtkSphericalConformalMapper.cxx)

add_executable(SurfaceToSurfaceMapper SurfaceToSurfaceMapper.cxx
  vtkSphericalConformalMapper.cxx
  vtkS2LandmarkMatcher.cxx
  vtkPlacePointsOnS2.cxx
  vtkMapInterpolator.cxx
  vtkSurfaceToSurfaceMapper.cxx
  math.cxx
  sparse_matrix.cxx)

add_executable(PlacePointsOnS2 PlacePointsOnS2.cxx
  vtkPlacePointsOnS2.cxx)

add_executable(CheckRotation CheckRotation.cxx
  vtkCheckRotation.cxx)

add_executable(MapInterpolator MapInterpolator.cxx
  vtkMakeParametricCube.cxx
  vtkMapInterpolator.cxx)

add_executable(MakeParametricCube MakeParametricCube.cxx
  vtkMakeParametricCube.cxx)

add_executable(SliceAndDice SliceAndDice.cxx
  vtkFindSeparateRegions.cxx
  vtkFindClosestGeodesicPoint.cxx
  vtkGeneralizedPolycube.cxx
  vtkGetBoundaryFaces.cxx
  vtkPolyDataSliceAndDiceFilter.cxx
  svGraph.cxx
  svGCell.cxx)

add_executable(PolyDataToNURBS PolyDataToNURBS.cxx
  vtkFindSeparateRegions.cxx
  vtkFindClosestGeodesicPoint.cxx
  vtkGeneralizedPolycube.cxx
  vtkGetBoundaryFaces.cxx
  vtkMapInterpolator.cxx
  vtkPlacePointsOnS2.cxx
  vtkPolyDataSliceAndDiceFilter.cxx
  vtkPolyDataToNURBSFilter.cxx
  vtkSphericalConformalMapper.cxx
  svGraph.cxx
  svGCell.cxx)

if(VTK_LIBRARIES)
  target_link_libraries(SphericalConformalMapper ${VTK_LIBRARIES})
  target_link_libraries(SurfaceToSurfaceMapper ${VTK_LIBRARIES})
  target_link_libraries(CheckRotation ${VTK_LIBRARIES})
  target_link_libraries(PlacePointsOnS2 ${VTK_LIBRARIES})
  target_link_libraries(MapInterpolator ${VTK_LIBRARIES})
  target_link_libraries(MakeParametricCube ${VTK_LIBRARIES})
  target_link_libraries(SliceAndDice ${VTK_LIBRARIES})
  target_link_libraries(PolyDataToNURBS ${VTK_LIBRARIES})
else()
  target_link_libraries(SphericalConformalMapper vtkHybrid)
  target_link_libraries(SurfaceToSurfaceMapper vtkHybrid)
  target_link_libraries(CheckRotation vtkHybrid)
  target_link_libraries(PlacePointsOnS2 vtkHybrid)
  target_link_libraries(MapInterpolator vtkHybrid)
  target_link_libraries(MakeParametricCube vtkHybrid)
  target_link_libraries(PolyDataToNURBS vtkHybrid)
endif()

