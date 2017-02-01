cmake_minimum_required(VERSION 2.8)

PROJECT(PlanarMapper)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(PlanarMapper PlanarMapper.cxx
  vtkBoundaryMapper.cxx
  vtkSquareBoundaryMapper.cxx
  vtkPlanarMapper.cxx
  svMath.cxx
  sparse_matrix.cxx)

add_executable(SuperSquarePlanarMapper SuperSquarePlanarMapper.cxx
  vtkBoundaryMapper.cxx
  vtkSuperSquareBoundaryMapper.cxx
  vtkPlanarMapper.cxx
  svMath.cxx
  sparse_matrix.cxx)

add_executable(SliceAndDice SliceAndDice.cxx
  vtkFindSeparateRegions.cxx
  vtkFindGeodesicPath.cxx
  vtkGeneralizedPolycube.cxx
  vtkGetBoundaryFaces.cxx
  vtkPolyDataSliceAndDiceFilter.cxx
  svGraph.cxx
  svGCell.cxx)

add_executable(PolyDataToNURBS PolyDataToNURBS.cxx
  vtkFindSeparateRegions.cxx
  vtkFindGeodesicPath.cxx
  vtkGeneralizedPolycube.cxx
  vtkGetBoundaryFaces.cxx
  vtkMapInterpolator.cxx
  vtkPolyDataSliceAndDiceFilter.cxx
  vtkPolyDataToNURBSFilter.cxx
  vtkPullApartPolyData.cxx
  vtkPlanarMapper.cxx
  vtkBoundaryMapper.cxx
  vtkSuperSquareBoundaryMapper.cxx
  sparse_matrix.cxx
  svMath.cxx
  svGraph.cxx
  svGCell.cxx)

if(VTK_LIBRARIES)
  target_link_libraries(PlanarMapper ${VTK_LIBRARIES})
  target_link_libraries(SuperSquarePlanarMapper ${VTK_LIBRARIES})
  target_link_libraries(SliceAndDice ${VTK_LIBRARIES})
  target_link_libraries(PolyDataToNURBS ${VTK_LIBRARIES})
else()
  target_link_libraries(PlanarMapper vtkHybrid)
  target_link_libraries(SuperSquarePlanarMapper vtkHybrid)
  target_link_libraries(SliceAndDice ${VTK_LIBRARIES})
  target_link_libraries(PolyDataToNURBS ${VTK_LIBRARIES})
endif()

