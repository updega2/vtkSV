cmake_minimum_required(VERSION 2.8)

PROJECT(LoftNURBSCurve)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(LoftNURBSCurve TestNURBSCurve.cxx
                              vtkLoftNURBSCurve.cxx
                              vtkControlGrid.cxx
                              vtkNURBSCurve.cxx
                              vtkNURBSUtils.cxx)

add_executable(LoftNURBSSurface TestNURBSSurface.cxx
	                        vtkLoftNURBSSurface.cxx
                                vtkControlGrid.cxx
				vtkNURBSSurface.cxx
                                vtkNURBSUtils.cxx
                                vtkGeneralizedPolyCube.cxx)

if(VTK_LIBRARIES)
  target_link_libraries(LoftNURBSCurve ${VTK_LIBRARIES})
  target_link_libraries(LoftNURBSSurface ${VTK_LIBRARIES})
else()
  target_link_libraries(LoftNURBSCurve vtkHybrid)
  target_link_libraries(LoftNURBSSurface vtkHybrid)
endif()
