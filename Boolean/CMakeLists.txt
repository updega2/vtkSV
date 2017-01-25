cmake_minimum_required(VERSION 2.8)

PROJECT(Intersect)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

option(WITH_MINE "" OFF)

#add_executable(Intersect Intersect.cxx)
	#predicates.cxx
	#vtkIntersectionPolyDataFilterOld.cxx
	#vtkDelaunay2D_60.cxx)

	#add_executable(Boolean Boolean.cxx)
	#	#predicates.cxx
	#	#vtkIntersectionPolyDataFilterOld.cxx
	#	#vtkBooleanOperationPolyDataFilterOld.cxx
	#	#vtkDelaunay2D_60.cxx)

	#add_executable(VesselIntersection VesselIntersection.cxx
	#	#predicates.cxx
	#	#vtkIntersectionPolyDataFilterOld.cxx
	#	#vtkBooleanOperationPolyDataFilterOld.cxx
	#	vtkMultiplePolyDataIntersectionFilter.cxx)
	#			 #vtkDelaunay2D_60.cxx)

add_executable(BooleanOld BooleanOld.cxx
	vtkIntersectionPolyDataFilterOld.cxx
	vtkBooleanOperationPolyDataFilterOld.cxx
	vtkDelaunay2D_60.cxx)

if(VTK_LIBRARIES)
	#target_link_libraries(Intersect ${VTK_LIBRARIES})
	#	target_link_libraries(Boolean ${VTK_LIBRARIES})
	#target_link_libraries(VesselIntersection ${VTK_LIBRARIES})
	target_link_libraries(BooleanOld ${VTK_LIBRARIES})
else()
	#target_link_libraries(Intersect vtkHybrid)
	#	target_link_libraries(Boolean vtkHybrid)
	#target_link_libraries(VesselIntersection vtkHybrid)
	target_link_libraries(BooleanOld vtkHybrid)
endif()

if(WITH_MINE)
	add_definitions(-DUSE_MINE)
endif()
