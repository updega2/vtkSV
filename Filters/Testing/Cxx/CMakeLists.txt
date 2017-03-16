vtksv_add_test_cxx(${vtk-module}CxxTests tests
  TestHausdorffDistance.cxx,NO_VALID)

vtksv_test_cxx_executable(${vtk-module}CxxTests tests)
