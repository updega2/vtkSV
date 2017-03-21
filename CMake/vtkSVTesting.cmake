# Copyright (c) 2014-2015 The Regents of the University of California.
# All Rights Reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject
# to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#-----------------------------------------------------------------------------
# ExternalData
set(ExternalData_BINARY_ROOT   "${CMAKE_BINARY_DIR}/ExternalData")
set(ExternalData_OBJECT_STORES "${ExternalData_OBJECT_STORES}" "${CMAKE_CURRENT_SOURCE_DIR}/Testing/Data")
set(ExternalData_URL_TEMPLATES "http://simvascular.stanford.edu/downloads/public/vtkSV/Testing/Data/%(algo)/%(hash)")
set(ExternalData_LINK_CONTENT  MD5)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# vtkSV Test
set(VTKSV_TEST_DIR                 "${CMAKE_BINARY_DIR}/Testing/Temporary")
set(VTKSV_TEST_DATA_DIR            "${CMAKE_CURRENT_SOURCE_DIR}/Testing/Data")
set(VTKSV_TEST_BASELINE_DIR        "${VTKSV_TEST_DATA_DIR}/Baseline")
set(VTKSV_TEST_OUTPUT_DATA_DIR     "${ExternalData_BINARY_ROOT}/Testing/Data")
set(VTKSV_TEST_OUTPUT_BASELINE_DIR "${VTKSV_TEST_OUTPUT_DATA_DIR}/Baseline")
set(VTKSV_TEST_OUTPUT_DIR          "${CMAKE_BINARY_DIR}/Testing/Temporary")
set(VTKSV_DATA_ROOT                "${VTKSV_TEST_OUTPUT_DATA_DIR}")
make_directory(${VTKSV_TEST_DIR})
enable_testing()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Test coverage (setup is for macos 10.10 with xcode
if(VTKSV_TEST_COVERAGE)
  set(COVERAGE_COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/llvm-cov" CACHE STRING "Forcing llvm-cov" FORCE)
  set(CTEST_COVERAGE_COMMAND "${COVERAGE_COMMAND}" CACHE STRING "Forcing llvm-cov" FORCE)
  set(COVERAGE_EXTRA_FLAGS "gcov -l" CACHE STRING "Force coverage flags" FORCE)
  set(CTEST_COVERAGE_EXTRA_FLAGS "${COVERAGE_EXTRA_FLAGS}" CACHE STRING "Force coverage flags" FORCE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -Wall -W -Wshadow -Wunused-variable \
    -Wunused-parameter -Wunused-function -Wunused -Wno-system-headers \
    -Wno-deprecated -Woverloaded-virtual -Wwrite-strings -fprofile-arcs -ftest-coverage")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0 -Wall -W -fprofile-arcs -ftest-coverage")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")
endif()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# vtkSV Test
include(CTest)
#-----------------------------------------------------------------------------
