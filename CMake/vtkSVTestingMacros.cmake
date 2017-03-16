
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

# -----------------------------------------------------------------------------
# _vtksv_test_parse_args(options source_ext args...)
#   INTERNAL: Parse arguments for testing functions.
#
#   Parses 'options' from the argument list into the 'options' variable in the
#   parent, Test instances found with the extension 'source_ext' are parsed
#   into the 'names' variable in the parent. Any comma-separated options after
#   the test instance is put into a '_${name}_options' variable for the test.
#   Any unrecognized arguments are put into the 'args' variable in the parent.
function(_vtksv_test_parse_args options source_ext)
  set(global_options)
  set(names)
  set(args)

  foreach(arg IN LISTS ARGN)
    set(handled 0)
    foreach(option IN LISTS options)
      if(arg STREQUAL option)
        list(APPEND global_options ${option})
        set(handled 1)
        break()
      endif()
    endforeach()
    if(handled)
      # Do nothing.
    elseif(source_ext AND "x${arg}" MATCHES "^x([^.]*)\\.${source_ext},?(.*)$")
      set(name "${CMAKE_MATCH_1}")
      string(REPLACE "," ";" _${name}_options "${CMAKE_MATCH_2}")
      list(APPEND names ${name})
    else()
      list(APPEND args ${arg})
    endif()
  endforeach()

  foreach(name IN LISTS names)
    set(_${name}_options "${_${name}_options}"
      PARENT_SCOPE)
  endforeach()
  set(options "${global_options}"
    PARENT_SCOPE)
  set(names "${names}"
    PARENT_SCOPE)
  set(args "${args}"
    PARENT_SCOPE)
endfunction()

# -----------------------------------------------------------------------------
# _vtksv_test_set_options(options prefix args...)
#   INTERNAL: Set variables related to options.
#
#   Looks in the arguments for options to set. Valid options are listed in the
#   'options' input list and the variables of the same name are set in the
#   parent scope to '1' if set, and '0' if they are not found. If 'prefix' is
#   non-empty, it is used as a prefix for the variable names to set and the
#   no-prefix variable is used as the unset value (rather than '0').
function(_vtksv_test_set_options options prefix)
  foreach(option IN LISTS options)
    set(default 0)
    if(prefix)
      set(default ${${option}})
    endif()
    set(${prefix}${option} ${default}
      PARENT_SCOPE)
  endforeach()
  foreach(option IN LISTS ARGN)
    set(${prefix}${option} 1
      PARENT_SCOPE)
  endforeach()
endfunction()

# -----------------------------------------------------------------------------
# _vtksv_test_parse_name(name)
#   INTERNAL: Parse the name of the test and the test file.
#
#   The 'name' argument must be of the form
#
#   [CustomTestName,]Test
#
#   where the CustomTestName followed by comma is optional. For example,
#   if 'name' has the value
#
#   Test1,Test
#
#   this function sets the variable 'test_name' to 'Test1' and
#   'test_file' to 'Test' in the parent scope. Note that the test file
#   does not include the file extension. For tests specified without a
#   custom name, .e.g.,
#
#   Test
#
#   both variables 'test_name' and 'test_file' will be set to the variable
#   'name' defined in the caller.
function(_vtksv_test_parse_name name)
  set(test_name ${name} PARENT_SCOPE)
  set(test_file ${name} PARENT_SCOPE)

  if(name AND "x${name}" MATCHES "^x([^,]*),(.*)$")
    set(test_name "${CMAKE_MATCH_1}" PARENT_SCOPE)
    set(test_file "${CMAKE_MATCH_2}" PARENT_SCOPE)
  endif()
endfunction()

# vtksv_add_test_cxx(exename tests [NO_DATA] [NO_VALID] [NO_OUTPUT]
#                  [test1.cxx...] [args...])
#   Adds C++ tests.
#
#   Adds tests using the 'exename' (which must be a CMake target) and the name
#   of the tests into the variable named by 'tests' in the parent scope. If the
#   NO_DATA option is specified, the test will not receive a -D argument (input file),
#   NO_VALID will suppress the -V argument (path to a baseline image), and
#   NO_OUTPUT will suppress the -T argument (output directory). Test-specific
#   arguments may be set to _${name}_ARGS. By default, the test name will be the part
#   of the source file before the '.cxx'. A custom test name can be specified by
#   giving a name followed by a comma before the test file name, .e.g.,
#
#   CustomTestName,TestSource.cxx
#
#   The 'vtk_test_prefix' variable may be set to create separate tests from a
#   single test name (e.g., running with different arguments), but should be
#   used only when required.
function(vtksv_add_test_cxx exename _tests)
  set(cxx_options
    NO_DATA
    NO_VALID
    NO_OUTPUT
    CUSTOM_BASELINES
    )
  _vtksv_test_parse_args("${cxx_options}" "cxx" ${ARGN})
  _vtksv_test_set_options("${cxx_options}" "" ${options})

  set(_vtksv_fail_regex "(\n|^)ERROR: " "instance(s)? still around")

  if(VTKSV_BASELINE_DIR)
    if(vtk-module)
      set(prefix ${vtk-module})
    elseif(vtk-example)
      set(prefix ${vtk-example})
    endif()
    set(baseline_dir ${VTKSV_BASELINE_DIR})
  elseif(vtk-module)
    set(prefix ${vtk-module})
    set(baseline_dir ${${vtk-module}_SOURCE_DIR}/Testing/Data/Baseline)
  elseif(vtk-example)
    set(prefix ${vtk-example})
    set(baseline_dir ${CMAKE_CURRENT_SOURCE_DIR}/Baseline)
  else()
    message(FATAL_ERROR "Neither vtk-module nor vtk-example is set!")
  endif()

  set(data_dir "${VTKSV_TEST_DATA_DIR}")
  if(${vtk-module}_DATA_DIR)
    set(data_dir "${${vtk-module}_DATA_DIR}")
  endif()

  set(externaldata_target VTKData)
  if(VTKSV_TEST_DATA_TARGET)
    set(externaldata_target ${VTKSV_TEST_DATA_TARGET})
  endif()

  foreach(name IN LISTS names)
    _vtksv_test_set_options("${cxx_options}" "local_" ${_${name}_options})
    _vtksv_test_parse_name(${name})

    set(_D "")
    if(NOT local_NO_DATA)
      set(_D -D ${data_dir})
    endif()

    set(_T "")
    if(NOT local_NO_OUTPUT)
      set(_T -T ${VTKSV_TEST_OUTPUT_DIR})
    endif()

    set(_V "")
    if(NOT local_NO_VALID)
      if(local_CUSTOM_BASELINES)
        set(_V -V "${data_dir}/Baseline")
      else()
        set(_V -V "DATA{${baseline_dir}/${test_name}.png,:}")
      endif()
    endif()

    ExternalData_add_test(${externaldata_target}
      NAME    ${prefix}Cxx-${vtk_test_prefix}${test_name}
      COMMAND $<TARGET_FILE:${exename}>
              ${test_file}
              ${args}
              ${${prefix}_ARGS}
              ${${name}_ARGS}
              ${_D} ${_T} ${_V})
    set_tests_properties(${prefix}Cxx-${vtk_test_prefix}${test_name}
      PROPERTIES
        LABELS "${${prefix}_TEST_LABELS}"
        FAIL_REGULAR_EXPRESSION "${_vtksv_fail_regex}"
      )

    list(APPEND ${_tests} "${test_file}")
  endforeach()

  set(${_tests} ${${_tests}} PARENT_SCOPE)
endfunction()

# -----------------------------------------------------------------------------
# vtksv_test_cxx_executable(exename, tests [RENDERING_FACTORY] [extra.cxx...])
#   Build a C++ test executable.
#
#   Creates a test executable for running the tests listed in the 'tests'
#   variable. If RENDERING_FACTORY is set, the rendering test driver will be
#   used instead. Any other sources found will be built into the executable as
#   well. Unrecognized arguments are ignored.
function(vtksv_test_cxx_executable exename _tests)
  set(exe_options
    RENDERING_FACTORY
    )
  _vtksv_test_parse_args("${exe_options}" "" ${ARGN})
  _vtksv_test_set_options("${exe_options}" "" ${options})

  if(NOT ${_tests})
    # No tests -> no need for an executable.
    return()
  endif()

  set(test_driver vtkTestDriver.h)
  if(RENDERING_FACTORY)
    include(vtkTestingRenderingDriver)
    set(test_driver ${vtkTestingRendering_SOURCE_DIR}/vtkTestingObjectFactory.h)
  endif()

  set(extra_sources ${args})

  if(vtk-module)
    set(CMAKE_TESTDRIVER_BEFORE_TESTMAIN
      "    vtksys::SystemInformation::SetStackTraceOnError(1);\n ${CMAKE_TESTDRIVER_BEFORE_TESTMAIN}")
  endif()

  create_test_sourcelist(test_sources ${exename}.cxx ${${_tests}}
    EXTRA_INCLUDE ${test_driver})

  if(vtk-module)
    vtk_module_test_executable(${exename} ${test_sources} ${extra_sources})
  elseif(vtk-example)
    add_executable(${exename} ${test_sources} ${extra_sources})
    target_link_libraries(${exename} ${VTK_LIBRARIES})
  else()
    message(FATAL_ERROR "Neither vtk-module nor vtk-example is set!")
  endif()
endfunction()
