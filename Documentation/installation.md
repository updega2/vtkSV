Installation {#installation}
==============

@tableofcontents

@section installation-required_dependencies Required Dependencies

@subsection installation-required_dependencies-cmake CMake

CMake is used as the build and test system for vtkSV, following a very similar format to VTK.

👉 CMake binaries: [https://cmake.org/download/](https://cmake.org/download/)

It is recommended to install the CMake command line tools, but it is also possible to use the CMake GUI to configure, generate, and build.

@subsection installation-required_dependencies-vtk VTK

The other dependency of vtkSV is vtk (whoa, what a surprise).

👉 VTK source code: [https://github.com/Kitware/VTK](https://github.com/Kitware/VTK)

👉 VTK binaries: [https://vtk.org/download/](https://vtk.org/download/)

vtkSV uses CMake to configure and include dependencies, so the CMake configure file VTKConfig.cmake will need to be located.
If VTK was built from source, the VTKConfig.cmake will be under the build directory.

```bash
VTK-source/build-directory/VTKConfig.cmake
```

If VTK was installed as binaries, the VTKConfig.cmake will be underneath the library directory.

```bash
VTK-install-location/lib/cmake/vtk-{vtk_version}/VTKConfig.cmake
```

@section installation-building Building

Building with default options should be straight forward.
Create a build directory, use CMake to configure the project and generate the Makefiles.
Then use make to build the project.

```bash
mkdir Build
cd Build
ccmake ../ -DVTK_DIR=<Path/to/VTKConfig.cmake>
```
Press 'c' to configure and then 'g' to generate.

```bash
make
make install
```

@section installation-build_options Build Options

@subsection installation-build_options-testing Testing

Testing is off by default. To turn on:

```bash
BUILD_TESTING=ON
```

To run the unit tests, run

```bash
make test
```

or

```bash
ctest
```

See the ctest documentation to run specific tests and run with different options.

👉 ctest documentation: [https://cmake.org/cmake/help/v3.12/manual/ctest.1.html](https://cmake.org/cmake/help/v3.12/manual/ctest.1.html)

@subsection installation-build_options-executables Executables

For some parts of the codebase, executables have been created to easily run and use the code.

Executable scripts are off by default. To turn on:

```bash
VTKSV_BUILD_EXES=ON
```

The executables are placed underneath the bin directory.

@subsection installation-build_options-modules Modules

All modules are on by default. To turn off:

```bash
VTKSV_BUILD_MODULE_{MODULE_NAME}=OFF
```

<span style="color:red">NOTE</span> some modules depend on other modules, and you may in turn off other modules if you turn off one of its dependencies. It is easiest to just build all the libraries.




