Introduction {#introduction}
==============

@tableofcontents

vtkSV is a set of VTK-based objects and filters that aid in geometric representation and manipulation.
It began as a series of filters that were developed for the software project SimVascular.

👉 SimVascular project page: [http://www.simvascular.org](http://www.simvascular.org)

There are a couple of toplevel directories with general macros and code, and then there is a directory containing various modules with the bulk of the code.
The following sections describe a bit about each of the directories and the code contained.
To see some examples, take a look at the documentation about examples (REF).

@section introduction-common Common

This is where all the general pound defines and general useful vtk functions are defined.
Other objects that are general to many of the modules should be added here.
Any function that is general and used in a variety of code should also become a function in *vtkSVGeneralUtils*.

@section introduction-thirdparty ThirdParty

This is where any potential third party code is placed.
Currently, that is just code from the Vascular Modeling Tool Kit (VMTK)

👉 VMTK project page:  [https://www.vmtk.org](https://www.vmtk.org)

@section introduction-testing Testing

@section introduction-travis Travis

Continuous integration is handled by travis CI.

@section introduction-modules Modules

The core of the code is contained with modules under this toplevel directory.

@subsection introduction-modules-boolean Boolean

@subsection introduction-modules-geometry Geometry

@subsection introduction-modules-nurbs NURBS

@subsection introduction-modules-parameterization Parameterization

@subsection introduction-modules-segmentation Segmentation

