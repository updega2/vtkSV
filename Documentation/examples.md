Examples {#examples}
==============

@tableofcontents

A variety of executables exist for the vtkSV functionality.
The following demonstrates some of the possibilities with the code.
In order to get help for a specific executable, use the -h flag.
This will demonstrate the typical usage for the executable as well as the options, definitions, and default values.

```bash
./vtkSVExample -h
```

@section examples-boolean Boolean

@section examples-geometry Misc

@section examples-parameterization Parameterization

@section examples-segmentation Segmentation

@subsection examples-segmentations-centerlines_extractor Centerlines Extractor

The centerlines extractor takes a triangulated surface and computes the centerline structure.
There are many options to this code to modify the output centerlines.
The most simple options are demonstrated here.
TODO: Longer explanation of how centerlines extractor works.

```bash
vtkSVCenterlinesExtractor -input simple_bifurcation.vtp -output simple_bifurcation_full_centerlines.vtp
```

@image html simple_bifurcation_centerline.png

```bash
vtkSVCenterlinesExtractor -input simple_bifurcation.vtp -output simple_bifurcation_full_centerlines.vtp -seedselector pickpoints
```
