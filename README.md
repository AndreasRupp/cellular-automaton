# Welcome to cellular-automaton

by **Andreas Rupp**, **Simon Zech**, and **Joona Lappalainen**.

It contains the C++ based library CAM implementing a simple [cellular automaton method (CAM)](
https://en.wikipedia.org/wiki/Cellular_automaton). In this CAM, single cells/pixels are allowed to
move within a [von Neumann neighborhood (VNN)](
https://en.wikipedia.org/wiki/Von_Neumann_neighborhood) of range `jump_parameter` and larger
particles, which are edge-connected sets of pixels, are allowed to move in a VNN depending on their
size, i.e., the number of the edge-connected pixels. The range of the VNN of larger particles can be
calculated as `jump_parameter` over the square-root of the larger particle's size. All particles,
i.e., single pixels and larger particles, move in such a way that the amount of their particle
neighbors is maximized.


# How to use cellular-automaton / CAM

There is two ways of using CAM:

- In `MATLAB`, one has to add the path of the `mex/CAM/` directory using MATLAB's `addpath`
  function. Then, one can run the CAM using `run_cam.m` as illustrated in `matlab_example.m`.
  Notably, for the MATLAB code to work smoothly in Ubuntu, you need to have `g++-10` installed,
  since it is currently the only compiler supporting both C++20 and MATLAB's supported libstdc++.
  In any case, you will also need `cmake`.
- In `C++`, one has to include the file `include/CAM/cellular_automaton.hxx`. Then, you can run the
  CAM as shown in `cpp_example.cxx` provided that you compile it using `-std=gnu++20`. A possible
  compilation command for the `cpp_example.cxx` is `clang++-12 -std=gnu++20 -Wall -Wextra -pedantic 
  -Iinclude -O3 examples_cpp/cpp_example.cxx -o test`. Here, `clang++-12` can be replaced by any
  suitable compiler implementing C++20.


# Copyright, License, and Contribution Policy

This directory contains the CAM library.

The CAM library is copyrighted by the authors of `cellular-automaton`. This term refers to the
people listed at the very top of this page.

The CAM library is free software; you can use it, redistribute it, and/or modify it under the terms
of the <b>GNU Lesser General Public License</b> as published by the Free Software Foundation; either
<b>version 2.1</b> of the License, or (at your option) any later version. The full text of the GNU
Lesser General Public version 2.1 is quoted in [License.txt](License.txt).


## Contributions

As a contributor to this project, you agree that all of your contributions be governed by the
<b>Developer Certificate of Origin version 1.1</b>. The CAM project does not require copyright
assignments for contributions. This means that the copyright for code contributions in the CAM
project is held by its respective contributors who have each agreed to release their contributed
code under a compatible open source license (LGPL v2.1 for library code). The full text of the 
Developer Certificate of Origin version 1.1 is quoted in [DeveloperCertificateOfOrigin.txt](
DeveloperCertificateOfOrigin.txt).


# Third party software

- The documentation for the MATLAB code is created using the submodule [doxymatlab](
  https://github.com/simgunz/doxymatlab). Please refer to  Simone Gaiarin's [GitHub page](
  https://github.com/simgunz/doxymatlab) and the [official file exchange page](
  https://se.mathworks.com/matlabcentral/fileexchange/25925-using-doxygen-with-matlab/)
  by Fabrice for additional information (e.g. considering licenses).


## Referencing the library

In addition to the terms imposed by the LGPL v2.1 or later, we ask for the following courtesy:

> Every publication presenting numerical results obtained with the help of CAM should state the name
> of the library and cite one or more of the following references  
> - N. Ray, A. Rupp, and A. Prechtel  
>   ***Discrete-continuum multiscale model for transport, biomass development and solid
    restructuring in porous media***  
>   Advances in Water Resources, doi: [10.1016/j.advwatres.2017.04.001](
    https://doi.org/10.1016/j.advwatres.2017.04.001)

This is the usual, fair way of giving credit to contributors to a scientific result. In addition, it
helps us justify our effort in developing CAM as an academic undertaking.


## Contact

For further questions regarding licensing and commercial use please contact Andreas Rupp directly
using [Email](mailto:info@rupp.ink).


## Links

- The license can be found in [License.txt](License.txt). It contains the [GNU Lesser General Public
License version 2.1](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html).
- The developer certificate of origin can be found in 
[DeveloperCertificateOfOrigin.txt](DeveloperCertificateOfOrigin.txt). It contains the [Developer 
Certificate of Origin version 1.1](https://developercertificate.org/).
