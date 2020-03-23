# mcd-flip

The flip is part of our toolkit to tackle the 2020 Computational Geometry:
Solving Hard Optimization Problems (CG:SHOP) challenge of finding a minimal
convex decomposition of a pointset.

# Obtaining the source code

Clone the git repository:
The logging class is included via a git submodule, so clone the source using

    git clone --recurse-submodules https://github.com/cgalab/mcd-flip

# Building

To build the tool, run cmake and make:

    mkdir build &&
    cd build &&
    cmake -DCMAKE_BUILD_TYPE=Release .. &&
    make

# Running the tool

`mcd-flip` accepts a number of optional parameteres, see the `--help` output
for a brief overview.

It needs a pointset and an existing decomposition as input, given as a n `.obj`
file which has both vertex coordinate and line segments.  You can obtain one
from a pointset using the sister tool `mcd-recurse` with the `--full-obj`
option.

# License

`mcd-flip` is free software.  You may redistribute it and/or modify
it under the terms of the GNU General Public License (v3).
