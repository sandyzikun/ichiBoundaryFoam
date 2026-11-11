# Cylinder validation matrix

This directory contains small OpenFOAM-com v2312 cases used to compare:

* `interFoam` on a body-fitted cylinder mesh;
* `ichiBoundaryFoam` on the same body-fitted mesh with IBM disabled;
* `ichiBoundaryFoam` on a Cartesian mesh with a static IBM cylinder;
* `ichiBoundaryFoam` with constant cylinder rotation;
* `ichiBoundaryFoam` with constant cylinder translation;
* `ichiBoundaryFoam` with a static IBM cylinder and Stokes II waves.

The first matrix uses uniform water and zero gravity.  It isolates the
velocity-pressure IBM coupling before adding a free surface or waves.

The cases exercise prescribed one-way WSI.  They do not yet solve the
cylinder dynamics from the computed fluid load, so they are not two-way
FSI cases.  Quantitative results and limitations are recorded in
`RESULTS.md`.

