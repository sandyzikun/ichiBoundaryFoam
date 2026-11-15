# Marker operator test utility

`ichiMarkerOperatorTest` reads an OpenFOAM mesh and the normal
`constant/ichiBoundaryGeometry` dictionary, constructs the production
`DirectForcingIbm` marker operator, and exits without solving a flow
equation.

The tests cover marker geometry, constant, affine-linear, and quadratic
field interpolation, zero-field consistency, the unweighted algebraic
transpose relation, analytic surface area, force and torque conservation,
volume/area-weighted adjointness, and the matrix-free marker constraint
solve. Use `-strict` to apply
`-linearTolerance` as a pass/fail criterion to affine-linear reproduction.

After building the utility and generating the supplied mesh, it is
intended to be run as follows:

```sh
make marker-operator-test
blockMesh -case tests/cases/markerOperatorBox
platforms/$WM_OPTIONS/bin/ichiMarkerOperatorTest \
    -case tests/cases/markerOperatorBox -test all
```

Parallel runs require a decomposed copy of the test case. The production
operator uses OpenFOAM's processor-coupled point-connected cell stencil;
serial and decomposed runs should therefore give decomposition-independent
interpolation and spreading results.
