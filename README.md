# ichiBoundaryFoam

`ichiBoundaryFoam` is an OpenFOAM-com v2312 `interFoam`-like VOF
solver with an analytic SDF geometry and an experimental
direct-forcing immersed-boundary backend.

The VOF and pressure-velocity algorithms follow the official v2312
solver.  Geometry, prescribed rigid-body motion, and IBM enforcement
are separate components.  IBM can be enabled for a case or left off
for an `interFoam` regression run.

Copy `examples/constant/ichiBoundaryGeometry` to the case `constant`
directory.  The first geometry implementation supports one finite
capped cylinder with either `static` or prescribed `harmonic6DoF`
motion.

Translation uses metres, frequency uses hertz, and rotation and phase
use radians.  Rotations are composed as `Rz*Ry*Rx`.

## Building

Load the OpenFOAM-com v2312 module and run:

```sh
make
```

The direct-forcing backend is compiled by default.  A fluid-only
build can exclude it with:

```sh
make ENABLE_DIRECT_FORCING=0
```

Immersed-surface slip diagnostics are excluded by default because
they traverse every cell during each residual report and add five
fields to every normal field write.  Enable them in a clean IBM
rebuild with:

```sh
make ibm ENABLE_SLIP_SAMPLING=1
```

Use `ENABLE_SLIP_SAMPLING=0` to compile out the fields, sampling loop,
parallel reductions, and associated log messages completely.  As
with the other preprocessor options, changing this option requires a
clean rebuild.

IBM constraint-impulse load accounting is compiled by default.  It
can be removed, together with its vector field and reductions, with:

```sh
make ibm ENABLE_LOAD_CALCULATION=0
```

For every final pressure-correction cycle, the solver accumulates the
actual IBM velocity increments as `rho*delta(U)` and reports the
equal-and-opposite body impulse, time-averaged force, angular impulse,
and moment about the instantaneous body centre.  The diagnostic does
not feed the load back to the prescribed body motion and therefore
does not alter the fluid equations.  It is an IBM constraint load; it
does not yet provide separate pressure and viscous surface integrals.

## Immersed turbulence-wall treatment

The SDF-based turbulence-wall treatment is compiled by default but is
disabled at run time unless the following case dictionary is present
and explicitly enabled:

```text
immersedTurbulenceWallTreatment
{
    enabled                     false;
    model                       kOmegaSSTLowRe;
    kField                      k;
    omegaField                  omega;
    nutField                    nut;
    bandWidthFactor             2;
    minimumWallDistanceFactor   0.05;
    kWallValue                  1e-12;
    beta1                       0.075;
    omegaViscousCoefficient     6;
    liquidThreshold             0.5;
    logDiagnostics              true;
    diagnosticInterval          100;
    writeDiagnosticFields       false;
}
```

The first implementation is a low-Reynolds-number treatment for a
`noSlip` immersed wall.  After the standard turbulence correction it
uses positive fluid-side SDF values and the local cubic cell length to
form a narrow wall band.  Within that band it smoothly damps `k` and
`nut` towards zero and imposes a viscous-scale lower bound on `omega`:

```text
omegaWall = omegaViscousCoefficient*nu/(beta1*yEff^2)
```

It deliberately rejects `freeSlip`, since applying this low-Re wall
law there would create an unphysical shear layer.  It also does not
replace the internal OpenFOAM `wallDist` used by SST blending
functions and is therefore an initial immersed-wall closure, not a
complete immersed wall-function implementation.

The lightweight diagnostics report treated-cell counts split by the
configured liquid threshold, SDF distance range, and maximum changes
to `k`, `omega`, and `nut`.  Parallel reductions occur only at
`diagnosticInterval`.  Optional fields
`immersedTurbulenceWallMask` and `immersedOmegaWallTarget` are written
only when `writeDiagnosticFields true`.

The feature can be compiled out completely with:

```sh
make ibm ENABLE_TURBULENCE_WALL_TREATMENT=0
```

Build the pure `interFoam` regression path with:

```sh
make pure-interfoam
```

This clean rebuild keeps the executable name `ichiBoundaryFoam`, but
does not construct the geometry, SDF, or IBM objects.  It directly uses
the OpenFOAM v2312 `UEqn.H` and `pEqn.H` implementations.  Rebuild the
regular development solver after the regression with:

```sh
make ibm
```

Both targets clean the old objects before compiling because `wmake`
does not track changes in preprocessor flags by itself.

The executable is installed under:

```text
platforms/<WM_OPTIONS>/bin/ichiBoundaryFoam
```

The top-level Makefile works around the compiler-wrapper ordering used
by the current Slurm module.  It delegates the OpenFOAM build itself to
`wmake`.

## One-way DPM scaffold

Passive Lagrangian particle tracking is compiled out by default.  It
can be included in an IBM build with:

```sh
make ibm ENABLE_ONE_WAY_DPM=1
```

The resulting executable constructs an OpenFOAM `basicKinematicCloud`.
The default cloud name is `kinematicCloud`; an alternative name can be
selected with the standard command-line form:

```sh
ichiBoundaryFoam -cloud sedimentCloud
```

The case must provide the matching standard cloud dictionary, for
example `constant/kinematicCloudProperties`.  Its `solution` dictionary
should set `coupled no`.  Carrier interpolation must include `rho`, `U`,
and `oneWayDpmMu`.  Particle injection, forces, wall interaction, and
cloud functions use the ordinary OpenFOAM kinematic-cloud models.
An editable inlet-injection example is provided at
`examples/constant/kinematicCloudProperties.oneWayDpm`.

This first implementation is deliberately one-way coupled:

* particles sample the instantaneous VOF mixture density, velocity,
  dynamic viscosity, and gravity;
* particles are advanced once after the final PIMPLE corrector of each
  fluid time step;
* particle momentum sources and occupied particle volume are never
  inserted into `UEqn`, `alpha.water`, or any other Eulerian equation;
* cloud data are written by the standard OpenFOAM cloud machinery.

The current scaffold has two important limitations.  Carrier density
and viscosity are the algebraic VOF mixture fields, so particles are
not yet restricted to water.  Also, standard mesh-patch interactions
do not see the SDF immersed surface, so particles can pass through an
IBM cylinder.  An SDF-aware particle interaction and phase-aware
tracking policy are separate follow-up developments.  This feature is
therefore suitable for passive-pathline experiments, not yet for
predicting scour or sediment feedback.

## VOF model

The fluid part follows the official OpenFOAM-com v2312 `interFoam`
algorithm and uses its maintained solver snippets.  It includes:

* algebraic VOF transport of `alpha.water` with MULES;
* interface compression;
* incompressible two-phase mixture properties;
* surface tension and gravity through `p_rgh`;
* PIMPLE pressure-velocity coupling;
* optional dynamic mesh and local time stepping support;
* the standard interfacial turbulence model selection.

An existing v2312 `interFoam` case therefore supplies the fluid-side
case structure.  The additional required input is
`constant/ichiBoundaryGeometry`.

## Validation convection scheme

Current cylinder and wave-validation cases use the bounded vector
convection scheme below as their project baseline:

```text
div(rhoPhi,U) Gauss limitedLinearV 1;
```

The earlier `Gauss linearUpwind grad(U)` setting was sensitive to a
loss of phase-fraction boundedness in cells next to the wave inlet.
Changing this scheme is treated as a separate numerical experiment and
must be reported explicitly.

## Geometry and prescribed motion

The sign convention is:

* negative inside the capped cylinder;
* zero on its surface;
* positive outside it.

The geometry evaluates the SDF, outward normal, and local rigid-body
velocity at cell and boundary-face centres.  It writes:

* `signedDistance`;
* `immersedSurfaceNormal`;
* `immersedSurfaceVelocity`.

The `harmonic6DoF` motion is prescribed independently for all three
translations and all three roll-pitch-yaw angles:

```text
q(t) = qRef + amplitude*sin(2*pi*frequency*(t-referenceTime)+phase)
```

This is currently prescribed one-way geometry motion.  It is not a
six-degree-of-freedom rigid-body dynamics solver.

## Direct-forcing IBM

The optional `ibm` dictionary is:

```text
ibm
{
    enabled                         false;
    couplingMethod                  cellBased;
    forcingMethod                   diffuseVelocity;
    transitionWidth                 0;
    sharpInterpolationLimit         1;
    markerSupportWidth              2;
    multiDirectForcing              false;
    nMultiDirectForcingCorrectors   1;
    slipTolerance                   0;
    slipSamplingBandWidth           0;
    slipLiquidThreshold             0.5;
    extendAlphaIntoSolid            false;
    nAlphaExtensionLayers           5;
    alphaExtensionWeightPower       2;
    boundExtendedAlpha              true;
}
```

IBM selection now has two levels. `couplingMethod` first chooses how
the Eulerian flow and immersed geometry are coupled. `forcingMethod`
then chooses the forcing implementation within that route.

* `cellBased` retains the existing cell-classification route and
  supports `diffuseVelocity` and `sharpReconstruction`.
* `markerBased` selects the marker interpolation/spreading route and
  currently supports `markerDirectForcing`.

If `couplingMethod` is omitted, it defaults to `cellBased` for backward
compatibility. Invalid coupling/forcing combinations stop during case
setup instead of silently selecting another implementation.

The forcing implementations are:

* `diffuseVelocity` preserves the original implementation.  It blends
  the predictor velocity with the rigid-body velocity using the smooth
  immersed-solid fraction.  `transitionWidth` controls the width of
  this diffuse forcing region.
* `sharpReconstruction` uses binary SDF classification.  Fluid cells
  remain unmodified, inner solid cells receive the rigid-body velocity,
  and the first fluid-side layer adjacent to the solid is treated as a
  forcing layer.  Its target velocity is currently obtained by bounded
  linear SDF interpolation between the immersed-surface velocity and
  point-connected fluid neighbours.
* `markerDirectForcing` is the first operator prototype for a
  Constant--Pinelli-like route.  It projects fluid-side forcing-cell
  centres onto the SDF zero surface, interpolates velocity to these
  markers with compact cosine weights, evaluates the marker residual,
  and spreads the correction back to the Eulerian cells.

The sharp backend is an experimental first reconstruction backend.  It
establishes the correct software and cell-classification path towards
the sharp DF-IBM of Ye et al., but it does not yet reproduce their
second-order Laplacian interpolation or immersed flux treatment.

The marker backend is also experimental and is disabled unless it is
selected explicitly.  Its stencils are currently processor-local and
the pressure equation is unchanged.  It is intended to validate the
marker interpolation/spreading operators before adding parallel halo
stencils, conservative adjoint spreading, and a coupled pressure--IBM
iteration.  Startup diagnostics report marker count, fallback stencil
count, and the maximum first-moment error.

The model writes:

* `immersedSolidFraction`, zero in fluid and one in solid;
* `immersedCellType`, with 0 fluid, 1 interface, and 2 solid;
* `immersedAcceleration`, the equivalent predictor-correction
  acceleration.

The IBM constrains `HbyA` to the prescribed rigid-body velocity before
`phiHbyA` is assembled.  The ordinary pressure equation is then solved
over the complete background mesh, including the fictitious solid
region.  Both `U` and `phi` are produced by this pressure projection;
neither field is overwritten after the projection.

When `multiDirectForcing` is enabled, the solver performs additional
forcing-projection correctors after the ordinary pressure equation.
Each corrector constrains the velocity with the same SDF weight, solves
a pressure-increment equation over the complete background domain, and
updates `U`, `phi`, and `p_rgh` consistently.  It does not create an
immersed pressure boundary, so the pressure field remains continuous.

`nMultiDirectForcingCorrectors` is the maximum number of additional
correctors.  A positive `slipTolerance` permits early termination when
the maximum weighted no-slip residual falls below the specified value.
The pressure-increment equation uses the standard `p_rghFinal` solver
entry from `system/fvSolution`.

The solver reports the maximum post-projection velocity slip in the
immersed region and the maximum discrete velocity divergence.  The
multi-direct-forcing loop reports the slip after every corrector.

When compiled with `ENABLE_SLIP_SAMPLING=1`, a positive
`slipSamplingBandWidth` enables a separate immersed-surface
diagnostic.  All cell centres satisfying
`abs(signedDistance) <= slipSamplingBandWidth` are sampled along the
complete finite cylinder, including its side and caps.  This is a
cell-centre SDF band, rather than exact interpolation onto the
zero-SDF surface.

At each residual report, the diagnostic logs:

* the total number of sampled cells;
* liquid and gas sample counts, separated by
  `slipLiquidThreshold`;
* maximum total slip for all, liquid, and gas samples;
* maximum normal and tangential slip magnitudes;
* the position, phase fraction, and signed distance of the largest
  sampled slip.

The relative velocity is computed as `U - immersedSurfaceVelocity`.
The following fields are retained in memory and written at ordinary
solver write times:

* `immersedSlipVelocity`;
* `immersedSlipMagnitude`;
* `immersedNormalSlip`;
* `immersedTangentialSlip`;
* `immersedSlipSamplingMask`.

The template wave-cylinder cases currently use a `0.0075 m` band and
an `alpha.water = 0.5` liquid threshold.  The band should be checked
against the local background-cell size: if no cell centre falls in
the band, the solver reports that condition instead of a maximum.

`transitionWidth 0` selects cell-centre sharp classification.  A
positive value applies a regularised transition over the specified
physical distance.

`extendAlphaIntoSolid` enables the optional phase-fraction extension
used by Ye et al. (2020).  After phase transport and before mixture
properties are corrected, values from cells on the fluid side of the
SDF are extended layer by layer into cells on the solid side.  Cells
sharing a mesh point form the interpolation stencil, with coupled-face
exchange at processor and cyclic boundaries.  The default
inverse-distance-square weight is selected by
`alphaExtensionWeightPower 2`.

This operation approximates a homogeneous normal-gradient condition
for `alpha.water` at the immersed surface.  Values produced inside the
solid are fictitious numerical extensions and must not be included in
the reported physical phase volume.  The option is disabled by default
to preserve the previous solver behaviour.  `boundExtendedAlpha`
clips only newly extended values to the physical interval `[0, 1]`.

## Current limitations

The first IBM backend constrains the momentum-predictor velocity, but
it does not yet:

* mask either fluid phase;
* receive hydrodynamic force feedback;
* separate its constraint load into pressure and viscous surface
  contributions;
* reconstruct newly uncovered fluid cells;
* provide a ghost-cell or cut-cell pressure treatment.

This implementation is suitable for algorithm development and the
next body-fitted comparison.  It is not yet a validated replacement
for the body-fitted monopile benchmark.

## Verification completed

The solver has been compiled with the cluster OpenFOAM-com v2312
module and run on the official laminar `damBreak` tutorial.  With IBM
disabled, at `t = 0.1 s`, against official `interFoam`:

* `U` was byte-identical and `max|delta U| = 0`;
* `max|delta p_rgh| = 0`;
* `max|delta alpha.water| = 2.95e-19`.

Two IBM checks used the same mesh with forcing enabled:

* a static cylinder produced zero velocity in all 82 full-solid
  cells;
* a translating cylinder produced `Ux = 0.0601363 m/s` in all 83
  full-solid cells, with zero slip at the written precision.

### Wave-model comparison

The official OpenFOAM-com v2312 `StokesI`, `StokesII`, and `StokesV`
wave tutorials were also run for one complete wave period with both
official `interFoam` and `ichiBoundaryFoam`.  Each run selected and
updated the wave model requested by `constant/waveProperties`.

The maximum absolute differences after one complete period were:

| Model | `alpha.water` | `mag(U)` [m/s] | `p_rgh` [Pa] |
| --- | ---: | ---: | ---: |
| StokesI | `1.05e-3` | `3.20e-4` | `5.38` |
| StokesII | `9.26e-4` | `3.95e-5` | `3.43` |
| StokesV | `1.00e-9` | `1.01e-8` | `1.00e-5` |

For Stokes II, the first visible difference was only about `1e-25`
in MULES clipping near zero.  The nonlinear time integration then
amplified this different floating-point trajectory over a full period.
The wave equations, boundary model, time stepping, and resolved wave
solution remained equivalent; bitwise identity is not expected once
the additional SDF calculations are present.

These are solver-to-solver regression tests with IBM disabled, not
validation against an analytical or experimental wave profile.  The
IBM-enabled wave-monopile solution still requires comparison with the
existing body-fitted case.

### Multi-direct-forcing smoke tests

An optional continuous-pressure multi-direct-forcing projection was
tested against the original single-forcing path.  The most important
results at the final simulated time were:

| Case | Single forcing slip | MDF slip | Single div(U) | MDF div(U) |
| --- | ---: | ---: | ---: | ---: |
| Static, full water | `1.30e-4` | `3.68e-5` | `2.83e-11` | `5.27e-12` |
| Translating, full water | `1.35e-4` | `8.48e-5` | `1.57e-11` | `1.16e-11` |
| Static, Stokes II | `1.77e-2` | `2.91e-3` | `1.30` | `1.42e-9` |

Four additional correctors were used for the static full-water test;
one was used for the translating and Stokes II tests.  In the Stokes II
smoke test, the observed maximum `alpha.water` decreased from `1.00065`
to `1.0`.  These results demonstrate algorithmic improvement, but they
are not a grid-converged wave-cylinder validation.
