#ifndef ICHIBOUNDARYFOAM_B361474C_D507_11F1_9472_A3B0B9D8CE68_DIRECTFORCINGIBM_H
#define ICHIBOUNDARYFOAM_B361474C_D507_11F1_9472_A3B0B9D8CE68_DIRECTFORCINGIBM_H
//
#include "fvCFD.H"
#include "immersedGeometry.H"

#ifndef ICHIBOUNDARYFOAM_ENABLE_SLIP_SAMPLING
#define ICHIBOUNDARYFOAM_ENABLE_SLIP_SAMPLING 0
#endif
#ifndef ICHIBOUNDARYFOAM_ENABLE_LOAD_CALCULATION
#define ICHIBOUNDARYFOAM_ENABLE_LOAD_CALCULATION 1
#endif
//
namespace Foam
{
namespace ichiBoundaryFoam
{
class DirectForcingIbm
{
    const fvMesh& mesh_;
    const ImmersedGeometry& geometry_;
    dictionary properties_;
    bool enabled_;
    word couplingMethod_;
    word forcingMethod_;
    word velocityBoundaryCondition_;
    scalar transitionWidth_;
    scalar sharpInterpolationLimit_;
    scalar markerSupportWidth_;
    scalar markerConstraintTolerance_;
    label markerConstraintMaxIterations_;
    bool requireMarkerConstraintConvergence_;
    bool multiDirectForcing_;
    label nMultiDirectForcingCorrectors_;
    scalar slipTolerance_;
    bool extendAlphaIntoSolid_;
    label nAlphaExtensionLayers_;
    scalar alphaExtensionWeightPower_;
    bool boundExtendedAlpha_;
    labelListList alphaExtensionStencil_;
    labelList markerSeedCells_;
    pointField markerPoints_;
    vectorField markerNormals_;
    vectorField markerVelocities_;
    scalarField markerAreaWeights_;
    labelListList markerStencils_;
    List<scalarList> markerWeights_;
    volScalarField solidFraction_;
    scalarField previousSolidFraction_;
    volScalarField cellType_;
#if ICHIBOUNDARYFOAM_ENABLE_LOAD_CALCULATION
    volVectorField constraintImpulseDensity_;
    label loadCorrectionCount_;
#endif
#if ICHIBOUNDARYFOAM_ENABLE_SLIP_SAMPLING
    scalar slipSamplingBandWidth_;
    scalar slipLiquidThreshold_;
    mutable volVectorField slipVelocity_;
    mutable volScalarField slipMagnitude_;
    mutable volScalarField normalSlip_;
    mutable volScalarField tangentialSlip_;
    mutable volScalarField slipSamplingMask_;
#endif
    //
    // [INFO] Selection of the constrained relative-velocity part
    vector constrainedRelativeVelocity
    (
        const vector& relativeVelocity,
        const label celli
    ) const;
    //
    // [INFO] Selection of a marker-relative velocity constraint
    vector constrainedMarkerRelativeVelocity
    (
        const vector& relativeVelocity,
        const vector& normal
    ) const;
    //
    // [INFO] Rebuilding of marker points and compact stencils
    void rebuildMarkerStencils();
    //
    // [INFO] Application of marker interpolation and spreading
    void constrainMarkers
    (
        volVectorField& predictor,
        volVectorField& correction
    ) const;
    //
    // [INFO] Interpolation from Eulerian cells to local markers
    vectorField interpolateMarkers
    (
        const vectorField& eulerianValues
    ) const;
    //
    // [INFO] Algebraic marker-to-cell spreading.  The normalized
    // form is the operator currently used by the production backend;
    // the unnormalized form is the transpose of interpolation.
    void spreadMarkers
    (
        const vectorField& markerValues,
        vectorField& eulerianValues,
        const bool normalizeOverlap
    ) const;

    // [INFO] Matrix-free application of J S in marker space
    vectorField applyMarkerSchur
    (
        const vectorField& markerValues
    ) const;
    //
    // [INFO] Application of diffuse velocity forcing
    void constrainDiffuse
    (
        volVectorField& predictor,
        volVectorField& correction
    ) const;
    //
    // [INFO] Application of sharp velocity reconstruction
    void constrainSharp
    (
        volVectorField& predictor,
        volVectorField& correction
    ) const;
    //
public:
    // [INFO] Construction of the direct-forcing IBM model
    DirectForcingIbm
    (
        const fvMesh& mesh,
        const ImmersedGeometry& geometry,
        const dictionary& properties
    );
    //
    // [INFO] Resetting of accumulated forcing for a new time step
    void beginTimeStep();
    //
    // [INFO] Resetting of load accounting for a pressure cycle
    void beginLoadCycle();
    //
    // [INFO] Update of fluid, forcing, and solid-cell classification
    void classify();
    //
    // [INFO] Extension of phase fraction into the immersed solid
    void extendAlpha(volScalarField& alpha) const;
    //
    // [INFO] Constraint of the momentum-predictor velocity
    void constrainPredictor
    (
        volVectorField& predictor,
        const volScalarField& density,
        const scalar deltaT
    );
    //
    // [INFO] Evaluation of the maximum velocity-constraint residual
    scalar maxSlip(const volVectorField& velocity) const;
    //
    // [INFO] Access to the multi-direct-forcing switch
    bool multiDirectForcing() const;
    //
    // [INFO] Access to the maximum number of forcing corrections
    label nMultiDirectForcingCorrectors() const;
    //
    // [INFO] Access to the velocity-constraint convergence tolerance
    scalar slipTolerance() const;
    //
    // [INFO] Reporting of the preliminary IBM force and moment
    void reportLoads(const scalar deltaT) const;
    //
    // [INFO] Reporting of phase volume in the accessible fluid region
    void reportPhaseVolume(const volScalarField& alpha) const;
    //
    // [INFO] Reporting of post-projection IBM residuals
    void reportResiduals
    (
        const volVectorField& velocity,
        const surfaceScalarField& flux,
        const volScalarField& alpha
    ) const;
    //
    // [INFO] Access to the IBM activation state
    bool enabled() const;
    //
    // [INFO] Access to the immersed coupling route
    const word& couplingMethod() const;
    //
    // [INFO] Access to the forcing implementation
    const word& forcingMethod() const;
    //
    // [INFO] Access to the immersed velocity boundary condition
    const word& velocityBoundaryCondition() const;
    //
    // [INFO] Local marker data exposed for operator verification
    const pointField& markerPoints() const;
    const vectorField& markerNormals() const;
    const vectorField& markerVelocities() const;
    const scalarField& markerAreaWeights() const;
    const labelListList& markerStencils() const;
    const List<scalarList>& markerWeights() const;
    //
    // [INFO] Production interpolation and spreading entry points for
    // the standalone marker-operator verification utility
    vectorField interpolateMarkerField
    (
        const vectorField& eulerianValues
    ) const;
    void spreadMarkerField
    (
        const vectorField& markerValues,
        vectorField& eulerianValues,
        const bool normalizeOverlap = true
    ) const;
    //
    // [INFO] Area/volume-weighted adjoint force spreading
    void spreadMarkerForceField
    (
        const vectorField& markerForces,
        vectorField& eulerianForceDensity
    ) const;
    //
    // [INFO] Access to the solid volume-fraction field
    const volScalarField& solidFraction() const;
};
}
}
//
#endif
