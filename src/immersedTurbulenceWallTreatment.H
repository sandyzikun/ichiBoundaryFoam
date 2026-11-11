#ifndef ICHIBOUNDARYFOAM_B361474C_D507_11F1_9472_A3B0B9D8CE68_IMMERSEDTURBULENCEWALLTREATMENT_H
#define ICHIBOUNDARYFOAM_B361474C_D507_11F1_9472_A3B0B9D8CE68_IMMERSEDTURBULENCEWALLTREATMENT_H
//
#include "fvCFD.H"
#include "immersedGeometry.H"
#include "directForcingIbm.H"
//
namespace Foam
{
namespace ichiBoundaryFoam
{
class ImmersedTurbulenceWallTreatment
{
    const fvMesh& mesh_;
    const ImmersedGeometry& geometry_;
    const DirectForcingIbm& ibm_;
    dictionary properties_;
    bool enabled_;
    word model_;
    word kFieldName_;
    word omegaFieldName_;
    word nutFieldName_;
    scalar bandWidthFactor_;
    scalar minimumWallDistanceFactor_;
    scalar kWallValue_;
    scalar beta1_;
    scalar omegaViscousCoefficient_;
    scalar liquidThreshold_;
    bool logDiagnostics_;
    bool writeDiagnosticFields_;
    label diagnosticInterval_;
    label correctionCount_;
    autoPtr<volScalarField> wallMask_;
    autoPtr<volScalarField> omegaWallTarget_;
    //
    // [INFO] Validation of the run-time model configuration
    void validate() const;
    //
    // [INFO] Allocation of optional immersed-wall diagnostic fields
    void createDiagnosticFields();
    //
public:
    // [INFO] Construction of the immersed turbulence-wall treatment
    ImmersedTurbulenceWallTreatment
    (
        const fvMesh& mesh,
        const ImmersedGeometry& geometry,
        const DirectForcingIbm& ibm,
        const dictionary& properties
    );
    //
    // [INFO] Constraint of internal near-wall turbulence fields
    void constrain
    (
        const volScalarField& alpha,
        const volScalarField& nu
    );
    //
    // [INFO] Access to the run-time activation state
    bool enabled() const;
};
}
}
//
#endif
