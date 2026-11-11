#include "immersedTurbulenceWallTreatment.H"
//
#ifndef ICHIBOUNDARYFOAM_ENABLE_TURBULENCE_WALL_TREATMENT
#define ICHIBOUNDARYFOAM_ENABLE_TURBULENCE_WALL_TREATMENT 1
#endif
//
#if ICHIBOUNDARYFOAM_ENABLE_TURBULENCE_WALL_TREATMENT
namespace Foam
{
namespace ichiBoundaryFoam
{
ImmersedTurbulenceWallTreatment::
ImmersedTurbulenceWallTreatment
(
    const fvMesh& mesh,
    const ImmersedGeometry& geometry,
    const DirectForcingIbm& ibm,
    const dictionary& properties
)
:
    mesh_(mesh),
    geometry_(geometry),
    ibm_(ibm),
    properties_(properties),
    enabled_(properties_.lookupOrDefault<Switch>("enabled", false)),
    model_(properties_.lookupOrDefault<word>
    (
        "model",
        "kOmegaSSTLowRe"
    )),
    kFieldName_(properties_.lookupOrDefault<word>("kField", "k")),
    omegaFieldName_
    (
        properties_.lookupOrDefault<word>("omegaField", "omega")
    ),
    nutFieldName_
    (
        properties_.lookupOrDefault<word>("nutField", "nut")
    ),
    bandWidthFactor_
    (
        properties_.lookupOrDefault<scalar>("bandWidthFactor", 2.0)
    ),
    minimumWallDistanceFactor_
    (
        properties_.lookupOrDefault<scalar>
        (
            "minimumWallDistanceFactor",
            0.05
        )
    ),
    kWallValue_
    (
        properties_.lookupOrDefault<scalar>("kWallValue", 1e-12)
    ),
    beta1_(properties_.lookupOrDefault<scalar>("beta1", 0.075)),
    omegaViscousCoefficient_
    (
        properties_.lookupOrDefault<scalar>
        (
            "omegaViscousCoefficient",
            6.0
        )
    ),
    liquidThreshold_
    (
        properties_.lookupOrDefault<scalar>("liquidThreshold", 0.5)
    ),
    logDiagnostics_
    (
        properties_.lookupOrDefault<Switch>("logDiagnostics", true)
    ),
    writeDiagnosticFields_
    (
        properties_.lookupOrDefault<Switch>
        (
            "writeDiagnosticFields",
            false
        )
    ),
    diagnosticInterval_
    (
        properties_.lookupOrDefault<label>("diagnosticInterval", 100)
    ),
    correctionCount_(0),
    wallMask_(nullptr),
    omegaWallTarget_(nullptr)
{
    validate();
    //
    if (enabled_)
    {
        createDiagnosticFields();
        Info<< "Immersed turbulence-wall treatment: enabled" << nl
            << "    model: " << model_ << nl
            << "    fields: " << kFieldName_ << ", "
            << omegaFieldName_ << ", " << nutFieldName_ << nl
            << "    SDF band width: " << bandWidthFactor_
            << " local cell lengths" << nl
            << "    diagnostic interval: " << diagnosticInterval_
            << " corrections" << nl
            << "    diagnostic fields: "
            << (writeDiagnosticFields_ ? "enabled" : "disabled")
            << nl << endl;
    }
    else
    {
        Info<< "Immersed turbulence-wall treatment: disabled at "
            << "run time" << nl << endl;
    }
}
//
void ImmersedTurbulenceWallTreatment::validate() const
{
    if (!enabled_)
    {
        return;
    }
    //
    if (!ibm_.enabled())
    {
        FatalErrorInFunction
            << "The immersed turbulence-wall treatment requires an "
            << "enabled IBM model." << exit(FatalError);
    }
    //
    if (ibm_.velocityBoundaryCondition() != "noSlip")
    {
        FatalErrorInFunction
            << "The kOmegaSSTLowRe immersed-wall treatment is only "
            << "valid with IBM velocityBoundaryCondition noSlip." << nl
            << "The configured condition is "
            << ibm_.velocityBoundaryCondition() << '.'
            << exit(FatalError);
    }
    //
    if (model_ != "kOmegaSSTLowRe")
    {
        FatalErrorInFunction
            << "Unsupported immersed turbulence-wall model: "
            << model_ << nl
            << "Supported model: kOmegaSSTLowRe"
            << exit(FatalError);
    }
    //
    if
    (
        bandWidthFactor_ <= 0
     || minimumWallDistanceFactor_ <= 0
     || kWallValue_ < 0
     || beta1_ <= 0
     || omegaViscousCoefficient_ <= 0
     || diagnosticInterval_ < 1
    )
    {
        FatalErrorInFunction
            << "Invalid immersed turbulence-wall parameters."
            << exit(FatalError);
    }
}
//
void ImmersedTurbulenceWallTreatment::createDiagnosticFields()
{
    if (!writeDiagnosticFields_)
    {
        return;
    }
    //
    wallMask_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "immersedTurbulenceWallMask",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
    omegaWallTarget_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "immersedOmegaWallTarget",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                dimless/dimTime,
                0.0
            )
        )
    );
}
//
void ImmersedTurbulenceWallTreatment::constrain
(
    const volScalarField& alpha,
    const volScalarField& nu
)
{
    if (!enabled_)
    {
        return;
    }
    //
    if
    (
        !mesh_.foundObject<volScalarField>(kFieldName_)
     || !mesh_.foundObject<volScalarField>(omegaFieldName_)
    )
    {
        FatalErrorInFunction
            << "The immersed turbulence-wall treatment cannot find "
            << "fields " << kFieldName_ << " and "
            << omegaFieldName_ << "." << nl
            << "Enable it only for a compatible k-omega model."
            << exit(FatalError);
    }
    //
    volScalarField& k =
        mesh_.lookupObjectRef<volScalarField>(kFieldName_);
    volScalarField& omega =
        mesh_.lookupObjectRef<volScalarField>(omegaFieldName_);
    volScalarField* nut = nullptr;
    if (mesh_.foundObject<volScalarField>(nutFieldName_))
    {
        nut = &mesh_.lookupObjectRef<volScalarField>(nutFieldName_);
    }
    scalarField& kValues = k.primitiveFieldRef();
    scalarField& omegaValues = omega.primitiveFieldRef();
    scalarField* nutValues = nullptr;
    if (nut)
    {
        nutValues = &nut->primitiveFieldRef();
    }
    const scalarField& alphaValues = alpha.primitiveField();
    const scalarField& nuValues = nu.primitiveField();
    const scalarField& distance =
        geometry_.signedDistance().primitiveField();
    //
    scalarField* maskValues = nullptr;
    scalarField* targetValues = nullptr;
    if (writeDiagnosticFields_)
    {
        maskValues = &wallMask_().primitiveFieldRef();
        targetValues = &omegaWallTarget_().primitiveFieldRef();
        *maskValues = scalar(0);
        *targetValues = scalar(0);
    }
    //
    label nTreated = 0;
    label nLiquid = 0;
    scalar minimumDistance = GREAT;
    scalar maximumDistance = 0;
    scalar maximumKChange = 0;
    scalar maximumOmegaChange = 0;
    scalar maximumOmegaTarget = 0;
    scalar maximumNutChange = 0;
    //
    forAll(distance, celli)
    {
        const scalar d = distance[celli];
        if (d <= 0)
        {
            continue;
        }
        //
        const scalar cellLength = Foam::pow
        (
            mesh_.V()[celli],
            1.0/3.0
        );
        const scalar bandWidth = bandWidthFactor_*cellLength;
        if (d >= bandWidth)
        {
            continue;
        }
        //
        const scalar weight = 1.0 - d/bandWidth;
        const scalar effectiveDistance = Foam::max
        (
            d,
            minimumWallDistanceFactor_*cellLength
        );
        const scalar omegaTarget =
            omegaViscousCoefficient_*nuValues[celli]
           /(beta1_*sqr(effectiveDistance));
        const scalar oldK = kValues[celli];
        const scalar oldOmega = omegaValues[celli];
        const scalar newK = Foam::max
        (
            kWallValue_,
            (1.0 - weight)*oldK + weight*kWallValue_
        );
        const scalar blendedOmega =
            (1.0 - weight)*oldOmega + weight*omegaTarget;
        const scalar newOmega = Foam::max(oldOmega, blendedOmega);
        scalar oldNut = 0;
        scalar newNut = 0;
        if (nutValues)
        {
            oldNut = (*nutValues)[celli];
            newNut = (1.0 - weight)*oldNut;
        }
        //
        kValues[celli] = newK;
        omegaValues[celli] = newOmega;
        if (nutValues)
        {
            (*nutValues)[celli] = newNut;
        }
        if (writeDiagnosticFields_)
        {
            (*maskValues)[celli] = weight;
            (*targetValues)[celli] = omegaTarget;
        }
        //
        ++nTreated;
        if (alphaValues[celli] >= liquidThreshold_)
        {
            ++nLiquid;
        }
        minimumDistance = Foam::min(minimumDistance, d);
        maximumDistance = Foam::max(maximumDistance, d);
        maximumKChange = Foam::max
        (
            maximumKChange,
            Foam::mag(newK - oldK)
        );
        maximumOmegaChange = Foam::max
        (
            maximumOmegaChange,
            Foam::mag(newOmega - oldOmega)
        );
        maximumOmegaTarget = Foam::max
        (
            maximumOmegaTarget,
            omegaTarget
        );
        maximumNutChange = Foam::max
        (
            maximumNutChange,
            Foam::mag(newNut - oldNut)
        );
    }
    //
    // Only internal cells are constrained above.  Do not re-evaluate
    // physical-patch wall functions here: model-specific temporary
    // objects, such as kOmegaSST:G, no longer exist after
    // turbulence->correct() returns.
    if (writeDiagnosticFields_)
    {
        wallMask_().correctBoundaryConditions();
        omegaWallTarget_().correctBoundaryConditions();
    }
    //
    ++correctionCount_;
    const bool reportNow =
        logDiagnostics_
     && (correctionCount_ == 1
      || correctionCount_ % diagnosticInterval_ == 0);
    if (!reportNow)
    {
        return;
    }
    //
    reduce(nTreated, sumOp<label>());
    reduce(nLiquid, sumOp<label>());
    reduce(minimumDistance, minOp<scalar>());
    reduce(maximumDistance, maxOp<scalar>());
    reduce(maximumKChange, maxOp<scalar>());
    reduce(maximumOmegaChange, maxOp<scalar>());
    reduce(maximumOmegaTarget, maxOp<scalar>());
    reduce(maximumNutChange, maxOp<scalar>());
    //
    if (nTreated == 0)
    {
        minimumDistance = 0;
    }
    Info<< "Immersed turbulence-wall diagnostics: cells = "
        << nTreated << " (liquid = " << nLiquid
        << ", gas/interface = " << nTreated - nLiquid << ")"
        << ", distance = [" << minimumDistance << ", "
        << maximumDistance << "]"
        << ", max|delta k| = " << maximumKChange
        << ", max|delta omega| = " << maximumOmegaChange
        << ", max omegaWall = " << maximumOmegaTarget
        << ", max|delta nut| = " << maximumNutChange
        << nl << endl;
}
//
bool ImmersedTurbulenceWallTreatment::enabled() const
{
    return enabled_;
}
}
}
#endif
