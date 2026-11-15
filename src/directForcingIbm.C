#include "directForcingIbm.H"
#include "centredCPCCellToCellStencilObject.H"
#include "mapDistributeBase.H"
#include "syncTools.H"
//
namespace Foam
{
namespace ichiBoundaryFoam
{
namespace
{
template<class Type>
Field<Type> distributeCellField
(
    const mapDistribute& cellMap,
    const GeometricField<Type, fvPatchField, volMesh>& field
)
{
    Field<Type> distributedValues
    (
        cellMap.constructSize(),
        pTraits<Type>::zero
    );
    forAll(field, celli)
    {
        distributedValues[celli] = field[celli];
    }
    forAll(field.boundaryField(), patchi)
    {
        const fvPatchField<Type>& patchField =
            field.boundaryField()[patchi];
        label compacti =
            patchField.patch().start()
          - field.mesh().nInternalFaces()
          + field.mesh().nCells();
        forAll(patchField, facei)
        {
            distributedValues[compacti++] = patchField[facei];
        }
    }
    cellMap.distribute(distributedValues);
    return distributedValues;
}

bool affineCorrectWeights
(
    const point& marker,
    const vectorField& centres,
    const labelList& stencil,
    scalarList& weights
)
{
    scalar augmented[4][5] = {{0.0}};
    forAll(stencil, i)
    {
        const vector offset = centres[stencil[i]] - marker;
        const scalar basis[4] =
        {
            1.0,
            offset.x(),
            offset.y(),
            offset.z()
        };
        for (label row = 0; row < 4; ++row)
        {
            for (label col = 0; col < 4; ++col)
            {
                augmented[row][col] +=
                    weights[i]*basis[row]*basis[col];
            }
        }
    }
    augmented[0][4] = 1.0;

    scalar matrixScale = 0.0;
    for (label row = 0; row < 4; ++row)
    {
        for (label col = 0; col < 4; ++col)
        {
            matrixScale = max(matrixScale, mag(augmented[row][col]));
        }
    }
    const scalar pivotTolerance = max(matrixScale, scalar(1))*1e-12;
    for (label pivot = 0; pivot < 4; ++pivot)
    {
        label pivotRow = pivot;
        scalar pivotMagnitude = mag(augmented[pivot][pivot]);
        for (label row = pivot + 1; row < 4; ++row)
        {
            if (mag(augmented[row][pivot]) > pivotMagnitude)
            {
                pivotMagnitude = mag(augmented[row][pivot]);
                pivotRow = row;
            }
        }
        if (pivotMagnitude <= pivotTolerance)
        {
            return false;
        }
        if (pivotRow != pivot)
        {
            for (label col = pivot; col < 5; ++col)
            {
                const scalar temporary = augmented[pivot][col];
                augmented[pivot][col] = augmented[pivotRow][col];
                augmented[pivotRow][col] = temporary;
            }
        }
        const scalar diagonal = augmented[pivot][pivot];
        for (label col = pivot; col < 5; ++col)
        {
            augmented[pivot][col] /= diagonal;
        }
        for (label row = 0; row < 4; ++row)
        {
            if (row == pivot)
            {
                continue;
            }
            const scalar factor = augmented[row][pivot];
            for (label col = pivot; col < 5; ++col)
            {
                augmented[row][col] -=
                    factor*augmented[pivot][col];
            }
        }
    }

    const scalar coefficients[4] =
    {
        augmented[0][4],
        augmented[1][4],
        augmented[2][4],
        augmented[3][4]
    };
    forAll(stencil, i)
    {
        const vector offset = centres[stencil[i]] - marker;
        weights[i] *=
            coefficients[0]
          + coefficients[1]*offset.x()
          + coefficients[2]*offset.y()
          + coefficients[3]*offset.z();
    }
    return true;
}

scalar globalMarkerDot
(
    const vectorField& first,
    const vectorField& second
)
{
    scalar localDot = 0.0;
    forAll(first, i)
    {
        localDot += first[i] & second[i];
    }
    return returnReduce(localDot, sumOp<scalar>());
}
}

DirectForcingIbm::DirectForcingIbm
(
    const fvMesh& mesh,
    const ImmersedGeometry& geometry,
    const dictionary& properties
)
:
    mesh_(mesh),
    geometry_(geometry),
    properties_(properties),
    enabled_(properties_.lookupOrDefault<bool>("enabled", false)),
    couplingMethod_
    (
        properties_.lookupOrDefault<word>
        (
            "couplingMethod",
            "cellBased"
        )
    ),
    forcingMethod_
    (
        properties_.lookupOrDefault<word>
        (
            "forcingMethod",
            "diffuseVelocity"
        )
    ),
    velocityBoundaryCondition_
    (
        properties_.lookupOrDefault<word>
        (
            "velocityBoundaryCondition",
            "noSlip"
        )
    ),
    transitionWidth_
    (
        properties_.lookupOrDefault<scalar>("transitionWidth", 0.0)
    ),
    sharpInterpolationLimit_
    (
        properties_.lookupOrDefault<scalar>
        (
            "sharpInterpolationLimit",
            1.0
        )
    ),
    markerSupportWidth_
    (
        properties_.lookupOrDefault<scalar>
        (
            "markerSupportWidth",
            2.0
        )
    ),
    markerConstraintTolerance_
    (
        properties_.lookupOrDefault<scalar>
        (
            "markerConstraintTolerance",
            1.0e-8
        )
    ),
    markerConstraintMaxIterations_
    (
        properties_.lookupOrDefault<label>
        (
            "markerConstraintMaxIterations",
            300
        )
    ),
    requireMarkerConstraintConvergence_
    (
        properties_.lookupOrDefault<bool>
        (
            "requireMarkerConstraintConvergence",
            false
        )
    ),
    multiDirectForcing_
    (
        properties_.lookupOrDefault<bool>
        (
            "multiDirectForcing",
            false
        )
    ),
    nMultiDirectForcingCorrectors_
    (
        properties_.lookupOrDefault<label>
        (
            "nMultiDirectForcingCorrectors",
            1
        )
    ),
    slipTolerance_
    (
        properties_.lookupOrDefault<scalar>
        (
            "slipTolerance",
            0.0
        )
    ),
    extendAlphaIntoSolid_
    (
        properties_.lookupOrDefault<bool>
        (
            "extendAlphaIntoSolid",
            false
        )
    ),
    nAlphaExtensionLayers_
    (
        properties_.lookupOrDefault<label>
        (
            "nAlphaExtensionLayers",
            5
        )
    ),
    alphaExtensionWeightPower_
    (
        properties_.lookupOrDefault<scalar>
        (
            "alphaExtensionWeightPower",
            2.0
        )
    ),
    boundExtendedAlpha_
    (
        properties_.lookupOrDefault<bool>
        (
            "boundExtendedAlpha",
            true
        )
    ),
    alphaExtensionStencil_(mesh_.nCells()),
    markerSeedCells_(),
    markerPoints_(),
    markerNormals_(),
    markerVelocities_(),
    markerAreaWeights_(),
    markerStencils_(),
    markerWeights_(),
    solidFraction_
    (
        IOobject
        (
            "immersedSolidFraction",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, Zero)
    ),
    previousSolidFraction_(mesh_.nCells(), 0.0),
    cellType_
    (
        IOobject
        (
            "immersedCellType",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("fluid", dimless, Zero)
    )
#if ICHIBOUNDARYFOAM_ENABLE_LOAD_CALCULATION
    ,
    constraintImpulseDensity_
    (
        IOobject
        (
            "immersedConstraintImpulseDensity",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimDensity*dimVelocity,
            Zero
        )
    ),
    loadCorrectionCount_(0)
#endif
#if ICHIBOUNDARYFOAM_ENABLE_SLIP_SAMPLING
    ,
    slipSamplingBandWidth_
    (
        properties_.lookupOrDefault<scalar>
        (
            "slipSamplingBandWidth",
            0.0
        )
    ),
    slipLiquidThreshold_
    (
        properties_.lookupOrDefault<scalar>
        (
            "slipLiquidThreshold",
            0.5
        )
    ),
    slipVelocity_
    (
        IOobject
        (
            "immersedSlipVelocity",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimVelocity, Zero)
    ),
    slipMagnitude_
    (
        IOobject
        (
            "immersedSlipMagnitude",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimVelocity, Zero)
    ),
    normalSlip_
    (
        IOobject
        (
            "immersedNormalSlip",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimVelocity, Zero)
    ),
    tangentialSlip_
    (
        IOobject
        (
            "immersedTangentialSlip",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimVelocity, Zero)
    ),
    slipSamplingMask_
    (
        IOobject
        (
            "immersedSlipSamplingMask",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, Zero)
    )
#endif
{
#if !ICHIBOUNDARYFOAM_ENABLE_DIRECT_FORCING
    if (enabled_)
    {
        FatalIOErrorInFunction(properties_)
            << "The direct-forcing backend was disabled at build time"
            << exit(FatalIOError);
    }
#endif
    if (transitionWidth_ < 0.0)
    {
        FatalIOErrorInFunction(properties_)
            << "transitionWidth must not be negative"
            << exit(FatalIOError);
    }
    if
    (
        couplingMethod_ != "cellBased"
     && couplingMethod_ != "markerBased"
    )
    {
        FatalIOErrorInFunction(properties_)
            << "Unsupported couplingMethod: "
            << couplingMethod_ << nl
            << "Supported methods: cellBased, markerBased"
            << exit(FatalIOError);
    }
    if
    (
        couplingMethod_ == "cellBased"
     && forcingMethod_ != "diffuseVelocity"
     && forcingMethod_ != "sharpReconstruction"
    )
    {
        FatalIOErrorInFunction(properties_)
            << "forcingMethod " << forcingMethod_
            << " is not available for couplingMethod cellBased" << nl
            << "Supported methods: diffuseVelocity, "
            << "sharpReconstruction"
            << exit(FatalIOError);
    }
    if
    (
        couplingMethod_ == "markerBased"
     && forcingMethod_ != "markerDirectForcing"
    )
    {
        FatalIOErrorInFunction(properties_)
            << "forcingMethod " << forcingMethod_
            << " is not available for couplingMethod markerBased" << nl
            << "Supported method: markerDirectForcing"
            << exit(FatalIOError);
    }
    if
    (
        velocityBoundaryCondition_ != "noSlip"
     && velocityBoundaryCondition_ != "freeSlip"
    )
    {
        FatalIOErrorInFunction(properties_)
            << "Unsupported velocityBoundaryCondition: "
            << velocityBoundaryCondition_ << nl
            << "Supported conditions: noSlip, freeSlip"
            << exit(FatalIOError);
    }
    if (markerSupportWidth_ <= 0.0)
    {
        FatalIOErrorInFunction(properties_)
            << "markerSupportWidth must be positive"
            << exit(FatalIOError);
    }
    if
    (
        sharpInterpolationLimit_ <= 0.0
     || sharpInterpolationLimit_ > 1.0
    )
    {
        FatalIOErrorInFunction(properties_)
            << "sharpInterpolationLimit must be in (0, 1]"
            << exit(FatalIOError);
    }
    if (nMultiDirectForcingCorrectors_ < 1)
    {
        FatalIOErrorInFunction(properties_)
            << "nMultiDirectForcingCorrectors must be at least one"
            << exit(FatalIOError);
    }
    if (markerConstraintTolerance_ <= 0.0)
    {
        FatalIOErrorInFunction(properties_)
            << "markerConstraintTolerance must be positive"
            << exit(FatalIOError);
    }
    if (markerConstraintMaxIterations_ < 1)
    {
        FatalIOErrorInFunction(properties_)
            << "markerConstraintMaxIterations must be at least one"
            << exit(FatalIOError);
    }
    if (slipTolerance_ < 0.0)
    {
        FatalIOErrorInFunction(properties_)
            << "slipTolerance must not be negative"
            << exit(FatalIOError);
    }
#if ICHIBOUNDARYFOAM_ENABLE_SLIP_SAMPLING
    if (slipSamplingBandWidth_ < 0.0)
    {
        FatalIOErrorInFunction(properties_)
            << "slipSamplingBandWidth must not be negative"
            << exit(FatalIOError);
    }
    if
    (
        slipLiquidThreshold_ < 0.0
     || slipLiquidThreshold_ > 1.0
    )
    {
        FatalIOErrorInFunction(properties_)
            << "slipLiquidThreshold must be in [0, 1]"
            << exit(FatalIOError);
    }
#endif
    if (nAlphaExtensionLayers_ < 1)
    {
        FatalIOErrorInFunction(properties_)
            << "nAlphaExtensionLayers must be at least one"
            << exit(FatalIOError);
    }
    if (alphaExtensionWeightPower_ <= 0.0)
    {
        FatalIOErrorInFunction(properties_)
            << "alphaExtensionWeightPower must be positive"
            << exit(FatalIOError);
    }
    if (enabled_)
    {
        const labelListList& cellPoints = mesh_.cellPoints();
        const labelListList& pointCells = mesh_.pointCells();
        forAll(alphaExtensionStencil_, celli)
        {
            labelHashSet neighbours;
            for (const label pointi : cellPoints[celli])
            {
                for (const label neighbour : pointCells[pointi])
                {
                    if (neighbour != celli)
                    {
                        neighbours.insert(neighbour);
                    }
                }
            }
            alphaExtensionStencil_[celli] = neighbours.toc();
        }
    }
    Info<< "Direct-forcing IBM: "
        << (enabled_ ? "enabled" : "disabled") << nl;
    if (enabled_)
    {
        Info<< "IBM coupling method: " << couplingMethod_ << nl;
        Info<< "IBM forcing method: " << forcingMethod_ << nl;
        Info<< "IBM velocity boundary condition: "
            << velocityBoundaryCondition_ << nl;
        if (velocityBoundaryCondition_ == "freeSlip")
        {
            Info<< "Free-slip forcing applies the kinematic "
                << "no-penetration constraint only" << nl;
        }
        if (forcingMethod_ == "sharpReconstruction")
        {
            Info<< "Sharp reconstruction is experimental; "
                << "interpolation limit = "
                << sharpInterpolationLimit_ << nl;
        }
        if (couplingMethod_ == "markerBased")
        {
            Info<< "Marker direct forcing is experimental; "
                << "compact support width = "
                << markerSupportWidth_
                << " local cell widths" << nl;
            Info<< "Marker stencils use processor-coupled "
                << "point-connected cell halos" << nl;
            Info<< "Marker constraint solver: matrix-free CG, tolerance = "
                << markerConstraintTolerance_ << ", maximum iterations = "
                << markerConstraintMaxIterations_ << nl;
        }
        Info<< "Multi-direct forcing: "
            << (multiDirectForcing_ ? "enabled" : "disabled")
            << ", correctors = "
            << nMultiDirectForcingCorrectors_
            << ", slip tolerance = "
            << slipTolerance_
            << nl;
#if ICHIBOUNDARYFOAM_ENABLE_MDF_PRESSURE_PROJECTION
        Info<< "MDF pressure projection: enabled after every pass"
            << nl;
#else
        Info<< "MDF pressure projection: disabled; predictor-only "
            << "passes with the regular PIMPLE projection"
            << nl;
#endif
        Info<< "Alpha extension into the solid: "
            << (extendAlphaIntoSolid_ ? "enabled" : "disabled")
            << ", layers = "
            << nAlphaExtensionLayers_
            << ", inverse-distance power = "
            << alphaExtensionWeightPower_
            << nl;
#if ICHIBOUNDARYFOAM_ENABLE_SLIP_SAMPLING
        Info<< "IBM slip sampling: compiled, band width = "
            << slipSamplingBandWidth_
            << ", liquid threshold = "
            << slipLiquidThreshold_
            << nl;
#endif
#if ICHIBOUNDARYFOAM_ENABLE_LOAD_CALCULATION
        Info<< "IBM constraint-impulse load calculation: compiled"
            << nl;
#else
        Info<< "IBM load calculation: disabled at build time" << nl;
#endif
    }
}
//
void DirectForcingIbm::extendAlpha(volScalarField& alpha) const
{
    if (!enabled_ || !extendAlphaIntoSolid_)
    {
        return;
    }
    //
    scalarField values(alpha.primitiveField());
    scalarField known(mesh_.nCells(), 0.0);
    boolList modified(mesh_.nCells(), false);
    const scalarField& distanceToSurface =
        geometry_.signedDistance().primitiveField();
    const vectorField& centres = mesh_.cellCentres();
    //
    forAll(known, celli)
    {
        if (distanceToSurface[celli] > 0.0)
        {
            known[celli] = 1.0;
        }
    }
    //
    label totalExtended = 0;
    for (label layer = 0; layer < nAlphaExtensionLayers_; ++layer)
    {
        scalarField weightedSum(mesh_.nCells(), 0.0);
        scalarField weightSum(mesh_.nCells(), 0.0);
        //
        forAll(known, celli)
        {
            if
            (
                known[celli] > 0.5
             || distanceToSurface[celli] > 0.0
            )
            {
                continue;
            }
            //
            for
            (
                const label neighbour
              : alphaExtensionStencil_[celli]
            )
            {
                if (known[neighbour] < 0.5)
                {
                    continue;
                }
                //
                const scalar distance = max
                (
                    mag(centres[celli] - centres[neighbour]),
                    ROOTVSMALL
                );
                const scalar weight = pow
                (
                    distance,
                    -alphaExtensionWeightPower_
                );
                weightedSum[celli] += weight*values[neighbour];
                weightSum[celli] += weight;
            }
        }
        //
        List<scalar> neighbourKnown;
        List<scalar> neighbourValues;
        List<point> neighbourCentres;
        syncTools::swapBoundaryCellList
        (
            mesh_,
            known,
            neighbourKnown
        );
        syncTools::swapBoundaryCellList
        (
            mesh_,
            values,
            neighbourValues
        );
        syncTools::swapBoundaryCellPositions
        (
            mesh_,
            centres,
            neighbourCentres
        );
        //
        for (const polyPatch& patch : mesh_.boundaryMesh())
        {
            if (!patch.coupled())
            {
                continue;
            }
            //
            forAll(patch, patchFacei)
            {
                const label facei = patch.start() + patchFacei;
                const label bFacei = facei - mesh_.nInternalFaces();
                const label celli = patch.faceCells()[patchFacei];
                if
                (
                    known[celli] > 0.5
                 || distanceToSurface[celli] > 0.0
                 || neighbourKnown[bFacei] < 0.5
                )
                {
                    continue;
                }
                //
                const scalar distance = max
                (
                    mag(centres[celli] - neighbourCentres[bFacei]),
                    ROOTVSMALL
                );
                const scalar weight = pow
                (
                    distance,
                    -alphaExtensionWeightPower_
                );
                weightedSum[celli] +=
                    weight*neighbourValues[bFacei];
                weightSum[celli] += weight;
            }
        }
        //
        label layerExtended = 0;
        forAll(known, celli)
        {
            if (known[celli] < 0.5 && weightSum[celli] > VSMALL)
            {
                values[celli] = weightedSum[celli]/weightSum[celli];
                if (boundExtendedAlpha_)
                {
                    values[celli] = min(max(values[celli], 0.0), 1.0);
                }
                modified[celli] = true;
                ++layerExtended;
            }
        }
        //
        forAll(known, celli)
        {
            if (modified[celli])
            {
                known[celli] = 1.0;
            }
        }
        //
        const label globalLayerExtended = returnReduce
        (
            layerExtended,
            sumOp<label>()
        );
        totalExtended += globalLayerExtended;
        if (globalLayerExtended == 0)
        {
            break;
        }
    }
    //
    alpha.primitiveFieldRef() = values;
    alpha.correctBoundaryConditions();
    Info<< "Alpha extension: extended "
        << totalExtended
        << " solid cells over at most "
        << nAlphaExtensionLayers_
        << " layers"
        << nl;
}
//
void DirectForcingIbm::classify()
{
    const scalarField& distance =
        geometry_.signedDistance().primitiveField();
    scalarField& fraction = solidFraction_.primitiveFieldRef();
    scalarField& type = cellType_.primitiveFieldRef();
    previousSolidFraction_ = fraction;
    if (forcingMethod_ == "sharpReconstruction")
    {
        forAll(fraction, celli)
        {
            if (!enabled_ || distance[celli] > 0.0)
            {
                fraction[celli] = 0.0;
                type[celli] = 0.0;
            }
            else
            {
                fraction[celli] = 1.0;
                type[celli] = 2.0;
            }
        }
        forAll(type, celli)
        {
            if (type[celli] != 0.0)
            {
                continue;
            }
            for
            (
                const label neighbour
              : alphaExtensionStencil_[celli]
            )
            {
                if (distance[neighbour] <= 0.0)
                {
                    type[celli] = 1.0;
                    break;
                }
            }
        }
    }
    else
    {
        forAll(fraction, celli)
        {
            if (!enabled_)
            {
                fraction[celli] = 0.0;
                type[celli] = 0.0;
            }
            else if (transitionWidth_ <= VSMALL)
            {
                fraction[celli] =
                    distance[celli] <= 0.0 ? 1.0 : 0.0;
                type[celli] = fraction[celli] > 0.5 ? 2.0 : 0.0;
            }
            else if (distance[celli] <= -transitionWidth_)
            {
                fraction[celli] = 1.0;
                type[celli] = 2.0;
            }
            else if (distance[celli] >= transitionWidth_)
            {
                fraction[celli] = 0.0;
                type[celli] = 0.0;
            }
            else
            {
                const scalar argument =
                    constant::mathematical::pi
                   *distance[celli]/transitionWidth_;
                fraction[celli] =
                    0.5*(1.0 - distance[celli]/transitionWidth_)
                  - 0.5*sin(argument)
                   /constant::mathematical::pi;
                type[celli] = 1.0;
            }
        }
        if (enabled_ && transitionWidth_ <= VSMALL)
        {
            const labelUList& owners = mesh_.owner();
            const labelUList& neighbours = mesh_.neighbour();
            forAll(neighbours, facei)
            {
                const label owner = owners[facei];
                const label neighbour = neighbours[facei];
                if
                (
                    (distance[owner] <= 0.0)
                 != (distance[neighbour] <= 0.0)
                )
                {
                    type[owner] = 1.0;
                    type[neighbour] = 1.0;
                }
            }
        }
    }
    label forcingCells = 0;
    forAll(type, celli)
    {
        if (type[celli] == 1.0)
        {
            ++forcingCells;
        }
    }
    solidFraction_.correctBoundaryConditions();
    cellType_.correctBoundaryConditions();
    if (enabled_ && couplingMethod_ == "markerBased")
    {
        rebuildMarkerStencils();
    }
    scalar coveredVolume = 0.0;
    scalar uncoveredVolume = 0.0;
    forAll(fraction, celli)
    {
        const scalar change =
            fraction[celli] - previousSolidFraction_[celli];
        coveredVolume += max(change, scalar(0))*mesh_.V()[celli];
        uncoveredVolume +=
            max(-change, scalar(0))*mesh_.V()[celli];
    }
    Info<< "IBM classification: solid volume = "
        << gSum(fraction*mesh_.V())
        << ", forcing cells = "
        << returnReduce(forcingCells, sumOp<label>())
        << ", covered/uncovered volume = "
        << returnReduce(coveredVolume, sumOp<scalar>())
        << '/'
        << returnReduce(uncoveredVolume, sumOp<scalar>())
        << nl;
}
//
void DirectForcingIbm::rebuildMarkerStencils()
{
    DynamicList<label> seedCells;
    DynamicList<point> markerPoints;
    DynamicList<vector> markerNormals;
    DynamicList<vector> markerVelocities;
    const scalarField& distance =
        geometry_.signedDistance().primitiveField();
    const scalarField& type = cellType_.primitiveField();
    const vectorField& normal =
        geometry_.surfaceNormal().primitiveField();
    const vectorField& centres = mesh_.cellCentres();
    //
    forAll(type, celli)
    {
        if (type[celli] != 1.0 || distance[celli] < 0.0)
        {
            continue;
        }
        const point marker =
            centres[celli] - distance[celli]*normal[celli];
        scalar markerDistance = 0.0;
        vector markerNormal(vector::zero);
        vector markerVelocity(vector::zero);
        geometry_.evaluatePoint
        (
            marker,
            markerDistance,
            markerNormal,
            markerVelocity
        );
        seedCells.append(celli);
        markerPoints.append(marker);
        markerNormals.append(markerNormal);
        markerVelocities.append(markerVelocity);
    }
    seedCells.shrink();
    markerPoints.shrink();
    markerNormals.shrink();
    markerVelocities.shrink();
    markerSeedCells_.transfer(seedCells);
    markerPoints_.transfer(markerPoints);
    markerNormals_.transfer(markerNormals);
    markerVelocities_.transfer(markerVelocities);
    markerAreaWeights_.setSize(markerPoints_.size(), 0.0);
    label localSideMarkers = 0;
    label localPositiveCapMarkers = 0;
    label localNegativeCapMarkers = 0;
    const vector cylinderAxis = geometry_.axis();
    for (const vector& markerNormal : markerNormals_)
    {
        const scalar axialComponent = markerNormal & cylinderAxis;
        if (axialComponent > 0.5)
        {
            ++localPositiveCapMarkers;
        }
        else if (axialComponent < -0.5)
        {
            ++localNegativeCapMarkers;
        }
        else
        {
            ++localSideMarkers;
        }
    }
    const label sideMarkers = returnReduce
    (
        localSideMarkers,
        sumOp<label>()
    );
    const label positiveCapMarkers = returnReduce
    (
        localPositiveCapMarkers,
        sumOp<label>()
    );
    const label negativeCapMarkers = returnReduce
    (
        localNegativeCapMarkers,
        sumOp<label>()
    );
    const scalar sideMarkerArea =
        4.0*constant::mathematical::pi
       *geometry_.radius()*geometry_.halfLength()
       /max(sideMarkers, label(1));
    const scalar capArea =
        constant::mathematical::pi*sqr(geometry_.radius());
    const scalar positiveCapMarkerArea =
        capArea/max(positiveCapMarkers, label(1));
    const scalar negativeCapMarkerArea =
        capArea/max(negativeCapMarkers, label(1));
    forAll(markerAreaWeights_, markeri)
    {
        const scalar axialComponent =
            markerNormals_[markeri] & cylinderAxis;
        markerAreaWeights_[markeri] = axialComponent > 0.5
          ? positiveCapMarkerArea
          : axialComponent < -0.5
          ? negativeCapMarkerArea
          : sideMarkerArea;
    }
    const centredCPCCellToCellStencilObject& parallelStencil =
        centredCPCCellToCellStencilObject::New(mesh_);
    const mapDistribute& cellMap = parallelStencil.map();
    const pointField distributedCentres
    (
        distributeCellField<vector>(cellMap, mesh_.C())
    );
    volScalarField localAvailability
    (
        IOobject
        (
            "markerCompactAvailability",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar(dimless, 1.0)
    );
    localAvailability.correctBoundaryConditions();
    const scalarField distributedAvailability
    (
        distributeCellField<scalar>(cellMap, localAvailability)
    );

    markerStencils_.setSize(markerSeedCells_.size());
    markerWeights_.setSize(markerSeedCells_.size());
    //
    label emptyStencils = 0;
    label affineFallbacks = 0;
    label negativeWeightMarkers = 0;
    scalar localMaximumMomentError = 0.0;
    forAll(markerSeedCells_, markeri)
    {
        const label seed = markerSeedCells_[markeri];
        const scalar cellWidth = pow
        (
            max(mesh_.V()[seed], ROOTVSMALL),
            scalar(1.0/3.0)
        );
        const scalar supportRadius =
            markerSupportWidth_*cellWidth;
        DynamicList<label> stencil;
        DynamicList<scalar> weights;
        //
        for (const label celli : parallelStencil.stencil()[seed])
        {
            if
            (
                celli < 0
             || celli >= distributedCentres.size()
             || distributedAvailability[celli] < 0.5
            )
            {
                FatalErrorInFunction
                    << "Marker stencil references unavailable compact "
                    << "cell data: marker = " << markeri
                    << ", seed cell = " << seed
                    << ", compact label = " << celli
                    << ", construct size = "
                    << distributedCentres.size()
                    << ", availability = "
                    <<
                    (
                        celli >= 0
                     && celli < distributedAvailability.size()
                      ? distributedAvailability[celli]
                      : -1.0
                    )
                    << exit(FatalError);
            }
            const scalar centreDistance = mag
            (
                distributedCentres[celli] - markerPoints_[markeri]
            );
            if (centreDistance >= supportRadius)
            {
                continue;
            }
            const scalar q = centreDistance/supportRadius;
            stencil.append(celli);
            weights.append
            (
                0.5*(1.0 + cos(constant::mathematical::pi*q))
            );
        }
        if (stencil.empty())
        {
            stencil.append(seed);
            weights.append(1.0);
            ++emptyStencils;
        }
        const scalarList compactWeights(weights);
        if
        (
            !affineCorrectWeights
            (
                markerPoints_[markeri],
                distributedCentres,
                stencil,
                weights
            )
        )
        {
            weights = compactWeights;
            scalar weightSum = 0.0;
            for (const scalar weight : weights)
            {
                weightSum += weight;
            }
            forAll(weights, i)
            {
                weights[i] /= max(weightSum, ROOTVSMALL);
            }
            ++affineFallbacks;
        }
        bool hasNegativeWeight = false;
        point interpolatedPoint(point::zero);
        forAll(weights, i)
        {
            hasNegativeWeight = hasNegativeWeight || weights[i] < 0.0;
            interpolatedPoint +=
                weights[i]*distributedCentres[stencil[i]];
        }
        negativeWeightMarkers += hasNegativeWeight ? 1 : 0;
        localMaximumMomentError = max
        (
            localMaximumMomentError,
            mag(interpolatedPoint - markerPoints_[markeri])
        );
        stencil.shrink();
        weights.shrink();
        markerStencils_[markeri].transfer(stencil);
        markerWeights_[markeri].transfer(weights);
    }
    const label markerCount = returnReduce
    (
        markerSeedCells_.size(),
        sumOp<label>()
    );
    const label fallbackCount = returnReduce
    (
        emptyStencils,
        sumOp<label>()
    );
    const label affineFallbackCount = returnReduce
    (
        affineFallbacks,
        sumOp<label>()
    );
    const label negativeWeightMarkerCount = returnReduce
    (
        negativeWeightMarkers,
        sumOp<label>()
    );
    const scalar maximumMomentError = returnReduce
    (
        localMaximumMomentError,
        maxOp<scalar>()
    );
    Info<< "IBM marker operator: markers = "
        << markerCount
        << ", fallback stencils = "
        << fallbackCount
        << ", affine fallbacks = "
        << affineFallbackCount
        << ", negative-weight markers = "
        << negativeWeightMarkerCount
        << ", max first-moment error = "
        << maximumMomentError
        << nl;
}
//
void DirectForcingIbm::beginTimeStep()
{
#if ICHIBOUNDARYFOAM_ENABLE_LOAD_CALCULATION
    beginLoadCycle();
#endif
}
//
void DirectForcingIbm::beginLoadCycle()
{
#if ICHIBOUNDARYFOAM_ENABLE_LOAD_CALCULATION
    constraintImpulseDensity_ = dimensionedVector
    (
        "zero",
        dimDensity*dimVelocity,
        vector::zero
    );
    loadCorrectionCount_ = 0;
#endif
}
//
vector DirectForcingIbm::constrainedRelativeVelocity
(
    const vector& relativeVelocity,
    const label celli
) const
{
    if (velocityBoundaryCondition_ == "freeSlip")
    {
        const vector& normal =
            geometry_.surfaceNormal().primitiveField()[celli];
        return (relativeVelocity & normal)*normal;
    }
    return relativeVelocity;
}
//
vector DirectForcingIbm::constrainedMarkerRelativeVelocity
(
    const vector& relativeVelocity,
    const vector& normal
) const
{
    if (velocityBoundaryCondition_ == "freeSlip")
    {
        return (relativeVelocity & normal)*normal;
    }
    return relativeVelocity;
}
//
void DirectForcingIbm::constrainDiffuse
(
    volVectorField& predictor,
    volVectorField& correction
) const
{
    vectorField& correctionValues = correction.primitiveFieldRef();
    const vectorField& predictorValues = predictor.primitiveField();
    const vectorField& bodyVelocity =
        geometry_.surfaceVelocity().primitiveField();
    const scalarField& fraction = solidFraction_.primitiveField();
    forAll(correctionValues, celli)
    {
        correctionValues[celli] = fraction[celli]
           *constrainedRelativeVelocity
            (
                bodyVelocity[celli] - predictorValues[celli],
                celli
            );
    }
    correction.correctBoundaryConditions();
}
//
void DirectForcingIbm::constrainSharp
(
    volVectorField& predictor,
    volVectorField& correction
) const
{
    vectorField& correctionValues = correction.primitiveFieldRef();
    const vectorField& predictorValues = predictor.primitiveField();
    const vectorField& bodyVelocity =
        geometry_.surfaceVelocity().primitiveField();
    const scalarField& distance =
        geometry_.signedDistance().primitiveField();
    const scalarField& type = cellType_.primitiveField();
    const vectorField& centres = mesh_.cellCentres();
    //
    correctionValues = vector::zero;
    forAll(correctionValues, celli)
    {
        if (type[celli] == 2.0)
        {
            correctionValues[celli] =
                constrainedRelativeVelocity
                (
                    bodyVelocity[celli] - predictorValues[celli],
                    celli
                );
            continue;
        }
        if (type[celli] != 1.0 || distance[celli] <= 0.0)
        {
            continue;
        }
        //
        vector weightedTarget(vector::zero);
        scalar weightSum = 0.0;
        for
        (
            const label neighbour
          : alphaExtensionStencil_[celli]
        )
        {
            if (distance[neighbour] <= distance[celli])
            {
                continue;
            }
            //
            const scalar ratio = min
            (
                distance[celli]
               /max(distance[neighbour], ROOTVSMALL),
                sharpInterpolationLimit_
            );
            const vector target =
                (1.0 - ratio)*bodyVelocity[celli]
              + ratio*predictorValues[neighbour];
            const scalar centreDistance = max
            (
                mag(centres[celli] - centres[neighbour]),
                ROOTVSMALL
            );
            const scalar weight = 1.0/sqr(centreDistance);
            weightedTarget += weight*target;
            weightSum += weight;
        }
        //
        const vector target = weightSum > VSMALL
          ? weightedTarget/weightSum
          : bodyVelocity[celli];
        correctionValues[celli] = constrainedRelativeVelocity
        (
            target - predictorValues[celli],
            celli
        );
    }
    correction.correctBoundaryConditions();
}
//
void DirectForcingIbm::constrainMarkers
(
    volVectorField& predictor,
    volVectorField& correction
) const
{
    const vectorField& predictorValues = predictor.primitiveField();
    vectorField& correctionValues = correction.primitiveFieldRef();
    const vectorField interpolatedVelocity
    (
        interpolateMarkers(predictorValues)
    );
    vectorField rightHandSide(markerPoints_.size(), vector::zero);
    forAll(rightHandSide, markeri)
    {
        rightHandSide[markeri] =
            constrainedMarkerRelativeVelocity
            (
                markerVelocities_[markeri]
              - interpolatedVelocity[markeri],
                markerNormals_[markeri]
            );
    }
    vectorField multiplier(markerPoints_.size(), vector::zero);
    vectorField residual(rightHandSide);
    scalarField inverseDiagonal(markerPoints_.size(), 0.0);
    forAll(inverseDiagonal, markeri)
    {
        scalar diagonal = 0.0;
        for (const scalar weight : markerWeights_[markeri])
        {
            diagonal += sqr(weight);
        }
        inverseDiagonal[markeri] = 1.0/max(diagonal, VSMALL);
    }
    vectorField preconditionedResidual(residual);
    forAll(preconditionedResidual, markeri)
    {
        preconditionedResidual[markeri] *= inverseDiagonal[markeri];
    }
    vectorField direction(preconditionedResidual);
    scalar residualNormSqr = globalMarkerDot(residual, residual);
    scalar residualPreconditionedDot =
        globalMarkerDot(residual, preconditionedResidual);
    const scalar initialNorm = sqrt(max(residualNormSqr, scalar(0)));
    const scalar targetNorm = markerConstraintTolerance_
       *max(initialNorm, scalar(1));
    label iteration = 0;
    bool converged = initialNorm <= targetNorm;
    for
    (
        iteration = 0;
        !converged && iteration < markerConstraintMaxIterations_;
        ++iteration
    )
    {
        const vectorField schurDirection(applyMarkerSchur(direction));
        const scalar denominator =
            globalMarkerDot(direction, schurDirection);
        if (denominator <= VSMALL)
        {
            break;
        }
        const scalar step = residualPreconditionedDot/denominator;
        multiplier += step*direction;
        residual -= step*schurDirection;
        const scalar nextResidualNormSqr =
            globalMarkerDot(residual, residual);
        converged = sqrt(max(nextResidualNormSqr, scalar(0)))
         <= targetNorm;
        if (!converged)
        {
            preconditionedResidual = residual;
            forAll(preconditionedResidual, markeri)
            {
                preconditionedResidual[markeri] *=
                    inverseDiagonal[markeri];
            }
            const scalar nextResidualPreconditionedDot =
                globalMarkerDot(residual, preconditionedResidual);
            const scalar beta = nextResidualPreconditionedDot
               /max(residualPreconditionedDot, VSMALL);
            direction = preconditionedResidual + beta*direction;
            residualPreconditionedDot =
                nextResidualPreconditionedDot;
        }
        residualNormSqr = nextResidualNormSqr;
    }
    const scalar finalNorm = sqrt(max(residualNormSqr, scalar(0)));
    Info<< "Marker constraint solve: iterations = " << iteration
        << ", initial/final residual = " << initialNorm << '/'
        << finalNorm << ", target = " << targetNorm
        << ", converged = " << (converged ? "yes" : "no") << nl;
    if (!converged && requireMarkerConstraintConvergence_)
    {
        FatalErrorInFunction
            << "Marker constraint solve failed to converge: residual = "
            << finalNorm << ", target = " << targetNorm
            << exit(FatalError);
    }
    spreadMarkers(multiplier, correctionValues, false);
    correction.correctBoundaryConditions();
}

vectorField DirectForcingIbm::applyMarkerSchur
(
    const vectorField& markerValues
) const
{
    vectorField eulerianValues;
    spreadMarkers(markerValues, eulerianValues, false);
    vectorField result(interpolateMarkers(eulerianValues));
    forAll(result, markeri)
    {
        result[markeri] = constrainedMarkerRelativeVelocity
        (
            result[markeri],
            markerNormals_[markeri]
        );
    }
    return result;
}
//
vectorField DirectForcingIbm::interpolateMarkers
(
    const vectorField& eulerianValues
) const
{
    if (eulerianValues.size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "Eulerian input size " << eulerianValues.size()
            << " does not match local cell count " << mesh_.nCells()
            << exit(FatalError);
    }
    const centredCPCCellToCellStencilObject& parallelStencil =
        centredCPCCellToCellStencilObject::New(mesh_);
    volVectorField eulerianField
    (
        IOobject
        (
            "markerInterpolationSource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector(dimless, vector::zero)
    );
    eulerianField.primitiveFieldRef() = eulerianValues;
    forAll(eulerianField.boundaryField(), patchi)
    {
        if (!eulerianField.boundaryField()[patchi].coupled())
        {
            eulerianField.boundaryFieldRef()[patchi] =
                eulerianField.boundaryField()[patchi].patchInternalField();
        }
    }
    eulerianField.correctBoundaryConditions();
    const vectorField distributedValues
    (
        distributeCellField<vector>(parallelStencil.map(), eulerianField)
    );
    vectorField markerValues(markerPoints_.size(), vector::zero);
    forAll(markerValues, markeri)
    {
        const labelList& stencil = markerStencils_[markeri];
        const scalarList& weights = markerWeights_[markeri];
        forAll(stencil, i)
        {
            markerValues[markeri] +=
                weights[i]*distributedValues[stencil[i]];
        }
    }
    return markerValues;
}
//
void DirectForcingIbm::spreadMarkers
(
    const vectorField& markerValues,
    vectorField& eulerianValues,
    const bool normalizeOverlap
) const
{
    if (markerValues.size() != markerPoints_.size())
    {
        FatalErrorInFunction
            << "Marker input size " << markerValues.size()
            << " does not match local marker count "
            << markerPoints_.size()
            << exit(FatalError);
    }
    const centredCPCCellToCellStencilObject& parallelStencil =
        centredCPCCellToCellStencilObject::New(mesh_);
    const mapDistribute& cellMap = parallelStencil.map();
    vectorField distributedValues(cellMap.constructSize(), vector::zero);
    scalarField spreadWeight(cellMap.constructSize(), 0.0);
    forAll(markerValues, markeri)
    {
        const labelList& stencil = markerStencils_[markeri];
        const scalarList& weights = markerWeights_[markeri];
        forAll(stencil, i)
        {
            const label celli = stencil[i];
            distributedValues[celli] +=
                weights[i]*markerValues[markeri];
            spreadWeight[celli] += weights[i];
        }
    }
    mapDistributeBase::distribute
    (
        UPstream::defaultCommsType, cellMap.schedule(), mesh_.nCells(),
        cellMap.constructMap(), cellMap.constructHasFlip(),
        cellMap.subMap(), cellMap.subHasFlip(), distributedValues,
        vector::zero, plusEqOp<vector>(), flipOp(),
        UPstream::msgType(), cellMap.comm()
    );
    mapDistributeBase::distribute
    (
        UPstream::defaultCommsType, cellMap.schedule(), mesh_.nCells(),
        cellMap.constructMap(), cellMap.constructHasFlip(),
        cellMap.subMap(), cellMap.subHasFlip(), spreadWeight,
        scalar(0), plusEqOp<scalar>(), flipOp(),
        UPstream::msgType(), cellMap.comm()
    );
    eulerianValues.transfer(distributedValues);
    if (!normalizeOverlap) return;
    forAll(eulerianValues, celli)
    {
        if (spreadWeight[celli] > VSMALL)
        {
            eulerianValues[celli] /= spreadWeight[celli];
        }
    }
}
//
void DirectForcingIbm::constrainPredictor
(
    volVectorField& predictor,
    const volScalarField& density,
    const scalar deltaT
)
{
    if (!enabled_ || deltaT <= VSMALL)
    {
        return;
    }
    volVectorField correction
    (
        IOobject
        (
            "immersedPredictorCorrection",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("zero", dimVelocity, vector::zero)
    );
    if (couplingMethod_ == "markerBased")
    {
        constrainMarkers(predictor, correction);
    }
    else if (forcingMethod_ == "sharpReconstruction")
    {
        constrainSharp(predictor, correction);
    }
    else
    {
        constrainDiffuse(predictor, correction);
    }
    predictor += correction;
    predictor.correctBoundaryConditions();
#if ICHIBOUNDARYFOAM_ENABLE_LOAD_CALCULATION
    constraintImpulseDensity_ += density*correction;
    constraintImpulseDensity_.correctBoundaryConditions();
    ++loadCorrectionCount_;
#endif
}
//
void DirectForcingIbm::reportResiduals
(
    const volVectorField& velocity,
    const surfaceScalarField& flux,
    const volScalarField& alpha
) const
{
    if (!enabled_)
    {
        return;
    }
    const volScalarField divergence
    (
        IOobject
        (
            "immersedContinuityResidual",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mag(fvc::div(flux))
    );
    Info<< "IBM projection residuals: max "
        << velocityBoundaryCondition_
        << " constraint residual = "
        << maxSlip(velocity)
        << ", max divergence = "
        << gMax(divergence.primitiveField())
        << nl;
#if ICHIBOUNDARYFOAM_ENABLE_SLIP_SAMPLING
    if (slipSamplingBandWidth_ <= VSMALL)
    {
        return;
    }
    const vectorField relativeVelocity
    (
        velocity.primitiveField()
      - geometry_.surfaceVelocity().primitiveField()
    );
    const vectorField& normal =
        geometry_.surfaceNormal().primitiveField();
    const scalarField& distance =
        geometry_.signedDistance().primitiveField();
    const scalarField& alphaValues = alpha.primitiveField();
    vectorField& slipValues = slipVelocity_.primitiveFieldRef();
    scalarField& magnitudeValues =
        slipMagnitude_.primitiveFieldRef();
    scalarField& normalValues = normalSlip_.primitiveFieldRef();
    scalarField& tangentialValues =
        tangentialSlip_.primitiveFieldRef();
    scalarField& maskValues =
        slipSamplingMask_.primitiveFieldRef();
    slipValues = vector::zero;
    magnitudeValues = 0.0;
    normalValues = 0.0;
    tangentialValues = 0.0;
    maskValues = 0.0;
    scalar localMaximum = -GREAT;
    scalar localAlpha = 0.0;
    scalar localDistance = 0.0;
    point localLocation(point::zero);
    scalar localLiquidMaximum = 0.0;
    scalar localGasMaximum = 0.0;
    scalar localNormalMaximum = 0.0;
    scalar localTangentialMaximum = 0.0;
    label localSampleCount = 0;
    label localLiquidCount = 0;
    label localGasCount = 0;
    forAll(relativeVelocity, celli)
    {
        if (mag(distance[celli]) > slipSamplingBandWidth_)
        {
            continue;
        }
        const scalar signedNormal =
            relativeVelocity[celli] & normal[celli];
        const vector tangential =
            relativeVelocity[celli] - signedNormal*normal[celli];
        const scalar slipMag = mag(relativeVelocity[celli]);
        const scalar normalMag = mag(signedNormal);
        const scalar tangentialMag = mag(tangential);
        slipValues[celli] = relativeVelocity[celli];
        magnitudeValues[celli] = slipMag;
        normalValues[celli] = normalMag;
        tangentialValues[celli] = tangentialMag;
        maskValues[celli] = 1.0;
        ++localSampleCount;
        localNormalMaximum = max(localNormalMaximum, normalMag);
        localTangentialMaximum = max
        (
            localTangentialMaximum,
            tangentialMag
        );
        if (alphaValues[celli] >= slipLiquidThreshold_)
        {
            ++localLiquidCount;
            localLiquidMaximum = max(localLiquidMaximum, slipMag);
        }
        else
        {
            ++localGasCount;
            localGasMaximum = max(localGasMaximum, slipMag);
        }
        if (slipMag > localMaximum)
        {
            localMaximum = slipMag;
            localAlpha = alphaValues[celli];
            localDistance = distance[celli];
            localLocation = mesh_.C()[celli];
        }
    }
    slipVelocity_.correctBoundaryConditions();
    slipMagnitude_.correctBoundaryConditions();
    normalSlip_.correctBoundaryConditions();
    tangentialSlip_.correctBoundaryConditions();
    slipSamplingMask_.correctBoundaryConditions();
    List<scalar> maximumByProc(Pstream::nProcs(), -GREAT);
    List<scalar> alphaByProc(Pstream::nProcs(), 0.0);
    List<scalar> distanceByProc(Pstream::nProcs(), 0.0);
    List<point> locationByProc
    (
        Pstream::nProcs(),
        point::zero
    );
    const label proc = Pstream::myProcNo();
    maximumByProc[proc] = localMaximum;
    alphaByProc[proc] = localAlpha;
    distanceByProc[proc] = localDistance;
    locationByProc[proc] = localLocation;
    Pstream::gatherList(maximumByProc);
    Pstream::gatherList(alphaByProc);
    Pstream::gatherList(distanceByProc);
    Pstream::gatherList(locationByProc);
    const label sampleCount = returnReduce
    (
        localSampleCount,
        sumOp<label>()
    );
    const label liquidCount = returnReduce
    (
        localLiquidCount,
        sumOp<label>()
    );
    const label gasCount = returnReduce
    (
        localGasCount,
        sumOp<label>()
    );
    const scalar liquidMaximum = returnReduce
    (
        localLiquidMaximum,
        maxOp<scalar>()
    );
    const scalar gasMaximum = returnReduce
    (
        localGasMaximum,
        maxOp<scalar>()
    );
    const scalar normalMaximum = returnReduce
    (
        localNormalMaximum,
        maxOp<scalar>()
    );
    const scalar tangentialMaximum = returnReduce
    (
        localTangentialMaximum,
        maxOp<scalar>()
    );
    if (Pstream::master())
    {
        if (sampleCount == 0)
        {
            Info<< "IBM surface-band slip: no cell centres lie in "
                << "the configured SDF band" << nl;
            return;
        }
        label maximumProc = 0;
        for (label proci = 1; proci < Pstream::nProcs(); ++proci)
        {
            if (maximumByProc[proci] > maximumByProc[maximumProc])
            {
                maximumProc = proci;
            }
        }
        Info<< "IBM surface-band slip: samples = "
            << sampleCount
            << ", liquid/gas = "
            << liquidCount
            << '/'
            << gasCount
            << ", max all/liquid/gas = "
            << maximumByProc[maximumProc]
            << '/'
            << liquidMaximum
            << '/'
            << gasMaximum
            << nl;
        Info<< "IBM surface-band components: max normal/tangential = "
            << normalMaximum
            << '/'
            << tangentialMaximum
            << ", location = "
            << locationByProc[maximumProc]
            << ", alpha = "
            << alphaByProc[maximumProc]
            << ", signed distance = "
            << distanceByProc[maximumProc]
            << nl;
    }
#endif
}
//
scalar DirectForcingIbm::maxSlip
(
    const volVectorField& velocity
) const
{
    if (!enabled_)
    {
        return 0.0;
    }
    const vectorField& velocityValues = velocity.primitiveField();
    if (couplingMethod_ == "markerBased")
    {
        const vectorField interpolatedVelocity
        (
            interpolateMarkers(velocityValues)
        );
        scalar localMaximum = 0.0;
        forAll(markerPoints_, markeri)
        {
            const vector residual =
                constrainedMarkerRelativeVelocity
                (
                    interpolatedVelocity[markeri]
                  - markerVelocities_[markeri],
                    markerNormals_[markeri]
                );
            localMaximum = max(localMaximum, mag(residual));
        }
        return returnReduce(localMaximum, maxOp<scalar>());
    }
    const scalarField& fraction = solidFraction_.primitiveField();
    const vectorField& bodyVelocity =
        geometry_.surfaceVelocity().primitiveField();
    scalar localMaximum = 0.0;
    forAll(velocityValues, celli)
    {
        const vector residual = constrainedRelativeVelocity
        (
            velocityValues[celli] - bodyVelocity[celli],
            celli
        );
        localMaximum = max
        (
            localMaximum,
            fraction[celli]*mag(residual)
        );
    }
    return returnReduce(localMaximum, maxOp<scalar>());
}
//
bool DirectForcingIbm::multiDirectForcing() const
{
    return enabled_ && multiDirectForcing_;
}
//
label DirectForcingIbm::nMultiDirectForcingCorrectors() const
{
    return nMultiDirectForcingCorrectors_;
}
//
scalar DirectForcingIbm::slipTolerance() const
{
    return slipTolerance_;
}
//
void DirectForcingIbm::reportLoads
(
    const scalar deltaT
) const
{
#if ICHIBOUNDARYFOAM_ENABLE_LOAD_CALCULATION
    if (!enabled_ || deltaT <= VSMALL)
    {
        return;
    }
    const vectorField cellImpulse
    (
        mesh_.V()*constraintImpulseDensity_.primitiveField()
    );
    const vector bodyImpulse = -gSum(cellImpulse);
    const vector bodyForce = bodyImpulse/deltaT;
    vectorField cellAngularImpulse(mesh_.nCells(), vector::zero);
    forAll(cellAngularImpulse, celli)
    {
        cellAngularImpulse[celli] =
            (mesh_.C()[celli] - geometry_.centre())
          ^ cellImpulse[celli];
    }
    const vector bodyAngularImpulse = -gSum(cellAngularImpulse);
    const vector bodyMoment = bodyAngularImpulse/deltaT;
    Info<< "IBM constraint load: force = "
        << bodyForce
        << ", moment about " << geometry_.centre() << " = "
        << bodyMoment
        << ", impulse = " << bodyImpulse
        << ", angular impulse = " << bodyAngularImpulse
        << ", forcing corrections = " << loadCorrectionCount_
        << nl;
#else
    (void)deltaT;
#endif
}
//
void DirectForcingIbm::reportPhaseVolume
(
    const volScalarField& alpha
) const
{
    if (!enabled_)
    {
        return;
    }
    const scalarField accessibleFraction
    (
        scalar(1) - solidFraction_.primitiveField()
    );
    const scalar accessibleVolume = gSum
    (
        accessibleFraction*mesh_.V()
    );
    const scalar phaseVolume = gSum
    (
        accessibleFraction
       *alpha.primitiveField()
       *mesh_.V()
    );
    Info<< "IBM accessible phase volume: fluid = "
        << accessibleVolume << ", phase-1 = "
        << phaseVolume << nl;
}
//
bool DirectForcingIbm::enabled() const
{
    return enabled_;
}
//
const word& DirectForcingIbm::couplingMethod() const
{
    return couplingMethod_;
}
//
const word& DirectForcingIbm::forcingMethod() const
{
    return forcingMethod_;
}
//
const word& DirectForcingIbm::velocityBoundaryCondition() const
{
    return velocityBoundaryCondition_;
}
//
const pointField& DirectForcingIbm::markerPoints() const
{
    return markerPoints_;
}
//
const vectorField& DirectForcingIbm::markerNormals() const
{
    return markerNormals_;
}
//
const vectorField& DirectForcingIbm::markerVelocities() const
{
    return markerVelocities_;
}
//
const scalarField& DirectForcingIbm::markerAreaWeights() const
{
    return markerAreaWeights_;
}
//
const labelListList& DirectForcingIbm::markerStencils() const
{
    return markerStencils_;
}
//
const List<scalarList>& DirectForcingIbm::markerWeights() const
{
    return markerWeights_;
}
//
vectorField DirectForcingIbm::interpolateMarkerField
(
    const vectorField& eulerianValues
) const
{
    return interpolateMarkers(eulerianValues);
}
//
void DirectForcingIbm::spreadMarkerField
(
    const vectorField& markerValues,
    vectorField& eulerianValues,
    const bool normalizeOverlap
) const
{
    spreadMarkers(markerValues, eulerianValues, normalizeOverlap);
}
//
void DirectForcingIbm::spreadMarkerForceField
(
    const vectorField& markerForces,
    vectorField& eulerianForceDensity
) const
{
    if (markerForces.size() != markerPoints_.size())
    {
        FatalErrorInFunction
            << "Marker force size " << markerForces.size()
            << " does not match local marker count "
            << markerPoints_.size() << exit(FatalError);
    }
    const centredCPCCellToCellStencilObject& parallelStencil =
        centredCPCCellToCellStencilObject::New(mesh_);
    const mapDistribute& cellMap = parallelStencil.map();
    vectorField distributedForce(cellMap.constructSize(), vector::zero);
    volScalarField cellVolumes
    (
        IOobject
        (
            "markerCellVolumes",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar(dimVolume, 0.0)
    );
    cellVolumes.primitiveFieldRef() = mesh_.V();
    forAll(cellVolumes.boundaryField(), patchi)
    {
        if (!cellVolumes.boundaryField()[patchi].coupled())
        {
            cellVolumes.boundaryFieldRef()[patchi] =
                cellVolumes.boundaryField()[patchi].patchInternalField();
        }
    }
    cellVolumes.correctBoundaryConditions();
    const scalarField distributedVolumes
    (
        distributeCellField<scalar>(cellMap, cellVolumes)
    );
    forAll(markerForces, markeri)
    {
        const labelList& stencil = markerStencils_[markeri];
        const scalarList& weights = markerWeights_[markeri];
        forAll(stencil, i)
        {
            const label celli = stencil[i];
            distributedForce[celli] +=
                weights[i]*markerForces[markeri]
               *markerAreaWeights_[markeri]/distributedVolumes[celli];
        }
    }
    mapDistributeBase::distribute
    (
        UPstream::defaultCommsType, cellMap.schedule(), mesh_.nCells(),
        cellMap.constructMap(), cellMap.constructHasFlip(),
        cellMap.subMap(), cellMap.subHasFlip(), distributedForce,
        vector::zero, plusEqOp<vector>(), flipOp(),
        UPstream::msgType(), cellMap.comm()
    );
    eulerianForceDensity.transfer(distributedForce);
}
//
const volScalarField& DirectForcingIbm::solidFraction() const
{
    return solidFraction_;
}
}
}
