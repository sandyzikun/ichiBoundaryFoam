#include "fvCFD.H"
#include "immersedGeometry.H"
#include "directForcingIbm.H"

using namespace Foam;
using namespace Foam::ichiBoundaryFoam;

namespace
{
scalar globalMaximum(const scalar localValue)
{
    return returnReduce(localValue, maxOp<scalar>());
}

scalar globalSum(const scalar localValue)
{
    return returnReduce(localValue, sumOp<scalar>());
}

bool selected(const word& requested, const word& testName)
{
    return requested == "all" || requested == testName;
}

void reportResult
(
    const word& name,
    const scalar error,
    const scalar tolerance,
    const bool checked,
    label& failures
)
{
    Info<< "markerOperatorTest " << name << ": error = " << error;
    if (!checked)
    {
        Info<< " (characterization only)" << nl;
        return;
    }

    const bool passed = error <= tolerance;
    Info<< ", tolerance = " << tolerance
        << ", result = " << (passed ? "PASS" : "FAIL") << nl;
    if (!passed)
    {
        ++failures;
    }
}

vector analyticLinearField(const point& location)
{
    return vector
    (
        location.x() + 2.0*location.y(),
        location.y() - location.z(),
        2.0*location.z() + location.x()
    );
}

vector analyticMarkerLoad(const point& location)
{
    return vector
    (
        1.0 + location.x(),
        -0.5 + location.y(),
        0.25 + location.z()
    );
}

vector analyticQuadraticField(const point& location)
{
    return vector
    (
        sqr(location.x()) + location.y()*location.z(),
        sqr(location.y()) - location.x()*location.z(),
        sqr(location.z()) + location.x()*location.y()
    );
}
}

int main(int argc, char* argv[])
{
    argList::addNote
    (
        "Verify the production marker interpolation/spreading operator "
        "without solving any flow equation"
    );
    argList::addOption
    (
        "test",
        "name",
        "geometry, constant, linear, quadratic, zero, transpose, "
        "conservation, weightedAdjoint, constraintSolve, or all"
    );
    argList::addBoolOption
    (
        "strict",
        "apply the configured tolerance to the linear-field test"
    );
    argList::addOption
    (
        "tolerance",
        "scalar",
        "tolerance for exact-consistency tests; default 1e-12"
    );
    argList::addOption
    (
        "geometryTolerance",
        "scalar",
        "surface projection and unit-normal tolerance; default 1e-10"
    );
    argList::addOption
    (
        "linearTolerance",
        "scalar",
        "linear interpolation tolerance used with -strict; default 1e-2"
    );
    argList::addOption
    (
        "constraintTolerance",
        "scalar",
        "post-correction marker residual tolerance; default 1e-8"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const word requested(args.getOrDefault<word>("test", "all"));
    const bool strict(args.found("strict"));
    const scalar tolerance
    (
        args.getOrDefault<scalar>("tolerance", 1e-12)
    );
    const scalar geometryTolerance
    (
        args.getOrDefault<scalar>("geometryTolerance", 1e-10)
    );
    const scalar linearTolerance
    (
        args.getOrDefault<scalar>("linearTolerance", 1e-2)
    );
    const scalar constraintTolerance
    (
        args.getOrDefault<scalar>("constraintTolerance", 1e-8)
    );

    if
    (
        requested != "all"
     && requested != "geometry"
     && requested != "constant"
     && requested != "linear"
     && requested != "quadratic"
     && requested != "zero"
     && requested != "transpose"
     && requested != "conservation"
     && requested != "weightedAdjoint"
     && requested != "constraintSolve"
    )
    {
        FatalErrorInFunction
            << "Unknown test '" << requested << "'. Supported tests: "
            << "all, geometry, constant, linear, quadratic, zero, "
            << "transpose, conservation, weightedAdjoint, constraintSolve"
            << exit(FatalError);
    }

    ImmersedGeometry geometry(mesh);
    geometry.update(runTime.value());
    IOdictionary properties
    (
        IOobject
        (
            "ichiBoundaryGeometry",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    DirectForcingIbm ibm
    (
        mesh,
        geometry,
        properties.subDict("ibm")
    );
    ibm.classify();

    const pointField& markers = ibm.markerPoints();
    const label globalMarkerCount = returnReduce
    (
        markers.size(),
        sumOp<label>()
    );
    Info<< "markerOperatorTest: local/global markers = "
        << markers.size() << '/' << globalMarkerCount << nl;
    if (globalMarkerCount == 0)
    {
        FatalErrorInFunction
            << "No markers were generated. Check the mesh, geometry, and "
            << "markerBased case configuration." << exit(FatalError);
    }

    label failures = 0;

    if (selected(requested, "geometry"))
    {
        scalar localDistanceError = 0.0;
        scalar localNormalError = 0.0;
        scalar localStoredNormalError = 0.0;
        forAll(markers, markeri)
        {
            scalar distance = 0.0;
            vector normal(vector::zero);
            vector velocity(vector::zero);
            geometry.evaluatePoint
            (
                markers[markeri],
                distance,
                normal,
                velocity
            );
            localDistanceError = max(localDistanceError, mag(distance));
            localNormalError = max
            (
                localNormalError,
                mag(mag(normal) - 1.0)
            );
            localStoredNormalError = max
            (
                localStoredNormalError,
                mag(normal - ibm.markerNormals()[markeri])
            );
        }
        reportResult
        (
            "geometry.surfaceDistance",
            globalMaximum(localDistanceError),
            geometryTolerance,
            true,
            failures
        );
        reportResult
        (
            "geometry.unitNormal",
            globalMaximum(localNormalError),
            geometryTolerance,
            true,
            failures
        );
        reportResult
        (
            "geometry.storedNormal",
            globalMaximum(localStoredNormalError),
            geometryTolerance,
            true,
            failures
        );
    }

    if (selected(requested, "constant"))
    {
        const vector expected(0.75, -1.25, 2.5);
        const vectorField cellValues(mesh.nCells(), expected);
        const vectorField markerValues
        (
            ibm.interpolateMarkerField(cellValues)
        );
        scalar localError = 0.0;
        for (const vector& value : markerValues)
        {
            localError = max(localError, mag(value - expected));
        }
        reportResult
        (
            "constantInterpolation",
            globalMaximum(localError),
            tolerance,
            true,
            failures
        );
    }

    if (selected(requested, "linear"))
    {
        vectorField cellValues(mesh.nCells(), vector::zero);
        forAll(cellValues, celli)
        {
            cellValues[celli] = analyticLinearField(mesh.C()[celli]);
        }
        const vectorField markerValues
        (
            ibm.interpolateMarkerField(cellValues)
        );
        scalar localMaximumError = 0.0;
        scalar localSquaredError = 0.0;
        forAll(markerValues, markeri)
        {
            const scalar error = mag
            (
                markerValues[markeri]
              - analyticLinearField(markers[markeri])
            );
            localMaximumError = max(localMaximumError, error);
            localSquaredError += sqr(error);
        }
        const scalar rmsError = Foam::sqrt
        (
            globalSum(localSquaredError)
           /scalar(max(globalMarkerCount, label(1)))
        );
        Info<< "markerOperatorTest linearInterpolation: RMS error = "
            << rmsError << nl;
        reportResult
        (
            "linearInterpolation.max",
            globalMaximum(localMaximumError),
            linearTolerance,
            strict,
            failures
        );
    }

    if (selected(requested, "zero"))
    {
        const vectorField zeroCells(mesh.nCells(), vector::zero);
        const vectorField interpolated
        (
            ibm.interpolateMarkerField(zeroCells)
        );
        scalar localInterpolationError = 0.0;
        for (const vector& value : interpolated)
        {
            localInterpolationError = max
            (
                localInterpolationError,
                mag(value)
            );
        }
        const vectorField zeroMarkers(markers.size(), vector::zero);
        vectorField spread;
        ibm.spreadMarkerField(zeroMarkers, spread, true);
        scalar localSpreadError = 0.0;
        for (const vector& value : spread)
        {
            localSpreadError = max(localSpreadError, mag(value));
        }
        reportResult
        (
            "zeroInterpolation",
            globalMaximum(localInterpolationError),
            tolerance,
            true,
            failures
        );
        reportResult
        (
            "zeroSpreading",
            globalMaximum(localSpreadError),
            tolerance,
            true,
            failures
        );
    }

    if (selected(requested, "constraintSolve"))
    {
        volVectorField velocity
        (
            IOobject
            (
                "markerConstraintTestVelocity",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedVector("zero", dimVelocity, vector::zero)
        );
        forAll(velocity.primitiveFieldRef(), celli)
        {
            velocity.primitiveFieldRef()[celli] =
                analyticQuadraticField(mesh.C()[celli]);
        }
        velocity.correctBoundaryConditions();
        volScalarField density
        (
            IOobject
            (
                "markerConstraintTestDensity",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("one", dimDensity, 1.0)
        );
        ibm.constrainPredictor(velocity, density, 1.0);
        reportResult
        (
            "constraintSolve.maxResidual",
            ibm.maxSlip(velocity),
            constraintTolerance,
            true,
            failures
        );
    }

    if (selected(requested, "quadratic"))
    {
        vectorField cellValues(mesh.nCells(), vector::zero);
        forAll(cellValues, celli)
        {
            cellValues[celli] = analyticQuadraticField(mesh.C()[celli]);
        }
        const vectorField markerValues
        (
            ibm.interpolateMarkerField(cellValues)
        );
        scalar localMaximumError = 0.0;
        scalar localSquaredError = 0.0;
        forAll(markerValues, markeri)
        {
            const scalar error = mag
            (
                markerValues[markeri]
              - analyticQuadraticField(markers[markeri])
            );
            localMaximumError = max(localMaximumError, error);
            localSquaredError += sqr(error);
        }
        const scalar rmsError = Foam::sqrt
        (
            globalSum(localSquaredError)
           /scalar(max(globalMarkerCount, label(1)))
        );
        Info<< "markerOperatorTest quadraticInterpolation: RMS error = "
            << rmsError << nl;
        reportResult
        (
            "quadraticInterpolation.max",
            globalMaximum(localMaximumError),
            0.0,
            false,
            failures
        );
    }

    if (selected(requested, "transpose"))
    {
        vectorField cellValues(mesh.nCells(), vector::zero);
        forAll(cellValues, celli)
        {
            cellValues[celli] = analyticLinearField(mesh.C()[celli]);
        }
        vectorField markerLoads(markers.size(), vector::zero);
        forAll(markerLoads, markeri)
        {
            markerLoads[markeri] = analyticMarkerLoad(markers[markeri]);
        }
        const vectorField interpolated
        (
            ibm.interpolateMarkerField(cellValues)
        );
        vectorField transposeSpread;
        ibm.spreadMarkerField(markerLoads, transposeSpread, false);

        scalar localMarkerProduct = 0.0;
        forAll(markerLoads, markeri)
        {
            localMarkerProduct +=
                interpolated[markeri] & markerLoads[markeri];
        }
        scalar localCellProduct = 0.0;
        forAll(cellValues, celli)
        {
            localCellProduct +=
                cellValues[celli] & transposeSpread[celli];
        }
        const scalar markerProduct = globalSum(localMarkerProduct);
        const scalar cellProduct = globalSum(localCellProduct);
        const scalar relativeError = mag(markerProduct - cellProduct)
           /max(max(mag(markerProduct), mag(cellProduct)), ROOTVSMALL);
        reportResult
        (
            "algebraicTranspose",
            relativeError,
            tolerance,
            true,
            failures
        );
        Info<< "markerOperatorTest: this is the unweighted algebraic "
            << "transpose check; the separately tested physical adjoint "
            << "uses marker areas and cell volumes." << nl;
    }

    if (selected(requested, "conservation"))
    {
        const scalarField& markerAreas = ibm.markerAreaWeights();
        const scalar expectedArea =
            4.0*constant::mathematical::pi
           *geometry.radius()*geometry.halfLength()
          + 2.0*constant::mathematical::pi*sqr(geometry.radius());
        scalar localArea = 0.0;
        vectorField markerForces(markers.size(), vector::zero);
        vector localMarkerForce(vector::zero);
        vector localMarkerTorque(vector::zero);
        scalar localAbsoluteMarkerForce = 0.0;
        forAll(markers, markeri)
        {
            markerForces[markeri] = analyticMarkerLoad(markers[markeri]);
            localArea += markerAreas[markeri];
            localMarkerForce +=
                markerAreas[markeri]*markerForces[markeri];
            localAbsoluteMarkerForce +=
                markerAreas[markeri]*mag(markerForces[markeri]);
            localMarkerTorque +=
                (markers[markeri] - geometry.centre())
              ^ (markerAreas[markeri]*markerForces[markeri]);
        }
        const scalar areaError = mag(globalSum(localArea) - expectedArea)
           /max(expectedArea, ROOTVSMALL);
        reportResult
        (
            "surfaceArea",
            areaError,
            tolerance,
            true,
            failures
        );

        vectorField forceDensity;
        ibm.spreadMarkerForceField(markerForces, forceDensity);
        vector localEulerianForce(vector::zero);
        vector localEulerianTorque(vector::zero);
        forAll(forceDensity, celli)
        {
            const vector cellForce =
                mesh.V()[celli]*forceDensity[celli];
            localEulerianForce += cellForce;
            localEulerianTorque +=
                (mesh.C()[celli] - geometry.centre()) ^ cellForce;
        }
        const vector markerForce = returnReduce
        (
            localMarkerForce,
            sumOp<vector>()
        );
        const vector markerTorque = returnReduce
        (
            localMarkerTorque,
            sumOp<vector>()
        );
        const vector eulerianForce = returnReduce
        (
            localEulerianForce,
            sumOp<vector>()
        );
        const vector eulerianTorque = returnReduce
        (
            localEulerianTorque,
            sumOp<vector>()
        );
        const scalar torqueScale = max
        (
            (geometry.radius() + geometry.halfLength())
           *globalSum(localAbsoluteMarkerForce),
            ROOTVSMALL
        );
        Info<< "markerOperatorTest conservation: marker/eulerian force = "
            << markerForce << '/' << eulerianForce
            << ", marker/eulerian torque = "
            << markerTorque << '/' << eulerianTorque << nl;
        reportResult
        (
            "forceConservation",
            mag(eulerianForce - markerForce)
           /max(mag(markerForce), ROOTVSMALL),
            tolerance,
            true,
            failures
        );
        reportResult
        (
            "torqueConservation",
            mag(eulerianTorque - markerTorque)/torqueScale,
            tolerance,
            true,
            failures
        );
    }

    if (selected(requested, "weightedAdjoint"))
    {
        vectorField cellValues(mesh.nCells(), vector::zero);
        forAll(cellValues, celli)
        {
            cellValues[celli] = analyticLinearField(mesh.C()[celli]);
        }
        vectorField markerForces(markers.size(), vector::zero);
        forAll(markerForces, markeri)
        {
            markerForces[markeri] = analyticMarkerLoad(markers[markeri]);
        }
        const vectorField interpolated
        (
            ibm.interpolateMarkerField(cellValues)
        );
        vectorField forceDensity;
        ibm.spreadMarkerForceField(markerForces, forceDensity);
        scalar localMarkerPower = 0.0;
        forAll(markerForces, markeri)
        {
            localMarkerPower +=
                ibm.markerAreaWeights()[markeri]
               *(interpolated[markeri] & markerForces[markeri]);
        }
        scalar localEulerianPower = 0.0;
        forAll(cellValues, celli)
        {
            localEulerianPower += mesh.V()[celli]
               *(cellValues[celli] & forceDensity[celli]);
        }
        const scalar markerPower = globalSum(localMarkerPower);
        const scalar eulerianPower = globalSum(localEulerianPower);
        reportResult
        (
            "weightedAdjoint",
            mag(eulerianPower - markerPower)
           /max(max(mag(eulerianPower), mag(markerPower)), ROOTVSMALL),
            tolerance,
            true,
            failures
        );
    }

    Info<< "markerOperatorTest summary: failures = " << failures << nl;
    return failures == 0 ? 0 : 1;
}
