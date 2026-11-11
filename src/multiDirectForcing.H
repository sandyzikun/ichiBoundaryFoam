#if ICHIBOUNDARYFOAM_ENABLE_MDF_PRESSURE_PROJECTION
if (ibm.multiDirectForcing())
{
    volScalarField pressureCorrection
    (
        IOobject
        (
            "immersedPressureCorrection",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        p_rgh
    );
    pressureCorrection = dimensionedScalar
    (
        "zero",
        p_rgh.dimensions(),
        0.0
    );
    mesh.setFluxRequired(pressureCorrection.name());
    //
    for
    (
        label ibmCorr = 0;
        ibmCorr < ibm.nMultiDirectForcingCorrectors();
        ++ibmCorr
    )
    {
        const scalar slipBefore = ibm.maxSlip(U);
        if
        (
            ibmCorr > 0
         && ibm.slipTolerance() > 0.0
         && slipBefore <= ibm.slipTolerance()
        )
        {
            Info<< "Multi-direct forcing converged after "
                << ibmCorr << " correctors, max slip = "
                << slipBefore << nl;
            break;
        }
        //
        volVectorField constrainedVelocity
        (
            IOobject
            (
                "immersedConstrainedVelocity",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            U
        );
        ibm.constrainPredictor
        (
            constrainedVelocity,
            rho,
            runTime.deltaTValue()
        );
        //
        surfaceScalarField constrainedFlux
        (
            "immersedConstrainedFlux",
            fvc::flux(constrainedVelocity)
        );
        MRF.makeRelative(constrainedFlux);
        if (pressureCorrection.needReference())
        {
            fvc::makeRelative
            (
                constrainedFlux,
                constrainedVelocity
            );
            adjustPhi
            (
                constrainedFlux,
                constrainedVelocity,
                pressureCorrection
            );
            fvc::makeAbsolute
            (
                constrainedFlux,
                constrainedVelocity
            );
        }
        //
        pressureCorrection = dimensionedScalar
        (
            "zero",
            p_rgh.dimensions(),
            0.0
        );
        constrainPressure
        (
            pressureCorrection,
            constrainedVelocity,
            constrainedFlux,
            rAUf,
            MRF
        );
        fvScalarMatrix pressureCorrectionEqn
        (
            fvm::laplacian(rAUf, pressureCorrection)
         == fvc::div(constrainedFlux)
        );
        pressureCorrectionEqn.setReference(pRefCell, 0.0);
        pressureCorrectionEqn.solve
        (
            "p_rghFinal"
        );
        //
        phi = constrainedFlux - pressureCorrectionEqn.flux();
        U = constrainedVelocity
          - rAU()*fvc::reconstruct
            (
                pressureCorrectionEqn.flux()/rAUf
            );
        U.correctBoundaryConditions();
        fvOptions.correct(U);
        p_rgh += pressureCorrection;
        p_rgh.correctBoundaryConditions();
        //
        Info<< "Multi-direct forcing corrector "
            << ibmCorr + 1 << ": max slip = "
            << ibm.maxSlip(U) << nl;
    }
}
#endif
