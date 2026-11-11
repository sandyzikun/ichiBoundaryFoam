{
    // This is an experimental alternating constraint-projection
    // prototype.  It is not yet a monolithic Schur-complement solve.
    ibmAwareProjection.auditOperators();
    ibm.beginLoadCycle();
    if (correctPhi)
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }
    //
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
    volVectorField HbyA
    (
        constrainHbyA(rAU()*UEqn.H(), U, p_rgh)
    );
    ibmAwareProjection.constrainVelocity
    (
        HbyA,
        rho,
        runTime.deltaTValue()
    );
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter
        (
            fvc::interpolate(rho*rAU())
           *fvc::ddtCorr(U, phi, Uf)
        )
    );
    MRF.makeRelative(phiHbyA);
    //
    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }
    //
    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );
    phiHbyA += phig;
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);
    //
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh)
         == fvc::div(phiHbyA)
        );
        p_rghEqn.setReference
        (
            pRefCell,
            getRefCellValue(p_rgh, pRefCell)
        );
        p_rghEqn.solve(p_rgh.select(pimple.finalInnerIter()));
        //
        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();
            p_rgh.relax();
            U = HbyA + rAU()*fvc::reconstruct
            (
                (phig - p_rghEqn.flux())/rAUf
            );
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }
    //
    volScalarField pressureCorrection
    (
        IOobject
        (
            "ibmAwarePressureCorrection",
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
    bool projectionConverged = false;
    scalar constraintResidual = ibm.maxSlip(U);
    scalar continuityResidual =
        ibmAwareProjection.continuityResidual(phi);
    for
    (
        label projectionIter = 1;
        projectionIter <= ibmAwareProjection.maximumIterations();
        ++projectionIter
    )
    {
        volVectorField constrainedVelocity
        (
            IOobject
            (
                "ibmAwareConstrainedVelocity",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            U
        );
        ibmAwareProjection.constrainVelocity
        (
            constrainedVelocity,
            rho,
            runTime.deltaTValue()
        );
        //
        surfaceScalarField constrainedFlux
        (
            "ibmAwareConstrainedFlux",
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
        pressureCorrectionEqn.solve("p_rghFinal");
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
        constraintResidual =
            ibmAwareProjection.constraintResidual(U);
        continuityResidual =
            ibmAwareProjection.continuityResidual(phi);
        ibmAwareProjection.reportIteration
        (
            projectionIter,
            constraintResidual,
            continuityResidual
        );
        projectionConverged = ibmAwareProjection.converged
        (
            projectionIter,
            constraintResidual,
            continuityResidual
        );
        if (projectionConverged)
        {
            break;
        }
    }
    ibmAwareProjection.finalise
    (
        projectionConverged,
        constraintResidual,
        continuityResidual
    );
    //
    #include "continuityErrs.H"
    fvc::correctUf(Uf, U, phi);
    fvc::makeRelative(phi, U);
    //
    p == p_rgh + rho*gh;
    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
    //
    if (!correctPhi)
    {
        rAU.clear();
    }
}
