if (ibm.multiDirectForcing())
{
    for
    (
        label ibmCorr = 0;
        ibmCorr < ibm.nMultiDirectForcingCorrectors();
        ++ibmCorr
    )
    {
        const scalar slipBefore = ibm.maxSlip(HbyA);
        if
        (
            ibmCorr > 0
         && ibm.slipTolerance() > 0.0
         && slipBefore <= ibm.slipTolerance()
        )
        {
            Info<< "Predictor MDF converged after "
                << ibmCorr << " correctors, max slip = "
                << slipBefore << nl;
            break;
        }
        //
        ibm.constrainPredictor
        (
            HbyA,
            rho,
            runTime.deltaTValue()
        );
        //
        Info<< "Predictor MDF corrector "
            << ibmCorr + 1 << ": max slip = "
            << ibm.maxSlip(HbyA) << nl;
    }
}
