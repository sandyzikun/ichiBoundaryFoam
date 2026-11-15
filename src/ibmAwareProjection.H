#ifndef ICHIBOUNDARYFOAM_B361474C_D507_11F1_9472_A3B0B9D8CE68_IBM_AWARE_PROJECTION_H
#define ICHIBOUNDARYFOAM_B361474C_D507_11F1_9472_A3B0B9D8CE68_IBM_AWARE_PROJECTION_H
//
#include "fvCFD.H"
#include "directForcingIbm.H"
//
namespace Foam
{
namespace ichiBoundaryFoam
{
class IbmAwareProjection
{
    DirectForcingIbm& ibm_;
    dictionary properties_;
    bool enabled_;
    word algorithm_;
    scalar constraintTolerance_;
    scalar continuityTolerance_;
    label minimumIterations_;
    label maximumIterations_;
    bool requireConvergence_;
    bool writeDiagnostics_;
    //
public:
    // [INFO] Construction of the experimental projection control
    IbmAwareProjection
    (
        const dictionary& properties,
        DirectForcingIbm& ibm
    )
    :
        ibm_(ibm),
        properties_(properties),
        enabled_
        (
            properties_.lookupOrDefault<bool>("enabled", false)
        ),
        algorithm_
        (
            properties_.lookupOrDefault<word>
            (
                "algorithm",
                "constantPinelliSubIteration"
            )
        ),
        constraintTolerance_
        (
            properties_.lookupOrDefault<scalar>
            (
                "constraintTolerance",
                1.0e-6
            )
        ),
        continuityTolerance_
        (
            properties_.lookupOrDefault<scalar>
            (
                "continuityTolerance",
                1.0e-8
            )
        ),
        minimumIterations_
        (
            properties_.lookupOrDefault<label>
            (
                "minimumIterations",
                1
            )
        ),
        maximumIterations_
        (
            properties_.lookupOrDefault<label>
            (
                "maximumIterations",
                8
            )
        ),
        requireConvergence_
        (
            properties_.lookupOrDefault<bool>
            (
                "requireConvergence",
                false
            )
        ),
        writeDiagnostics_
        (
            properties_.lookupOrDefault<bool>
            (
                "writeDiagnostics",
                false
            )
        )
    {
        if (enabled_ && !ibm_.enabled())
        {
            FatalIOErrorInFunction(properties_)
                << "IBM-aware projection requires enabled IBM"
                << exit(FatalIOError);
        }
        if
        (
            algorithm_ != "constantPinelliSubIteration"
         && algorithm_ != "alternatingConstraintProjection"
        )
        {
            FatalIOErrorInFunction(properties_)
                << "Unsupported IBM-aware algorithm: "
                << algorithm_ << nl
                << "Supported algorithm: constantPinelliSubIteration"
                << exit(FatalIOError);
        }
        if (enabled_ && ibm_.couplingMethod() != "markerBased")
        {
            FatalIOErrorInFunction(properties_)
                << "Constant-Pinelli sub-iteration currently requires "
                << "ibm.couplingMethod markerBased"
                << exit(FatalIOError);
        }
        if (constraintTolerance_ <= 0.0)
        {
            FatalIOErrorInFunction(properties_)
                << "constraintTolerance must be positive"
                << exit(FatalIOError);
        }
        if (continuityTolerance_ <= 0.0)
        {
            FatalIOErrorInFunction(properties_)
                << "continuityTolerance must be positive"
                << exit(FatalIOError);
        }
        if (minimumIterations_ < 1)
        {
            FatalIOErrorInFunction(properties_)
                << "minimumIterations must be at least one"
                << exit(FatalIOError);
        }
        if (maximumIterations_ < minimumIterations_)
        {
            FatalIOErrorInFunction(properties_)
                << "maximumIterations must be no smaller than "
                << "minimumIterations"
                << exit(FatalIOError);
        }
        Info<< "IBM-aware pressure projection: "
            << (enabled_ ? "requested" : "disabled") << nl;
        if (enabled_)
        {
            Info<< "IBM-aware algorithm: "
                << algorithm_ << nl
                << "Constraint/continuity tolerances: "
                << constraintTolerance_ << '/'
                << continuityTolerance_ << nl
                << "Projection iteration range: "
                << minimumIterations_ << '-'
                << maximumIterations_ << nl
                << "Require convergence: "
                << (requireConvergence_ ? "yes" : "no") << nl;
        }
    }
    //
    // [INFO] Application of one immersed velocity correction
    void constrainVelocity
    (
        volVectorField& velocity,
        const volScalarField& density,
        const scalar deltaT
    )
    {
        ibm_.constrainPredictor(velocity, density, deltaT);
    }
    //
    // [INFO] Evaluation of the immersed constraint residual
    scalar constraintResidual
    (
        const volVectorField& velocity
    ) const
    {
        return ibm_.maxSlip(velocity);
    }
    //
    // [INFO] Evaluation of the finite-volume continuity residual
    scalar continuityResidual
    (
        const surfaceScalarField& flux
    ) const
    {
        const tmp<volScalarField> tResidual
        (
            mag(fvc::div(flux))
        );
        return gMax(tResidual().primitiveField());
    }
    //
    // [INFO] Testing of the coupled residual stopping condition
    bool converged
    (
        const label iteration,
        const scalar constraintResidualValue,
        const scalar continuityResidualValue
    ) const
    {
        return
            iteration >= minimumIterations_
         && constraintResidualValue <= constraintTolerance_
         && continuityResidualValue <= continuityTolerance_;
    }
    //
    // [INFO] Reporting of an IBM-aware projection iteration
    void reportIteration
    (
        const label iteration,
        const scalar constraintResidualValue,
        const scalar continuityResidualValue
    ) const
    {
        if (writeDiagnostics_)
        {
            Info<< "IBM-aware projection iteration "
                << iteration << ": constraint residual = "
                << constraintResidualValue
                << ", continuity residual = "
                << continuityResidualValue << nl;
        }
    }
    //
    // [INFO] Finalisation of an IBM-aware projection cycle
    void finalise
    (
        const bool hasConverged,
        const scalar constraintResidualValue,
        const scalar continuityResidualValue
    ) const
    {
        if (hasConverged)
        {
            return;
        }
        if (requireConvergence_)
        {
            FatalErrorInFunction
                << "IBM-aware projection failed to converge: "
                << "constraint residual = "
                << constraintResidualValue
                << ", continuity residual = "
                << continuityResidualValue
                << exit(FatalError);
        }
        WarningInFunction
            << "IBM-aware projection reached its iteration limit: "
            << "constraint residual = "
            << constraintResidualValue
            << ", continuity residual = "
            << continuityResidualValue << nl;
    }
    //
    // [INFO] Access to the run-time activation state
    bool enabled() const
    {
        return enabled_;
    }
    //
    // [INFO] Access to the maximum correction count
    label maximumIterations() const
    {
        return maximumIterations_;
    }
};
}
}
//
#endif
