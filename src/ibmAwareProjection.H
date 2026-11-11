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
    word operatorModel_;
    scalar constraintTolerance_;
    scalar continuityTolerance_;
    scalar adjointTolerance_;
    label minimumIterations_;
    label maximumIterations_;
    bool requireConvergence_;
    bool writeDiagnostics_;
    mutable bool operatorsAudited_;
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
                "alternatingConstraintProjection"
            )
        ),
        operatorModel_
        (
            properties_.lookupOrDefault<word>
            (
                "operatorModel",
                "cellVolumeWeightedDiffuse"
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
        adjointTolerance_
        (
            properties_.lookupOrDefault<scalar>
            (
                "adjointTolerance",
                1.0e-12
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
        ),
        operatorsAudited_(false)
    {
        if (enabled_ && !ibm_.enabled())
        {
            FatalIOErrorInFunction(properties_)
                << "IBM-aware projection requires enabled IBM"
                << exit(FatalIOError);
        }
        if (algorithm_ != "alternatingConstraintProjection")
        {
            FatalIOErrorInFunction(properties_)
                << "Unsupported IBM-aware algorithm: "
                << algorithm_ << nl
                << "Supported prototype algorithm: "
                << "alternatingConstraintProjection"
                << exit(FatalIOError);
        }
        if (operatorModel_ != "cellVolumeWeightedDiffuse")
        {
            FatalIOErrorInFunction(properties_)
                << "Unsupported IBM-aware operator model: "
                << operatorModel_ << nl
                << "Supported operator model: "
                << "cellVolumeWeightedDiffuse"
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
        if (adjointTolerance_ <= 0.0)
        {
            FatalIOErrorInFunction(properties_)
                << "adjointTolerance must be positive"
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
            Info<< "IBM-aware prototype algorithm: "
                << algorithm_ << nl
                << "IBM-aware constraint operator: "
                << operatorModel_ << nl
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
    // [INFO] Application of the linear immersed constraint operator
    void applyConstraint
    (
        const vectorField& velocityIncrement,
        vectorField& constraintValue
    ) const
    {
        const scalarField& weight =
            ibm_.solidFraction().primitiveField();
        if (velocityIncrement.size() != weight.size())
        {
            FatalErrorInFunction
                << "Constraint input size differs from the mesh size"
                << exit(FatalError);
        }
        constraintValue.setSize(weight.size());
        forAll(weight, celli)
        {
            constraintValue[celli] =
                weight[celli]*velocityIncrement[celli];
        }
    }
    //
    // [INFO] Application of the adjoint force-spreading operator
    void spreadMultiplier
    (
        const vectorField& multiplier,
        vectorField& spreadValue
    ) const
    {
        const scalarField& weight =
            ibm_.solidFraction().primitiveField();
        if (multiplier.size() != weight.size())
        {
            FatalErrorInFunction
                << "Multiplier input size differs from the mesh size"
                << exit(FatalError);
        }
        spreadValue.setSize(weight.size());
        forAll(weight, celli)
        {
            spreadValue[celli] = weight[celli]*multiplier[celli];
        }
    }
    //
    // [INFO] Auditing of the discrete constraint-spreading pair
    void auditOperators() const
    {
        if (!enabled_ || operatorsAudited_)
        {
            return;
        }
        const fvMesh& mesh = ibm_.solidFraction().mesh();
        const vectorField& centre = mesh.C().primitiveField();
        const scalarField& volume = mesh.V();
        vectorField testVelocity(mesh.nCells(), vector::zero);
        vectorField testMultiplier(mesh.nCells(), vector::zero);
        forAll(testVelocity, celli)
        {
            const vector& c = centre[celli];
            testVelocity[celli] = vector
            (
                1.0 + 0.17*c.x(),
                -0.31 + 0.11*c.y(),
                0.23 - 0.07*c.z()
            );
            testMultiplier[celli] = vector
            (
                -0.19 + 0.13*c.z(),
                0.29 - 0.05*c.x(),
                0.37 + 0.09*c.y()
            );
        }
        vectorField constrainedVelocity;
        vectorField spreadMultiplierValue;
        applyConstraint(testVelocity, constrainedVelocity);
        spreadMultiplier(testMultiplier, spreadMultiplierValue);
        scalar localConstraintWork = 0.0;
        scalar localSpreadWork = 0.0;
        forAll(volume, celli)
        {
            localConstraintWork += volume[celli]
               *(constrainedVelocity[celli] & testMultiplier[celli]);
            localSpreadWork += volume[celli]
               *(testVelocity[celli]
               & spreadMultiplierValue[celli]);
        }
        const scalar constraintWork = returnReduce
        (
            localConstraintWork,
            sumOp<scalar>()
        );
        const scalar spreadWork = returnReduce
        (
            localSpreadWork,
            sumOp<scalar>()
        );
        const scalar relativeError = mag
        (
            constraintWork - spreadWork
        )/max
        (
            max(mag(constraintWork), mag(spreadWork)),
            VSMALL
        );
        Info<< "IBM-aware operator audit: <Jv,lambda> = "
            << constraintWork << ", <v,S(lambda)> = "
            << spreadWork << ", relative error = "
            << relativeError << nl;
        if (relativeError > adjointTolerance_)
        {
            FatalErrorInFunction
                << "IBM constraint and spreading operators fail the "
                << "discrete-adjoint audit: relative error = "
                << relativeError << ", tolerance = "
                << adjointTolerance_
                << exit(FatalError);
        }
        operatorsAudited_ = true;
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
