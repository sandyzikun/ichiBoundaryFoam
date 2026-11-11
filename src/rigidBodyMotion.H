#ifndef ICHIBOUNDARYFOAM_B361474C_D507_11F1_9472_A3B0B9D8CE68_RIGIDBODYMOTION_H
#define ICHIBOUNDARYFOAM_B361474C_D507_11F1_9472_A3B0B9D8CE68_RIGIDBODYMOTION_H
//
#include "fvCFD.H"
//
namespace Foam
{
namespace ichiBoundaryFoam
{
class RigidBodyMotion
{
    dictionary properties_;
    word type_;
    point referenceCentre_;
    vector referenceRotation_;
    vector translationAmplitude_;
    vector translationFrequency_;
    vector translationPhase_;
    vector rotationAmplitude_;
    vector rotationFrequency_;
    vector rotationPhase_;
    vector translationVelocity_;
    vector rotationRate_;
    scalar referenceTime_;
    //
    // [INFO] Evaluation of one harmonic coordinate
    scalar harmonicValue
    (
        const scalar reference,
        const scalar amplitude,
        const scalar frequency,
        const scalar phase,
        const scalar timeValue
    ) const;
    //
    // [INFO] Evaluation of one harmonic coordinate derivative
    scalar harmonicDerivative
    (
        const scalar amplitude,
        const scalar frequency,
        const scalar phase,
        const scalar timeValue
    ) const;
    //
public:
    // [INFO] Construction of prescribed rigid-body motion
    RigidBodyMotion
    (
        const dictionary& properties,
        const point& referenceCentre,
        const vector& referenceRotation
    );
    //
    // [INFO] Conversion of roll-pitch-yaw angles to a rotation tensor
    static tensor rotationTensor(const vector& angles);
    //
    // [INFO] Evaluation of the rigid-body centre
    point centre(const scalar timeValue) const;
    //
    // [INFO] Evaluation of the roll-pitch-yaw angles
    vector angles(const scalar timeValue) const;
    //
    // [INFO] Evaluation of the translational velocity
    vector linearVelocity(const scalar timeValue) const;
    //
    // [INFO] Evaluation of the spatial angular velocity
    vector angularVelocity(const scalar timeValue) const;
};
}
}
//
#endif
