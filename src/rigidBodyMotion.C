#include "rigidBodyMotion.H"
//
namespace Foam
{
namespace ichiBoundaryFoam
{
scalar RigidBodyMotion::harmonicValue
(
    const scalar reference,
    const scalar amplitude,
    const scalar frequency,
    const scalar phase,
    const scalar timeValue
) const
{
    const scalar argument =
        constant::mathematical::twoPi
       *frequency*(timeValue - referenceTime_) + phase;
    return reference + amplitude*sin(argument);
}
//
scalar RigidBodyMotion::harmonicDerivative
(
    const scalar amplitude,
    const scalar frequency,
    const scalar phase,
    const scalar timeValue
) const
{
    const scalar omega = constant::mathematical::twoPi*frequency;
    const scalar argument =
        omega*(timeValue - referenceTime_) + phase;
    return amplitude*omega*cos(argument);
}
//
RigidBodyMotion::RigidBodyMotion
(
    const dictionary& properties,
    const point& referenceCentre,
    const vector& referenceRotation
)
:
    properties_(properties),
    type_(properties_.lookupOrDefault<word>("type", "static")),
    referenceCentre_(referenceCentre),
    referenceRotation_(referenceRotation),
    translationAmplitude_(vector::zero),
    translationFrequency_(vector::zero),
    translationPhase_(vector::zero),
    rotationAmplitude_(vector::zero),
    rotationFrequency_(vector::zero),
    rotationPhase_(vector::zero),
    translationVelocity_(vector::zero),
    rotationRate_(vector::zero),
    referenceTime_
    (
        properties_.lookupOrDefault<scalar>("referenceTime", 0.0)
    )
{
    if (type_ == "constant6DoF")
    {
        translationVelocity_ = properties_.lookupOrDefault<vector>
        (
            "translationVelocity",
            vector::zero
        );
        rotationRate_ = properties_.lookupOrDefault<vector>
        (
            "rotationRate",
            vector::zero
        );
    }
    else if (type_ == "harmonic6DoF")
    {
        translationAmplitude_ =
            properties_.get<vector>("translationAmplitude");
        translationFrequency_ =
            properties_.get<vector>("translationFrequency");
        translationPhase_ = properties_.lookupOrDefault<vector>
        (
            "translationPhase",
            vector::zero
        );
        rotationAmplitude_ =
            properties_.get<vector>("rotationAmplitude");
        rotationFrequency_ =
            properties_.get<vector>("rotationFrequency");
        rotationPhase_ = properties_.lookupOrDefault<vector>
        (
            "rotationPhase",
            vector::zero
        );
    }
    else if (type_ != "static")
    {
        FatalIOErrorInFunction(properties_)
            << "Unsupported motion type: " << type_ << nl
            << "Supported motion types: static, constant6DoF, "
            << "harmonic6DoF"
            << exit(FatalIOError);
    }
}
//
tensor RigidBodyMotion::rotationTensor(const vector& angles)
{
    const scalar cx = cos(angles.x());
    const scalar sx = sin(angles.x());
    const scalar cy = cos(angles.y());
    const scalar sy = sin(angles.y());
    const scalar cz = cos(angles.z());
    const scalar sz = sin(angles.z());
    const tensor rotationX
    (
        1.0, 0.0, 0.0,
        0.0, cx, -sx,
        0.0, sx, cx
    );
    const tensor rotationY
    (
        cy, 0.0, sy,
        0.0, 1.0, 0.0,
        -sy, 0.0, cy
    );
    const tensor rotationZ
    (
        cz, -sz, 0.0,
        sz, cz, 0.0,
        0.0, 0.0, 1.0
    );
    return rotationZ & rotationY & rotationX;
}
//
point RigidBodyMotion::centre(const scalar timeValue) const
{
    point value(referenceCentre_);
    if (type_ == "static")
    {
        return value;
    }
    if (type_ == "constant6DoF")
    {
        return value
             + translationVelocity_*(timeValue - referenceTime_);
    }
    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        value[cmpt] = harmonicValue
        (
            referenceCentre_[cmpt],
            translationAmplitude_[cmpt],
            translationFrequency_[cmpt],
            translationPhase_[cmpt],
            timeValue
        );
    }
    return value;
}
//
vector RigidBodyMotion::angles(const scalar timeValue) const
{
    vector value(referenceRotation_);
    if (type_ == "static")
    {
        return value;
    }
    if (type_ == "constant6DoF")
    {
        return value + rotationRate_*(timeValue - referenceTime_);
    }
    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        value[cmpt] = harmonicValue
        (
            referenceRotation_[cmpt],
            rotationAmplitude_[cmpt],
            rotationFrequency_[cmpt],
            rotationPhase_[cmpt],
            timeValue
        );
    }
    return value;
}
//
vector RigidBodyMotion::linearVelocity(const scalar timeValue) const
{
    vector value(vector::zero);
    if (type_ == "static")
    {
        return value;
    }
    if (type_ == "constant6DoF")
    {
        return translationVelocity_;
    }
    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        value[cmpt] = harmonicDerivative
        (
            translationAmplitude_[cmpt],
            translationFrequency_[cmpt],
            translationPhase_[cmpt],
            timeValue
        );
    }
    return value;
}
//
vector RigidBodyMotion::angularVelocity(const scalar timeValue) const
{
    if (type_ == "static")
    {
        return vector::zero;
    }
    vector rates(vector::zero);
    if (type_ == "constant6DoF")
    {
        rates = rotationRate_;
    }
    else
    {
        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            rates[cmpt] = harmonicDerivative
            (
                rotationAmplitude_[cmpt],
                rotationFrequency_[cmpt],
                rotationPhase_[cmpt],
                timeValue
            );
        }
    }
    const vector valueAngles = angles(timeValue);
    const scalar roll = valueAngles.x();
    const scalar pitch = valueAngles.y();
    const scalar yawRate = rates.z();
    const scalar pitchRate = rates.y();
    const scalar rollRate = rates.x();
    return vector
    (
        rollRate - yawRate*sin(pitch),
        pitchRate*cos(roll) + yawRate*sin(roll)*cos(pitch),
        -pitchRate*sin(roll) + yawRate*cos(roll)*cos(pitch)
    );
}
}
}
