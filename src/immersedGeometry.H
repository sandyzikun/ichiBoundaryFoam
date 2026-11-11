#ifndef ICHIBOUNDARYFOAM_B361474C_D507_11F1_9472_A3B0B9D8CE68_IMMERSEDGEOMETRY_H
#define ICHIBOUNDARYFOAM_B361474C_D507_11F1_9472_A3B0B9D8CE68_IMMERSEDGEOMETRY_H
//
#include "fvCFD.H"
#include "rigidBodyMotion.H"
//
namespace Foam
{
namespace ichiBoundaryFoam
{
class ImmersedGeometry
{
    const fvMesh& mesh_;
    IOdictionary properties_;
    word geometryType_;
    point referenceCentre_;
    vector referenceRotation_;
    scalar radius_;
    scalar halfLength_;
    RigidBodyMotion motion_;
    point centre_;
    tensor rotation_;
    vector linearVelocity_;
    vector angularVelocity_;
    volScalarField signedDistance_;
    volVectorField surfaceNormal_;
    volVectorField surfaceVelocity_;
    //
    // [INFO] Evaluation of the capped-cylinder signed distance
    scalar localSignedDistance(const point& localPoint) const;
    //
    // [INFO] Evaluation of the capped-cylinder SDF normal
    vector localNormal(const point& localPoint) const;
    //
    // [INFO] Evaluation of geometry fields at one point
    void evaluate
    (
        const point& globalPoint,
        scalar& distance,
        vector& normal,
        vector& velocity
    ) const;
    //
public:
    // [INFO] Construction of an analytic immersed geometry
    explicit ImmersedGeometry(const fvMesh& mesh);
    //
    // [INFO] Update of pose and geometry fields
    void update(const scalar timeValue);
    //
    // [INFO] Access to the signed-distance field
    const volScalarField& signedDistance() const;
    //
    // [INFO] Access to the outward surface-normal field
    const volVectorField& surfaceNormal() const;
    //
    // [INFO] Access to the rigid-body surface-velocity field
    const volVectorField& surfaceVelocity() const;
    //
    // [INFO] Access to the current body centre
    const point& centre() const;
    //
    // [INFO] Writing of geometry fields
    bool write() const;
};
}
}
//
#endif
