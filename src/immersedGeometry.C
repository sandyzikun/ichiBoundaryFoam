#include "immersedGeometry.H"
//
namespace Foam
{
namespace ichiBoundaryFoam
{
scalar ImmersedGeometry::localSignedDistance
(
    const point& localPoint
) const
{
    const scalar radial =
        sqrt(sqr(localPoint.x()) + sqr(localPoint.y()));
    const scalar radialDistance = radial - radius_;
    const scalar axialDistance = mag(localPoint.z()) - halfLength_;
    const scalar outsideRadial = max(radialDistance, scalar(0));
    const scalar outsideAxial = max(axialDistance, scalar(0));
    const scalar outsideDistance =
        sqrt(sqr(outsideRadial) + sqr(outsideAxial));
    const scalar insideDistance =
        min(max(radialDistance, axialDistance), scalar(0));
    return outsideDistance + insideDistance;
}
//
vector ImmersedGeometry::localNormal(const point& localPoint) const
{
    const scalar radial =
        sqrt(sqr(localPoint.x()) + sqr(localPoint.y()));
    const scalar radialDistance = radial - radius_;
    const scalar axialDistance = mag(localPoint.z()) - halfLength_;
    vector radialNormal(1.0, 0.0, 0.0);
    if (radial > VSMALL)
    {
        radialNormal = vector
        (
            localPoint.x()/radial,
            localPoint.y()/radial,
            0.0
        );
    }
    const vector axialNormal
    (
        0.0,
        0.0,
        localPoint.z() >= 0.0 ? 1.0 : -1.0
    );
    if (radialDistance > 0.0 && axialDistance > 0.0)
    {
        return normalised
        (
            radialDistance*radialNormal
          + axialDistance*axialNormal
        );
    }
    return radialDistance >= axialDistance
        ? radialNormal : axialNormal;
}
//
void ImmersedGeometry::evaluate
(
    const point& globalPoint,
    scalar& distance,
    vector& normal,
    vector& velocity
) const
{
    const point localPoint =
        rotation_.T() & (globalPoint - centre_);
    distance = localSignedDistance(localPoint);
    normal = rotation_ & localNormal(localPoint);
    velocity =
        linearVelocity_
      + (angularVelocity_ ^ (globalPoint - centre_));
}
//
ImmersedGeometry::ImmersedGeometry(const fvMesh& mesh)
:
    mesh_(mesh),
    properties_
    (
        IOobject
        (
            "ichiBoundaryGeometry",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    geometryType_(properties_.subDict("geometry").get<word>("type")),
    referenceCentre_
    (
        properties_.subDict("geometry").get<point>("centre")
    ),
    referenceRotation_
    (
        properties_.subDict("geometry").lookupOrDefault<vector>
        (
            "referenceRotation",
            vector::zero
        )
    ),
    radius_(properties_.subDict("geometry").get<scalar>("radius")),
    halfLength_
    (
        0.5*properties_.subDict("geometry").get<scalar>("length")
    ),
    motion_
    (
        properties_.subDict("motion"),
        referenceCentre_,
        referenceRotation_
    ),
    centre_(referenceCentre_),
    rotation_(tensor::I),
    linearVelocity_(vector::zero),
    angularVelocity_(vector::zero),
    signedDistance_
    (
        IOobject
        (
            properties_.lookupOrDefault<word>
            (
                "fieldName",
                "signedDistance"
            ),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimLength, Zero)
    ),
    surfaceNormal_
    (
        IOobject
        (
            "immersedSurfaceNormal",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimless, Zero)
    ),
    surfaceVelocity_
    (
        IOobject
        (
            "immersedSurfaceVelocity",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimVelocity, Zero)
    )
{
    if (geometryType_ != "cappedCylinder")
    {
        FatalIOErrorInFunction(properties_)
            << "Unsupported geometry type: " << geometryType_ << nl
            << "Supported geometry types: cappedCylinder"
            << exit(FatalIOError);
    }
    if (radius_ <= 0.0 || halfLength_ <= 0.0)
    {
        FatalIOErrorInFunction(properties_)
            << "Cylinder radius and length must be positive"
            << exit(FatalIOError);
    }
}
//
void ImmersedGeometry::update(const scalar timeValue)
{
    centre_ = motion_.centre(timeValue);
    const vector bodyAngles = motion_.angles(timeValue);
    rotation_ = RigidBodyMotion::rotationTensor(bodyAngles);
    linearVelocity_ = motion_.linearVelocity(timeValue);
    angularVelocity_ = motion_.angularVelocity(timeValue);
    scalarField& distance = signedDistance_.primitiveFieldRef();
    vectorField& normal = surfaceNormal_.primitiveFieldRef();
    vectorField& velocity = surfaceVelocity_.primitiveFieldRef();
    forAll(distance, celli)
    {
        evaluate
        (
            mesh_.C()[celli],
            distance[celli],
            normal[celli],
            velocity[celli]
        );
    }
    forAll(signedDistance_.boundaryFieldRef(), patchi)
    {
        fvPatchScalarField& patchDistance =
            signedDistance_.boundaryFieldRef()[patchi];
        fvPatchVectorField& patchNormal =
            surfaceNormal_.boundaryFieldRef()[patchi];
        fvPatchVectorField& patchVelocity =
            surfaceVelocity_.boundaryFieldRef()[patchi];
        const vectorField& faceCentres =
            mesh_.Cf().boundaryField()[patchi];
        forAll(patchDistance, facei)
        {
            scalar valueDistance = 0.0;
            vector valueNormal(vector::zero);
            vector valueVelocity(vector::zero);
            evaluate
            (
                faceCentres[facei],
                valueDistance,
                valueNormal,
                valueVelocity
            );
            patchDistance[facei] = valueDistance;
            patchNormal[facei] = valueNormal;
            patchVelocity[facei] = valueVelocity;
        }
    }
    Info<< "SDF geometry: centre = " << centre_
        << ", rotation = " << bodyAngles
        << ", min/max = " << gMin(distance)
        << '/' << gMax(distance) << nl;
}
//
const volScalarField& ImmersedGeometry::signedDistance() const
{
    return signedDistance_;
}
//
const volVectorField& ImmersedGeometry::surfaceNormal() const
{
    return surfaceNormal_;
}
//
const volVectorField& ImmersedGeometry::surfaceVelocity() const
{
    return surfaceVelocity_;
}
//
const point& ImmersedGeometry::centre() const
{
    return centre_;
}
//
bool ImmersedGeometry::write() const
{
    return
        signedDistance_.write()
     && surfaceNormal_.write()
     && surfaceVelocity_.write();
}
}
}
