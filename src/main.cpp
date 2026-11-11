#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#ifndef ICHIBOUNDARYFOAM_ENABLE_MDF_PRESSURE_PROJECTION
#define ICHIBOUNDARYFOAM_ENABLE_MDF_PRESSURE_PROJECTION 0
#endif
#ifndef ICHIBOUNDARYFOAM_ENABLE_IBM_AWARE_PROJECTION
#define ICHIBOUNDARYFOAM_ENABLE_IBM_AWARE_PROJECTION 0
#endif
#ifndef ICHIBOUNDARYFOAM_ENABLE_ONE_WAY_DPM
#define ICHIBOUNDARYFOAM_ENABLE_ONE_WAY_DPM 0
#endif
#ifndef ICHIBOUNDARYFOAM_ENABLE_TURBULENCE_WALL_TREATMENT
#define ICHIBOUNDARYFOAM_ENABLE_TURBULENCE_WALL_TREATMENT 1
#endif
#if ICHIBOUNDARYFOAM_ENABLE_ONE_WAY_DPM
#include "basicKinematicCloud.H"
#endif
#if !ICHIBOUNDARYFOAM_PURE_INTERFOAM
#include "immersedGeometry.H"
#include "directForcingIbm.H"
#if ICHIBOUNDARYFOAM_ENABLE_TURBULENCE_WALL_TREATMENT
#include "immersedTurbulenceWallTreatment.H"
#endif
#if ICHIBOUNDARYFOAM_ENABLE_IBM_AWARE_PROJECTION
#include "ibmAwareProjection.H"
#endif
//
using namespace Foam;
#endif
//
int main(int argc, char* argv[])
{
#if ICHIBOUNDARYFOAM_ENABLE_ONE_WAY_DPM
    argList::addOption
    (
        "cloud",
        "name",
        "specify the passive cloud name; default is kinematicCloud"
    );
#endif
    argList::addNote
    (
        "[INFO] "
#if ICHIBOUNDARYFOAM_PURE_INTERFOAM
        "Pure interFoam regression path"
#else
        "Two-phase VOF solver with SDF-based direct-forcing IBM"
#endif
    );
    //
    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"
    //
#if ICHIBOUNDARYFOAM_ENABLE_ONE_WAY_DPM
    const word oneWayDpmCloudName
    (
        args.getOrDefault<word>("cloud", "kinematicCloud")
    );
    volScalarField oneWayDpmMu
    (
        IOobject
        (
            "oneWayDpmMu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mixture.mu()
    );
    Info<< "Constructing passive one-way cloud "
        << oneWayDpmCloudName << endl;
    basicKinematicCloud oneWayDpmCloud
    (
        oneWayDpmCloudName,
        rho,
        U,
        oneWayDpmMu,
        g
    );
    Info<< "[WARN] The passive cloud has no SDF collision model; "
        << "particles can cross immersed geometry." << endl;
    //
#endif
#if !ICHIBOUNDARYFOAM_PURE_INTERFOAM
    Foam::ichiBoundaryFoam::ImmersedGeometry geometry(mesh);
    IOdictionary ichiBoundaryProperties
    (
        IOobject
        (
            "ichiBoundaryGeometry",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    Foam::ichiBoundaryFoam::DirectForcingIbm ibm
    (
        mesh,
        geometry,
        ichiBoundaryProperties.subOrEmptyDict("ibm")
    );
#if ICHIBOUNDARYFOAM_ENABLE_TURBULENCE_WALL_TREATMENT
    Foam::ichiBoundaryFoam::ImmersedTurbulenceWallTreatment
        immersedTurbulenceWallTreatment
    (
        mesh,
        geometry,
        ibm,
        ichiBoundaryProperties.subOrEmptyDict
        (
            "immersedTurbulenceWallTreatment"
        )
    );
#endif
#if ICHIBOUNDARYFOAM_ENABLE_IBM_AWARE_PROJECTION
    Foam::ichiBoundaryFoam::IbmAwareProjection ibmAwareProjection
    (
        ichiBoundaryProperties.subOrEmptyDict
        (
            "ibmAwareProjection"
        ),
        ibm
    );
#endif
    geometry.update(runTime.value());
    ibm.classify();
    geometry.write();
    //
#endif
    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }
    //
    Info<< "\nStarting time loop\n" << endl;
    //
    while (runTime.run())
    {
        #include "readDyMControls.H"
        //
        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }
        //
        ++runTime;
        Info<< "Time = " << runTime.timeName() << nl << endl;
#if ICHIBOUNDARYFOAM_ENABLE_ONE_WAY_DPM
        oneWayDpmCloud.storeGlobalPositions();
#endif
#if !ICHIBOUNDARYFOAM_PURE_INTERFOAM
        ibm.beginTimeStep();
#endif
        //
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();
                //
                if (mesh.changing())
                {
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }
                    //
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;
                    MRF.update();
                    //
                    if (correctPhi)
                    {
                        phi = mesh.Sf() & Uf();
                        #include "correctPhi.H"
                        fvc::makeRelative(phi, U);
                        mixture.correct();
                    }
                    //
                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }
            //
#if !ICHIBOUNDARYFOAM_PURE_INTERFOAM
            geometry.update(runTime.value());
            ibm.classify();
#endif
            //
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"
#if !ICHIBOUNDARYFOAM_PURE_INTERFOAM
            ibm.extendAlpha(alpha1);
#endif
            mixture.correct();
            //
            if (pimple.frozenFlow())
            {
                continue;
            }
            //
#if ICHIBOUNDARYFOAM_PURE_INTERFOAM
            #include "UEqn.H"
#else
            MRF.correctBoundaryVelocity(U);
            fvVectorMatrix UEqn
            (
                fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
              + MRF.DDt(rho, U)
              + turbulence->divDevRhoReff(rho, U)
             ==
                fvOptions(rho, U)
            );
            UEqn.relax();
            fvOptions.constrain(UEqn);
            if (pimple.momentumPredictor())
            {
                solve
                (
                    UEqn
                 ==
                    fvc::reconstruct
                    (
                        (
                            mixture.surfaceTensionForce()
                          - ghf*fvc::snGrad(rho)
                          - fvc::snGrad(p_rgh)
                        ) * mesh.magSf()
                    )
                );
                fvOptions.correct(U);
            }
#endif
            //
            while (pimple.correct())
            {
#if ICHIBOUNDARYFOAM_PURE_INTERFOAM
                #include <pEqn.H>
#else
#if ICHIBOUNDARYFOAM_ENABLE_IBM_AWARE_PROJECTION
                if (ibmAwareProjection.enabled())
                {
                    #include "ibmAwarePEqn.H"
                }
                else
                {
                    #include "pEqn.H"
                }
#else
                #include "pEqn.H"
#endif
#endif
            }
#if !ICHIBOUNDARYFOAM_PURE_INTERFOAM
            ibm.reportResiduals(U, phi, alpha1);
            ibm.reportLoads(runTime.deltaTValue());
            ibm.reportPhaseVolume(alpha1);
#endif
            //
            if (pimple.turbCorr())
            {
                turbulence->correct();
#if !ICHIBOUNDARYFOAM_PURE_INTERFOAM \
 && ICHIBOUNDARYFOAM_ENABLE_TURBULENCE_WALL_TREATMENT
                const tmp<volScalarField> tNu(mixture.nu());
                immersedTurbulenceWallTreatment.constrain
                (
                    alpha1,
                    tNu()
                );
#endif
            }
        }
        //
#if ICHIBOUNDARYFOAM_ENABLE_ONE_WAY_DPM
        oneWayDpmMu = mixture.mu();
        oneWayDpmMu.correctBoundaryConditions();
        Info<< "Evolving passive one-way cloud "
            << oneWayDpmCloud.name() << endl;
        oneWayDpmCloud.evolve();
#endif
        //
        runTime.write();
        runTime.printExecutionTime(Info);
    }
    //
    Info<< "End\n" << endl;
    return 0;
}
