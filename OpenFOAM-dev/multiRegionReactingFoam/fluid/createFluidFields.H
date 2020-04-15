// Initialise fluid field pointer lists
PtrList<psiReactionThermo> thermoFluid(fluidRegions.size());
PtrList<CombustionModel<psiReactionThermo>> reactionFluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
PtrList<volScalarField> CpFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<compressible::momentumTransportModel>
    turbulenceFluid(fluidRegions.size());
PtrList<compressible::thermophysicalTransportModel>
    thermophysicalTransportFluid(fluidRegions.size());
PtrList<radiationModel> radiation(fluidRegions.size());
PtrList<volScalarField> KFluid(fluidRegions.size());
PtrList<volScalarField> dpdtFluid(fluidRegions.size());
PtrList<volScalarField> QdotFluid(fluidRegions.size());

List<scalar> initialMassFluid(fluidRegions.size());
List<bool> frozenFlowFluid(fluidRegions.size(), false);
PtrList<multivariateSurfaceInterpolationScheme<scalar>::fieldTable> fieldsFluid(fluidRegions.size());

PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());


#include "createRDeltaT.H"

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{

    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;
        
        
    Info<< "Reading thermophysical properties\n" << endl;

    thermoFluid.set
    (
        i,
          autoPtr<psiReactionThermo>(psiReactionThermo::New(fluidRegions[i])).ptr()
    );

    Info<< "    Adding to rhoFluid\n" << endl;
    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermoFluid[i].rho()

        )
    );

    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(rhoFluid[i]*UFluid[i])
                & fluidRegions[i].Sf()
        )
    );

    Info<< "    Adding to gFluid\n" << endl;
    gFluid.set
    (
        i,
        new uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                runTime.constant(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    CpFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Cp",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermoFluid[i].Cp()

        )
    );
    
    Info<< "    Adding to turbulenceFluid\n" << endl;
    turbulenceFluid.set
    (
        i,
        compressible::momentumTransportModel::New
        (
            rhoFluid[i],
            UFluid[i],
            phiFluid[i],
            thermoFluid[i]
        ).ptr()
    );

    Info<< "    Adding to thermophysicalTransport\n" << endl;
    thermophysicalTransportFluid.set
    (
        i,
        compressible::thermophysicalTransportModel::New
        (
            turbulenceFluid[i]
        ).ptr()
    );
    
    Info<< "Creating reaction model\n" << endl;
 
    reactionFluid.set
    (
        i,
          autoPtr<CombustionModel<psiReactionThermo>>
          (
          CombustionModel<psiReactionThermo>::New
           (
            thermoFluid[i],
            turbulenceFluid[i]          
           )
          )
    );


      
    radiation.set
    (
        i,
        radiationModel::New(thermoFluid[i].T())
    );

    initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

    Info<< "    Adding to KFluid\n" << endl;
    KFluid.set
    (
        i,
        new volScalarField
        (
            "K",
            0.5*magSqr(UFluid[i])
        )
    );

    Info<< "    Adding to dpdtFluid\n" << endl;
    dpdtFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "dpdt",
                runTime.timeName(),
                fluidRegions[i]
            ),
            fluidRegions[i],
            dimensionedScalar
            (
                "dpdt",
                reactionFluid[i].thermo().p().dimensions()/dimTime,
                0
            )
        )
    );

   Info<< "    Adding to Qdot\n" << endl;
   QdotFluid.set
   (
        i,
        new volScalarField
        (
            IOobject
            (
                "Qdot",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
                fluidRegions[i],
                dimensionedScalar
            (
                "Qdot",
                dimEnergy/dimVolume/dimTime,
                0
            )
        )
    );


    const dictionary& pimpleDict =
        fluidRegions[i].solutionDict().subDict("PIMPLE");
    pimpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);


    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    Info<< "    Adding fvOptions\n" << endl;

    fluidFvOptions.set
     (
        i,
        new fv::options(fluidRegions[i])
     );
 
 
    fieldsFluid.set
    	(
    	    i,
			new multivariateSurfaceInterpolationScheme<scalar>::fieldTable 
		);

    forAll(thermoFluid[i].composition().Y(), j)
    {
        fieldsFluid[i].add(thermoFluid[i].composition().Y()[j]);
    }

    fieldsFluid[i].add(thermoFluid[i].he());

		
    fieldsFluid[i].add(reactionFluid[i].thermo().he());

    turbulenceFluid[i].validate();    
}


