// Initialise fluid field pointer lists
PtrList<combustionModels::psiCombustionModel> reactionFluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
PtrList<compressible::turbulenceModel> turbulence(fluidRegions.size());
PtrList<radiation::radiationModel> radiation(fluidRegions.size());
PtrList<volScalarField> QdotFluid(fluidRegions.size());
List<scalar> initialMassFluid(fluidRegions.size());
List<bool> frozenFlowFluid(fluidRegions.size(), false);

PtrList<dimensionedScalar> rhoMax(fluidRegions.size());
PtrList<dimensionedScalar> rhoMin(fluidRegions.size());

//OF4 Update
PtrList<multivariateSurfaceInterpolationScheme<scalar>::fieldTable> fieldsFluid(fluidRegions.size());
PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    Info<< "Creating reaction model\n" << endl;

    reactionFluid.set
    (
        i,
        autoPtr<combustionModels::psiCombustionModel>(combustionModels::psiCombustionModel::New(fluidRegions[i])).ptr()
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
            reactionFluid[i].thermo().rho()
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


    Info<< "    Adding to turbulence\n" << endl;
    turbulence.set
    (
        i,
        compressible::turbulenceModel::New
        (
            rhoFluid[i],
            UFluid[i],
            phiFluid[i],
            reactionFluid[i].thermo()
        ).ptr()
    );

   reactionFluid[i].setTurbulence(turbulence[i]);

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

    radiation.set
    (
        i,
        radiation::radiationModel::New(reactionFluid[i].thermo().T())
    );

    initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

    const dictionary& simpleDict =
        fluidRegions[i].solutionDict().subDict("SIMPLE");

    simpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);

    rhoMax.set
    (
        i,
        new dimensionedScalar
        (
            dimensionedScalar::lookupOrDefault
            (
                "rhoMax",
                simpleDict,
                dimDensity,
                GREAT
            )
        )
    );

    rhoMin.set
    (
        i,
        new dimensionedScalar
        (
            dimensionedScalar::lookupOrDefault
            (
                "rhoMin",
                simpleDict,
                dimDensity,
                0
            )
        )
    );

    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    Info<< "    Adding fvOptions\n" << endl;
//OF4 Update
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

    basicMultiComponentMixture& composition = reactionFluid[i].thermo().composition();	
		
	
	forAll(composition.Y(), j)
	{
		    fieldsFluid[i].add(composition.Y()[i]);
	}


		
    fieldsFluid[i].add(reactionFluid[i].thermo().he());
    
    turbulence[i].validate();
}
