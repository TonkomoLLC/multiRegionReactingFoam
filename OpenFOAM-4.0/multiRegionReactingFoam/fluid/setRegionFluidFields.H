    fvMesh& mesh = fluidRegions[i];
 
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields = fieldsFluid[i];
 
    combustionModels::psiCombustionModel& reaction = reactionFluid[i];
    psiReactionThermo& thermo = reaction.thermo();
    PtrList<volScalarField>& Y = reactionFluid[i].thermo().composition().Y(); 
    thermo.validate(args.executable(), "h", "e");
    
    volScalarField& rho = rhoFluid[i];
    volVectorField& U = UFluid[i];
    const dimensionedVector& g = gFluid[i];
    surfaceScalarField& phi = phiFluid[i];

    compressible::turbulenceModel& turb = turbulence[i];
    
    reaction.setTurbulence(turb);
    
    volScalarField& K = KFluid[i];
    volScalarField& dpdt = dpdtFluid[i];

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();

    volScalarField& dQ = dQFluid[i];

    radiation::radiationModel& rad = radiation[i];

    IOMRFZoneList& MRF = MRFfluid[i];
    fv::options& fvOptions = fluidFvOptions[i];

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );

    const word inertSpecie(thermo.lookup("inertSpecie"));
    const bool frozenFlow = frozenFlowFluid[i];
