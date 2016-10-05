//  bool LTS = fv::localEulerDdt::enabled(mesh);
    bool LTS = fv::localEulerDdt::enabled(fluidRegions[0]);

tmp<volScalarField> trDeltaT;

if (LTS)
{
    Info<< "Using LTS" << endl;

    trDeltaT = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                fv::localEulerDdt::rDeltaTName,
                runTime.timeName(),
                fluidRegions[0],
//                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
              fluidRegions[0],
//            mesh,
            dimensionedScalar("one", dimless/dimTime, 1),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}
