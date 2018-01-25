/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

{
    fvMesh& mesh = fluidRegions[0];
    CombustionModel<psiReactionThermo>& reaction = reactionFluid[0];
    surfaceScalarField& phi = phiFluid[0];
    volScalarField& rho = rhoFluid[0];    
    PtrList<volScalarField>& Y = reactionFluid[0].thermo().composition().Y();
    basicMultiComponentMixture& composition = reactionFluid[0].thermo().composition();	
    psiReactionThermo& thermo = reaction.thermo();
    
    const word inertSpecie(thermo.lookup("inertSpecie"));
    if (!composition.species().found(inertSpecie))
    {
        FatalIOErrorIn(args.executable().c_str(), thermo)
            << "Inert specie " << inertSpecie << " not found in available species "
            << composition.species()
            << exit(FatalIOError);
    }
    const label inertIndex(composition.species()[inertSpecie]);   

    volScalarField& rDeltaT = trDeltaT.ref();

    const dictionary& pimpleDict = mesh.solutionDict().subDict("PIMPLE");

    // Maximum flow Courant number
    scalar maxCo(readScalar(pimpleDict.lookup("maxCo")));

    // Maximum time scale
    scalar maxDeltaT(pimpleDict.lookupOrDefault<scalar>("maxDeltaT", great));

    // Smoothing parameter (0-1) when smoothing iterations > 0
    scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.1)
    );

    // Damping coefficient (1-0)
    scalar rDeltaTDampingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTDampingCoeff", 1.0)
    );

    // Maximum change in cell temperature per iteration
    // (relative to previous value)
    scalar alphaTemp(pimpleDict.lookupOrDefault("alphaTemp", 0.05));

    // Maximum change in cell concentration per iteration
    // (relative to reference value)
    scalar alphaY(pimpleDict.lookupOrDefault("alphaY", 1.0));

    Info<< "Time scales min/max:" << endl;

    // Cache old reciprocal time scale field
    volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    // Flow time scale
    {
        rDeltaT.ref() =
        (
            fvc::surfaceSum(mag(phi))()()
           /((2*maxCo)*mesh.V()*rho())
        );

        // Limit the largest time scale
        rDeltaT.max(1/maxDeltaT);

        Info<< "    Flow        = "
            << 1/gMax(rDeltaT.primitiveField()) << ", "
            << 1/gMin(rDeltaT.primitiveField()) << endl;
    }

    // Heat release rate time scale
    if (alphaTemp < 1)
    {
        volScalarField::Internal rDeltaTT
        (
            mag(reaction.Qdot())/(alphaTemp*rho*reactionFluid[0].thermo().Cp()*reaction.thermo().T())
        );

        Info<< "    Temperature = "
            << 1/(gMax(rDeltaTT.field()) + vSmall) << ", "
            << 1/(gMin(rDeltaTT.field()) + vSmall) << endl;

        rDeltaT.ref() = max(rDeltaT(), rDeltaTT);
    }

    // Reaction rate time scale
    if (alphaY < 1)
    {
        dictionary Yref(pimpleDict.subDict("Yref"));

        volScalarField::Internal rDeltaTY
        (
            IOobject
            (
                "rDeltaTY",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("rDeltaTY", rDeltaT.dimensions(), 0)
        );

        bool foundY = false;

        forAll(Y, i)
        {
            if (i != inertIndex && composition.active(i))
            {
                volScalarField& Yi = Y[i];

                if (Yref.found(Yi.name()))
                {
                    foundY = true;
                    scalar Yrefi = readScalar(Yref.lookup(Yi.name()));

                    rDeltaTY.field() = max
                    (
                        mag
                        (
                            reaction.R(Yi)().source()
                           /((Yrefi*alphaY)*(rho*mesh.V()))
                        ),
                        rDeltaTY
                    );
                }
            }
        }

        if (foundY)
        {
            Info<< "    Composition = "
                << 1/(gMax(rDeltaTY.field()) + vSmall) << ", "
                << 1/(gMin(rDeltaTY.field()) + vSmall) << endl;

            rDeltaT.ref() = max(rDeltaT(), rDeltaTY);
        }
        else
        {
            IOWarningIn(args.executable().c_str(), Yref)
                << "Cannot find any active species in Yref " << Yref
                << endl;
        }
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    // Spatially smooth the time scale field
    if (rDeltaTSmoothingCoeff < 1)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        rDeltaTDampingCoeff < 1
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        rDeltaT = max
        (
            rDeltaT,
            (scalar(1) - rDeltaTDampingCoeff)*rDeltaT0
        );
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    Info<< "    Overall     = "
        << 1/gMax(rDeltaT.primitiveField())
        << ", " << 1/gMin(rDeltaT.primitiveField()) << endl;
}


// ************************************************************************* //
