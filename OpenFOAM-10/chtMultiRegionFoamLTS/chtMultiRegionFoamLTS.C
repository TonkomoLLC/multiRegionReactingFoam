/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

Application
    chtMultiRegionFoam

Description
    Solver for steady or transient fluid flow and solid heat conduction, with
    conjugate heat transfer between regions, buoyancy effects, turbulence,
    reactions and radiation modelling.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "compressibleMomentumTransportModels.H"
#include "fluidReactionThermophysicalTransportModel.H"
#include "fluidReactionThermo.H"
#include "combustionModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "compressibleCourantNo.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "coordinateSystem.H"
#include "pimpleMultiRegionControl.H"
#include "pressureReference.H"
#include "hydrostaticInitialisation.H"

// LTS
#include "pimpleControl.H"
#include "multivariateScheme.H"
#include "fvcSmooth.H"
#include "localEulerDdtScheme.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    pimpleMultiRegionControl pimples(fluidRegions, solidRegions);
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createFluidPressureControls.H"
    #include "createTimeControls.H"
    #include "readSolidTimeControls.H"
    #include "compressibleMultiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"
    #include "setInitialMultiRegionDeltaT.H"

    if (!LTS)
    {
        #include "compressibleMultiRegionCourantNo.H"
        #include "setMultiRegionDeltaT.H"
    }

    while (pimples.run(runTime))
    {
        //#include "readTimeControls.H"
        #include "createTimeControls.H"
        #include "readSolidTimeControls.H"

       if (!LTS)
       {
          #include "compressibleMultiRegionCourantNo.H"
          #include "solidRegionDiffusionNo.H"
          #include "setMultiRegionDeltaT.H"
       }


        runTime++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        // Optional number of energy correctors
        const int nEcorr = pimples.dict().lookupOrDefault<int>
        (
            "nEcorrectors",
            1
        );

        // --- PIMPLE loop
        while (pimples.loop())
        {
            List<tmp<fvVectorMatrix>> UEqns(fluidRegions.size());

            for(int Ecorr=0; Ecorr<nEcorr; Ecorr++)
            {
                forAll(solidRegions, i)
                {
                    Info<< "\nSolving for solid region "
                        << solidRegions[i].name() << endl;
                    #include "setRegionSolidFields.H"
                    #include "solveSolid.H"
                }
                forAll(fluidRegions, i)
                {
                    if (LTS)
                    {
                        #include "setRDeltaTRegion.H" //MULTIREGION LTS ADDITION
                    }
                    Info<< "\nSolving for fluid region "
                        << fluidRegions[i].name() << endl;
                    #include "setRegionFluidFields.H"
                    #include "solveFluid.H"
                }
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
