/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    multiRegionReactingFoam

Description
    Transient solver for laminar or turbulent fluid flow and solid heat
    conduction with conjugate heat transfer between solid and fluid regions,
    plus combustion with chemical reactions.

    It handles secondary fluid or solid circuits which can be coupled
    thermally with the main fluid region. i.e radiators, etc.


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
//#include "psiCombustionModel.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "turbulentFluidThermoModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "pressureControl.H"
#include "compressibleCourantNo.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "coordinateSystem.H"
#include "fvOptions.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "multivariateScheme.H"
#include "fvcSmooth.H"
#include "localEulerDdtScheme.H"
#include "pimpleMultiRegionControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"

    #include "createMeshes.H"
    #include "createControl.H"
    #include "createFields.H"

    pimpleMultiRegionControl pimples(fluidRegions, solidRegions);
    #include "createFluidPressureControls.H"
    #include "initContinuityErrs.H"
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

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    
    while (pimples.run(runTime))
    {
//        #include "createTimeControls.H"
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
 
          {
              #include "compressibleMultiRegionCourantNo.H"
              #include "solidRegionDiffusionNo.H"
              #include "setMultiRegionDeltaT.H"
          }
       
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        // --- PIMPLE loop
        while (pimples.loop())
        {
            forAll(fluidRegions, i)
            {
                Info<< "\nSolving for fluid region "
                    << fluidRegions[i].name() << endl;
                #include "setRegionFluidFields.H"
                #include "solveFluid.H"
            }

            forAll(solidRegions, i)
            {
                Info<< "\nSolving for solid region "
                    << solidRegions[i].name() << endl;
                #include "setRegionSolidFields.H"
                #include "solveSolid.H"
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
