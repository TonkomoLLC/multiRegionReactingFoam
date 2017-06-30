rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

if (nCorr <= 1)
{
    tUEqn.clear();
}

if (transonic)
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            fvc::flux(HbyA)
          + rhorAUf*fvc::ddtCorr(rho, U, phi)/fvc::interpolate(rho)
        )
    );

    MRF.makeRelative(fvc::interpolate(psi), phid);

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - fvm::laplacian(rhorAUf, p)
         ==
            fvOptions(psi, p, rho.name())
        );

//        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
          pEqn.solve
            (
                mesh.solver
                (
                    p.select
                    (
                        (
                           oCorr == nOuterCorr-1
                        && corr == nCorr-1
                        && nonOrth == nNonOrthCorr
                        )
                    )
                )
            );


        if (nonOrth == nNonOrthCorr)
        {
            phi == pEqn.flux();
        }
    }
}
else
{
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::flux(rho*HbyA)
          + rhorAUf*fvc::ddtCorr(rho, U, phi)
        )
    );

    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phiHbyA)
          - fvm::laplacian(rhorAUf, p)
         ==
            fvOptions(psi, p, rho.name())
        );

//        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

          pEqn.solve
            (
                mesh.solver
                (
                    p.select
                    (
                        (
                           oCorr == nOuterCorr-1
                        && corr == nCorr-1
                        && nonOrth == nNonOrthCorr
                        )
                    )
                )
            );


        if (nonOrth == nNonOrthCorr)
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}

#include "rhoEqn.H"
//#include "compressibleContinuityErrs.H"
#include "compressibleContinuityErrors.H"


// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

if (pressureControl.limit(p))
{
    p.correctBoundaryConditions();
    rho = thermo.rho();
}

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}
