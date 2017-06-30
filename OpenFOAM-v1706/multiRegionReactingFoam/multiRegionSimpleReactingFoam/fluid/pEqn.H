// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
tUEqn.clear();
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::interpolate(rho)*fvc::flux(HbyA)
);

MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);

//    while (simple.correctNonOrthogonal())
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pEqn
    (
        fvc::div(phiHbyA)
      - fvm::laplacian(rhorAUf, p)
     ==
        fvOptions(psi, p, rho.name())
    );

    pEqn.solve();

//        if (simple.finalNonOrthogonalIter())
        if (nonOrth == nNonOrthCorr)
    {
        phi = phiHbyA + pEqn.flux();
    }
}

p.relax();

// Thermodynamic density update
thermo.correctRho(psi*p - psip0);

// //    #include "compressibleContinuityErrs.H"
    #include "continuityErrs.H"

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

rho = thermo.rho();
rho = max(rho, rhoMin[i]);
rho = min(rho, rhoMax[i]);
rho.relax();

Info<< "p min/max = " << min(p).value() << ", " << max(p).value() << endl;
