{
    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution - done in 2 parts. Part 1:
//    thermo.rho() -= psi*p;

//OF4 Update
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    tUEqn.clear();    
//OF4 Update End
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(rho)*(fvc::interpolate(HbyA) & mesh.Sf())
    );


    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

//    while (simple.correctNonOrthogonal())
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rho*rAU, p)
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

    // Second part of thermodynamic density update
//    thermo.rho() += psi*p;

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
}
