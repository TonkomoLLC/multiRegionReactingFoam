if (finalIter)
{
    mesh.data().setFinalIteration(true);
}

if (frozenFlow)
{
    #include "EEqn.H"
}
else
{
    p.storePrevIter();
    rho.storePrevIter();

    if (oCorr == 0)
    {
        #include "rhoEqn.H"
    }

    #include "UEqn.H"
    #include "YEqn.H"
    #include "EEqn.H"

    // --- PISO loop
    for (int corr=0; corr<nCorr; corr++)
    {
          if (consistent)
          {
               #include "pcEqn.H"
          }
          else
          {
               #include "pEqn.H"
          }
    }

    turb.correct();

    rho = thermo.rho();
}

if (finalIter)
{
    mesh.data().setFinalIteration(false);
}
