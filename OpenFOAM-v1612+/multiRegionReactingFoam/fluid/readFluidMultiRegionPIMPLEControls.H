    const dictionary& pimple = mesh.solutionDict().subDict("PIMPLE");

    const int nCorr =
        pimple.lookupOrDefault<int>("nCorrectors", 1);

    const int nNonOrthCorr =
        pimple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    const bool momentumPredictor =
        pimple.lookupOrDefault("momentumPredictor", true);
        
    const bool transonic =
        pimple.lookupOrDefault("transonic", false);
        
    const bool consistent =
        pimple.lookupOrDefault("consistent", false);
