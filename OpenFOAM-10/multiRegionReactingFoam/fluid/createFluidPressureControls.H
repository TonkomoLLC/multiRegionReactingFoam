PtrList<pressureReference> pressureReferenceFluid(fluidRegions.size());

forAll(fluidRegions, i)
{
    pressureReferenceFluid.set
    (
        i,
        new pressureReference
        (
            thermoFluid[i].p(),
            pimples.pimple(i).dict(),
            false
        )
    );
}
/*
    hydrostaticInitialisation
    (
        pFluid[i],
        thermoFluid[i].p(),
        rhoFluid[i],
        UFluid[i],
        ghFluid[i],
        ghfFluid[i],
        pRefFluid[i],
        thermoFluid[i],
        pimples.pimple(i).dict()
    );

}
*/

