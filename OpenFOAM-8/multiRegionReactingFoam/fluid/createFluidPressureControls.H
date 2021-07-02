PtrList<pressureControl> pressureControlFluid(fluidRegions.size());

forAll(fluidRegions, i)
{
    pressureControlFluid.set
    (
        i,
        new pressureControl
        (
            reactionFluid[i].thermo().p(),
            rhoFluid[i],
            pimples.pimple(i).dict(),
            false
        )
    );
}
