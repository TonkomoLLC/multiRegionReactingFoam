# multiRegionReactingFoam
OpenFOAM transient solver for laminar or turbulent fluid flow and solid heat conduction with conjugate heat transfer between solid and fluid regions, plus combustion with chemical reactions (psi thermo model)

Notes
1. OpenFOAM-v1706 solver also compiles with OpenFOAM-v1712
2. OpenFOAM-v1806 solver also compiles with OpenFOAM-v1906
3. Release includes transient and steady state solvers (multiRegionReactingFoam and multiRegionSimpleReactingFoam)
4. Unsteady state solver benchmarked against DETCHEM, FLUENT and reactingFoam.
5. The main difference with chtMultiRegionReactingFoam is the form of PEqn.H.  PEqn in chtMultiRegionReactingFoam
   is based on PEqn from chtMultiRegionFoam, while PEqn in multiRegionReactingFoam is based on reactingFoam. One may find
   multiRegionReactingFoam to be better than chtMultiRegionReactingFoam for some situations, especially where there
   are large pressure changes during the simulation.
6. LTS is based on a one fluid region (region 0). Therefore, if there are multiple fluid regions, LTS time steps will
   be chosen based on the characteristic times of the fluid assigned to region 0. Therefore, in a multiple fluid region
   application where LTS will be used, one may wish to assign the most "time step sensitive" fluid region to fluid region
   0 (e.g., by placing the name of this more sensitive fluid region first in under "fluid" in constant/regionProperties.
7. The OpenFOAM-6 solver does not include a steady state solver (similar to the chtMultiRegionFoam solver in OpenFOAM-6).

