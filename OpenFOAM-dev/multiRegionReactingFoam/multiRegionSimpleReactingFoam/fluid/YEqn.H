tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi)")
    )
);


{
    reaction.correct();
    Qdot = reaction.Qdot();
    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
    
        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turb.muEff(), Yi)
             ==
                reaction.R(Yi)
                + 
                fvOptions(rho, Yi)
            );

           YiEqn.relax();
           YiEqn.relax(mesh.equationRelaxationFactor("Yi"));

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
        else
        {
            inertIndex = i;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}

