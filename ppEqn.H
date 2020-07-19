{
    // The inverse of the momentum equation "A" operator (the matrix diagonal) and
    // its projection from cell centers to cell faces.  Because it is a diagonal matrix,
    // The inverse is just the diagonal matrix of the reciprocals, hence the "r" in the
    // name "rAU".
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    // The momentum equation "H" operator is the off-diagonal times the last known U.
    // In the equations it is always multiplied with inv(A).  This is HbyA.
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    // This is the Boussinesq buoyancy term multipled by the inverse of the A operator.
    Boussinesq.updateBuoyancyTerm();
    surfaceScalarField phig(rAUf * Boussinesq.buoyancyTerm() * mesh.magSf());

    // Project HbyA to cell faces and apply a correction for time stepping scheme.
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + rAUf*fvc::ddtCorr(U, phi)
    );

    // RhieChow flux correction for pressure gradient
    if(RhieChow)
    {
        tmp<surfaceScalarField> gradpDiff_f
        (   (
                (  (fvc::interpolate(fvc::grad(p_rgh)) & ed)
                  - faceGradient.snGrad(p_rgh)
                ) * ed
            )  & mesh.Sf()
        );

        phiHbyA += rAUf * gradpDiff_f();

        gradpDiff_f.clear();
    }

    // adjust mass flux for variable flow boundaries
    #include "adjustMassFlux.H"

    // include buoyancy term in the boundary flux phiHbyA
    phiHbyA += phig;

    if(fixedFluxBC)
    {
        #include "updateFixedFluxPressureBC.H"
    }

    // update phi2
    phi2 = (phiHbyA - rAUf*fvc::snGrad(p_rgh)*mesh.magSf());

    {
        // find face velocity difference between phi and phi2
        // resulting from  RhieChow-corr, ddtCorr and adjustPhi operations
        surfaceScalarField faceVelDiff(
            "faceVelDiff",
            // (fvc::flux(U) - (phiHbyA - rAUf*fvc::snGrad(p_rgh)*mesh.magSf()))/mesh.magSf()
            (fvc::flux(U) - phi2)/mesh.magSf()
        );

        Info << "faceVelDiff  min: " << gMin(faceVelDiff) << tab
             << "max: " << gMax(faceVelDiff) << endl;
    }


    // Non-orthogonal corrector loop.
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix ppEqn
        (
            fvm::laplacian(rAUf, pp) == fvc::div(phi2)
        );

        scalar ppRefValue = pRefValue - getRefCellValue(p_rgh, pRefCell);
        Info << "ppRefValue: " << ppRefValue << endl;

        ppEqn.setReference(pRefCell, ppRefValue);
        ppEqn.solve(mesh.solver(pp.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            // phi2 -= ppEqn.flux();
            phi = phi2 - ppEqn.flux();

            surfaceScalarField& ppEqnFlux = ppEqn.flux().ref();
            Info<< "min(ppEqn.flux()): " << gMin(ppEqnFlux) << tab
                << "max(ppEqn.flux()): " << gMax(ppEqnFlux) << endl;

            // Explicitly relax pressure for momentum corrector
            pp.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            // U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U -= rAU*fvc::grad(pp);
            U.correctBoundaryConditions();
            fvOptions.correct(U);

            p_rgh += pp;
            p_rgh.correctBoundaryConditions();
        }
    }

    #include "adjustPressureLevel.H"
}
