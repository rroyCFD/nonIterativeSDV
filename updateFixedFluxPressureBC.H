// Update the fixedFluxPressure BCs to ensure flux consistency
setSnGrad<fixedFluxPressureFvPatchScalarField>
(
    p_rgh.boundaryFieldRef(),
    (
        phiHbyA.boundaryField()
      - phiFixedFlux.boundaryField()
    )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
);
p_rgh.correctBoundaryConditions();

forAll(p_rgh.boundaryField(), patchi)
{
    if(isType<fixedFluxPressureFvPatchScalarField> (p_rgh.boundaryField()[patchi]))
    {
        fixedFluxPressureFvPatchScalarField& p_rghGradPatch =
            refCast<fixedFluxPressureFvPatchScalarField>(p_rgh.boundaryFieldRef()[patchi]);

        // Info<< mesh.boundary()[patchi].name()  << "\tp_rgh"   << tab
        //     << "minGrad: " << gMin(p_rghGradPatch.gradient()) << tab
        //     << "maxGrad: " << gMax(p_rghGradPatch.gradient()) << endl;

        // assign the same gradient to pp
        fixedGradientFvPatchScalarField& ppGradPatch        =
            refCast<fixedGradientFvPatchScalarField>(pp.boundaryFieldRef()[patchi]);

        ppGradPatch.gradient() = p_rghGradPatch.gradient();
        Info<< mesh.boundary()[patchi].name() << "\tpp"    << tab
            << "minGrad: " << gMin(ppGradPatch.gradient()) << tab
            << "maxGrad: " << gMax(ppGradPatch.gradient()) << endl;

    }
}