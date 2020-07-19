// This part balances global mass flux.  It does it in a temporary field, and then
// applies the correction indirectly by setting the pressure gradient to be used in
// the fixedFluxPressure boundary condition on p_rgh, or directly if the zeroGradient
// boundary condition on p_rgh is used.
surfaceScalarField phiFixedFlux = phiHbyA;

adjustPhi(phiFixedFlux, U, p_rgh);

forAll(p_rgh.boundaryField(), patchi)
{
    if (isA<zeroGradientFvPatchScalarField>(p_rgh.boundaryField()[patchi]))
    {
        phiHbyA.boundaryFieldRef()[patchi] = phiFixedFlux.boundaryField()[patchi];
    }
}
