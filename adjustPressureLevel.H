// Compute the full pressure and adjust the pressure level.
Boussinesq.updateBackgroundPressure();

// Check whether the background pressure has changed
Info << "backgroundPressure at refCell: "
    << getRefCellValue(Boussinesq.backgroundPressure(), pRefCell) << endl;

p = p_rgh + Boussinesq.backgroundPressure();

if (p_rghNeedRef)
{
    Info<< "p_rgh[pRefCell] (before adjustment): "
    	<< getRefCellValue(p_rgh, pRefCell) << endl;

    dimensionedScalar pAdjust(
    	"pAdjust",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell));

    p     += pAdjust;
    p_rgh += pAdjust;

    Info<< "p_rgh[pRefCell] (after  adjustment): "
    	<< getRefCellValue(p_rgh, pRefCell) << endl;
}
