# nuX_WeakRates

`nuX_WeakRates` is a C++/GPU port of the analytic weak-interaction rates used
by the legacy THC `WeakRates` thorn.  The library accepts thermodynamic and
nuclear-composition data explicitly and does not read an EOS table itself.
Callers can therefore use the same GPU-resident EOS data as the hydrodynamics.

The original Fortran implementation remains the validation reference.  It is
not a runtime dependency of this thorn.
