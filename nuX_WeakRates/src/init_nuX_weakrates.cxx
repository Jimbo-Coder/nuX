#include <AMReX.H>

#include <cctk.h>
#include <cctk_Arguments.h>

#include "nuX_weakrates.hxx"

namespace nuX_WeakRates {

WeakRates *global_weakrates = nullptr;

extern "C" void nuX_WeakRates_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_WeakRates_Setup;

  global_weakrates =
      static_cast<WeakRates *>(amrex::The_Managed_Arena()->alloc(
          sizeof(*global_weakrates)));
  new (global_weakrates) WeakRates;
  global_weakrates->init();
}

} // namespace nuX_WeakRates
