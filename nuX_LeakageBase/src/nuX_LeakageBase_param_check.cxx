#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_LeakageBase {

extern "C" void nuX_LeakageBase_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_LeakageBase_ParamCheck;
  DECLARE_CCTK_PARAMETERS;

  if (atmo_rho <= 0.0) {
    CCTK_PARAMWARN("atmo_rho must be positive");
  }

  if (CCTK_Equals(rates_lib, "FakeRates") &&
      !CCTK_IsThornActive("nuX_FakeRates")) {
    CCTK_PARAMWARN("nuX_LeakageBase::rates_lib=FakeRates requires the "
                   "nuX_FakeRates thorn to be active");
  } else if (CCTK_Equals(rates_lib, "WeakRates") &&
             !CCTK_IsThornActive("nuX_WeakRates")) {
    CCTK_PARAMWARN("nuX_LeakageBase::rates_lib=WeakRates requires the "
                   "nuX_WeakRates thorn to be active");
  } else if (CCTK_Equals(rates_lib, "NuRates") &&
             !CCTK_IsThornActive("nuX_NuRates")) {
    CCTK_PARAMWARN("nuX_LeakageBase::rates_lib=NuRates requires the "
                   "nuX_NuRates thorn to be active");
  }

  if (CCTK_Equals(init_method, "spherical")) {
    CCTK_PARAMWARN(
        "Spherical optical-depth initialization is not implemented yet");
  }
}

} // namespace nuX_LeakageBase
