#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "nuX_M0_device.hxx"
#include "nuX_M0_schedule.hxx"

namespace nuX_M0 {

extern "C" void nuX_M0_InitData(CCTK_ARGUMENTS);

extern "C" void nuX_M0_Switch(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_Switch;
  DECLARE_CCTK_PARAMETERS;

  if (!CarpetX_AllLevelsSynchronized() ||
      !iteration_is_due(cctk_iteration, compute_every)) {
    return;
  }

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_M0_Switch");
  }

  if (wait_until_time > 0.0 && cctk_time < wait_until_time) {
    return;
  }

  bool const was_on = load_device_scalar(nuX_M0_is_on) != 0;
  bool is_on;

  if (bns_sep_threshold > 0.0) {
    CCTK_REAL const *bns_sep_tot = static_cast<CCTK_REAL const *>(
        CCTK_VarDataPtr(cctkGH, 0, "BNSTrackerGen::bns_sep_tot"));
    is_on = bns_sep_tot &&
            load_device_scalar(bns_sep_tot) < bns_sep_threshold;
  } else {
    is_on = true;
  }
  store_device_scalar(nuX_M0_is_on, CCTK_INT(is_on));

  if (!was_on && is_on) {
    store_device_scalar(nuX_M0_time, CCTK_REAL(cctk_time));
    if (verbose && CCTK_MyProc(cctkGH) == 0) {
      CCTK_INFO("nuX_M0_Switch: M0 is on");
    }
  } else if (was_on && !is_on) {
    nuX_M0_InitData(CCTK_PASS_CTOC);
    if (verbose && CCTK_MyProc(cctkGH) == 0) {
      CCTK_INFO("nuX_M0_Switch: M0 is off");
    }
  }
}

} // namespace nuX_M0
