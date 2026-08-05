#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_M0 {

extern "C" void nuX_M0_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_ParamCheck;
  DECLARE_CCTK_PARAMETERS;

  int const ntheta = static_cast<int>(std::round(std::sqrt(nray / 2.0)));
  if (2 * ntheta * ntheta != nray) {
    CCTK_PARAMWARN("nray must equal 2*ntheta*ntheta for an integer ntheta");
  }

  if (bns_sep_threshold > 0.0 && !CCTK_IsThornActive("BNSTrackerGen")) {
    CCTK_PARAMWARN(
        "bns_sep_threshold requires the BNSTrackerGen thorn to be active");
  }

  if (excision) {
    CCTK_PARAMWARN("M0 ray excision is not implemented yet");
  }

  if (use_enedep_opacity) {
    CCTK_PARAMWARN(
        "Energy-dependent M0 absorption opacity is not implemented yet");
  }
}

} // namespace nuX_M0
