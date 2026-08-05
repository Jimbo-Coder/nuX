#include <loop_device.hxx>

#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_LeakageBase {

using namespace Loop;

extern "C" void nuX_LeakageBase_InitOpticalDepthSimple(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_LeakageBase_InitOpticalDepthSimple;
  DECLARE_CCTK_PARAMETERS;
  using std::isfinite;

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_LeakageBase_InitOpticalDepthSimple");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        CCTK_REAL rho_local = rho[ijk];
        if (!isfinite(rho_local) || rho_local < atmo_rho) {
          rho_local = atmo_rho;
        }
        CCTK_REAL const tau = rho_local / atmo_rho - 1.0;
        optd_0_nue[ijk] = tau;
        optd_0_nua[ijk] = tau;
        optd_0_nux[ijk] = tau;
        optd_1_nue[ijk] = tau;
        optd_1_nua[ijk] = tau;
        optd_1_nux[ijk] = tau;
        optd_0_nue_p[ijk] = tau;
        optd_0_nua_p[ijk] = tau;
        optd_0_nux_p[ijk] = tau;
        optd_1_nue_p[ijk] = tau;
        optd_1_nua_p[ijk] = tau;
        optd_1_nux_p[ijk] = tau;
        optd_0_nue_p_p[ijk] = tau;
        optd_0_nua_p_p[ijk] = tau;
        optd_0_nux_p_p[ijk] = tau;
        optd_1_nue_p_p[ijk] = tau;
        optd_1_nua_p_p[ijk] = tau;
        optd_1_nux_p_p[ijk] = tau;
      });
}

} // namespace nuX_LeakageBase
