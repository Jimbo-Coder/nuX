#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_LeakageBase {

using namespace Loop;

extern "C" void nuX_LeakageBase_InitData(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_LeakageBase_InitData;
  DECLARE_CCTK_PARAMETERS;

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_LeakageBase_InitData");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        kappa_0_nue[ijk] = 0.0;
        kappa_0_nua[ijk] = 0.0;
        kappa_0_nux[ijk] = 0.0;
        kappa_1_nue[ijk] = 0.0;
        kappa_1_nua[ijk] = 0.0;
        kappa_1_nux[ijk] = 0.0;
        abs_0_nue[ijk] = 0.0;
        abs_0_nua[ijk] = 0.0;
        abs_0_nux[ijk] = 0.0;

        R_free_nue[ijk] = 0.0;
        R_free_nua[ijk] = 0.0;
        R_free_nux[ijk] = 0.0;
        Q_free_nue[ijk] = 0.0;
        Q_free_nua[ijk] = 0.0;
        Q_free_nux[ijk] = 0.0;

        R_eff_nue[ijk] = 0.0;
        R_eff_nua[ijk] = 0.0;
        R_eff_nux[ijk] = 0.0;
        Q_eff_nue[ijk] = 0.0;
        Q_eff_nua[ijk] = 0.0;
        Q_eff_nux[ijk] = 0.0;

        abs_number[ijk] = 0.0;
        abs_energy[ijk] = 0.0;
        luminosity_nue[ijk] = 0.0;
        luminosity_nua[ijk] = 0.0;
        luminosity_nux[ijk] = 0.0;

        optd_0_nue[ijk] = 0.0;
        optd_0_nua[ijk] = 0.0;
        optd_0_nux[ijk] = 0.0;
        optd_1_nue[ijk] = 0.0;
        optd_1_nua[ijk] = 0.0;
        optd_1_nux[ijk] = 0.0;

        optd_0_nue_p[ijk] = 0.0;
        optd_0_nua_p[ijk] = 0.0;
        optd_0_nux_p[ijk] = 0.0;
        optd_1_nue_p[ijk] = 0.0;
        optd_1_nua_p[ijk] = 0.0;
        optd_1_nux_p[ijk] = 0.0;

        optd_0_nue_p_p[ijk] = 0.0;
        optd_0_nua_p_p[ijk] = 0.0;
        optd_0_nux_p_p[ijk] = 0.0;
        optd_1_nue_p_p[ijk] = 0.0;
        optd_1_nua_p_p[ijk] = 0.0;
        optd_1_nux_p_p[ijk] = 0.0;
      });
}

} // namespace nuX_LeakageBase
