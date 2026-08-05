#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_LeakageBase_rates.hxx"
#include "nuX_fakerates.hxx"

namespace nuX_LeakageBase {

using namespace Loop;

extern "C" void nuX_LeakageBase_CalcOpacityFakeRates(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_LeakageBase_CalcOpacityFakeRates;
  DECLARE_CCTK_PARAMETERS;

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_LeakageBase_CalcOpacityFakeRates");
  }

  const auto fakerates = nuX_FakeRates::global_fakerates;
  if (!fakerates) {
    CCTK_ERROR("nuX_LeakageBase_CalcOpacityFakeRates requires nuX_FakeRates");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        const auto coeffs = fakerates->ComputeFakeOpacities(rho[ijk]);

        kappa_0_nue[ijk] = coeffs.kappa_0_a[0] + coeffs.kappa_s[0];
        kappa_0_nua[ijk] = coeffs.kappa_0_a[1] + coeffs.kappa_s[1];
        kappa_0_nux[ijk] = coeffs.kappa_0_a[2] + coeffs.kappa_s[2];
        kappa_1_nue[ijk] = coeffs.kappa_a[0] + coeffs.kappa_s[0];
        kappa_1_nua[ijk] = coeffs.kappa_a[1] + coeffs.kappa_s[1];
        kappa_1_nux[ijk] = coeffs.kappa_a[2] + coeffs.kappa_s[2];
        abs_0_nue[ijk] = coeffs.kappa_0_a[0];
        abs_0_nua[ijk] = coeffs.kappa_0_a[1];
        abs_0_nux[ijk] = coeffs.kappa_0_a[2];
      });
}

extern "C" void nuX_LeakageBase_RatesFakeRates(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_LeakageBase_RatesFakeRates;
  DECLARE_CCTK_PARAMETERS;

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_LeakageBase_RatesFakeRates");
  }

  const auto fakerates = nuX_FakeRates::global_fakerates;
  if (!fakerates) {
    CCTK_ERROR("nuX_LeakageBase_RatesFakeRates requires nuX_FakeRates");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        const auto coeffs = fakerates->ComputeFakeOpacities(rho[ijk]);

        CCTK_REAL number[3], energy[3];
        fakerates->FakeNeutrinoDens(rho[ijk], number[0], number[1], number[2],
                                    energy[0], energy[1], energy[2]);

        const CCTK_REAL r_free[3] = {
            coeffs.eta_0[0], coeffs.eta_0[1], coeffs.eta_0[2]};
        const CCTK_REAL q_free[3] = {
            coeffs.eta[0], coeffs.eta[1], coeffs.eta[2]};
        const CCTK_REAL kappa_0[3] = {
            kappa_0_nue[ijk], kappa_0_nua[ijk], kappa_0_nux[ijk]};
        const CCTK_REAL kappa_1[3] = {
            kappa_1_nue[ijk], kappa_1_nua[ijk], kappa_1_nux[ijk]};
        const CCTK_REAL tau_0[3] = {
            optd_0_nue[ijk], optd_0_nua[ijk], optd_0_nux[ijk]};
        const CCTK_REAL tau_1[3] = {
            optd_1_nue[ijk], optd_1_nua[ijk], optd_1_nux[ijk]};

        R_free_nue[ijk] = store_free_rates ? r_free[0] : 0.0;
        R_free_nua[ijk] = store_free_rates ? r_free[1] : 0.0;
        R_free_nux[ijk] = store_free_rates ? r_free[2] : 0.0;
        Q_free_nue[ijk] = store_free_rates ? q_free[0] : 0.0;
        Q_free_nua[ijk] = store_free_rates ? q_free[1] : 0.0;
        Q_free_nux[ijk] = store_free_rates ? q_free[2] : 0.0;

        R_eff_nue[ijk] =
            calc_eff_rate(r_free[0], number[0], kappa_0[0], tau_0[0], DiffFact);
        R_eff_nua[ijk] =
            calc_eff_rate(r_free[1], number[1], kappa_0[1], tau_0[1], DiffFact);
        R_eff_nux[ijk] =
            calc_eff_rate(r_free[2], number[2], kappa_0[2], tau_0[2], DiffFact);
        Q_eff_nue[ijk] =
            calc_eff_rate(q_free[0], energy[0], kappa_1[0], tau_1[0], DiffFact);
        Q_eff_nua[ijk] =
            calc_eff_rate(q_free[1], energy[1], kappa_1[1], tau_1[1], DiffFact);
        Q_eff_nux[ijk] =
            calc_eff_rate(q_free[2], energy[2], kappa_1[2], tau_1[2], DiffFact);
      });
}

} // namespace nuX_LeakageBase
