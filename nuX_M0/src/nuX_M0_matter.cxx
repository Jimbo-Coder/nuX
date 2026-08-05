#include <loop_device.hxx>

#include <cassert>
#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "aster_utils.hxx"
#include "nuX_baryon_mass.hxx"

namespace nuX_M0 {

using namespace Arith;
using namespace AsterUtils;
using namespace Loop;

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
apply_semi_implicit_source(const CCTK_REAL old, const CCTK_REAL source) {
  if (old == 0.0)
    return source;
  const CCTK_REAL denominator = 1.0 - source / old;
  if (denominator == 0.0)
    return old + source;
  return old / denominator;
}

extern "C" void nuX_M0_CalcMatterRHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_nuX_M0_CalcMatterRHS;
  DECLARE_CCTK_PARAMETERS;

  if (verbose && CCTK_MyProc(cctkGH) == 0)
    CCTK_INFO("nuX_M0_CalcMatterRHS");

  const CCTK_REAL mb = nuX_Utils::AverageBaryonMass(particle_mass);
  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{
      gxx, gxy, gxz, gyy, gyz, gzz};

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        if (rho(p.I) <= density_cutoff)
          return;

        const CCTK_REAL alp_cc = calc_avg_v2c(alp, p);
        const smat<CCTK_REAL, 3> g([&](int i, int j) ARITH_INLINE {
          return calc_avg_v2c(gf_g(i, j), p);
        });
        const CCTK_REAL detg = calc_det(g);
        assert(detg > 0.0);
        const CCTK_REAL sqrt_detg = sqrt(detg);

        const vec<CCTK_REAL, 3> v_up{
            velx(p.I), vely(p.I), velz(p.I)};
        const vec<CCTK_REAL, 3> v_low = calc_contraction(g, v_up);
        const CCTK_REAL w_lorentz = calc_wlorentz(v_low, v_up);

        const CCTK_REAL R =
            abs_number(p.I) - mb * (R_eff_nue(p.I) - R_eff_nua(p.I));
        const CCTK_REAL Q =
            abs_energy(p.I) -
            (Q_eff_nue(p.I) + Q_eff_nua(p.I) + Q_eff_nux(p.I));
        const CCTK_REAL factor = alp_cc * sqrt_detg;

        DYe_rhs(p.I) += factor * R;
        momxrhs(p.I) += factor * w_lorentz * v_low(0) * Q;
        momyrhs(p.I) += factor * w_lorentz * v_low(1) * Q;
        momzrhs(p.I) += factor * w_lorentz * v_low(2) * Q;
        taurhs(p.I) += factor * w_lorentz * Q;

        assert(isfinite(DYe_rhs(p.I)));
        assert(isfinite(momxrhs(p.I)));
        assert(isfinite(momyrhs(p.I)));
        assert(isfinite(momzrhs(p.I)));
        assert(isfinite(taurhs(p.I)));
      });
}

extern "C" void nuX_M0_ApplyMatterSource(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_nuX_M0_ApplyMatterSource;
  DECLARE_CCTK_PARAMETERS;

  if (verbose && CCTK_MyProc(cctkGH) == 0)
    CCTK_INFO("nuX_M0_ApplyMatterSource");

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const CCTK_REAL mb = nuX_Utils::AverageBaryonMass(particle_mass);
  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{
      gxx, gxy, gxz, gyy, gyz, gzz};

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        if (rho(p.I) <= density_cutoff)
          return;

        const CCTK_REAL alp_cc = calc_avg_v2c(alp, p);
        const smat<CCTK_REAL, 3> g([&](int i, int j) ARITH_INLINE {
          return calc_avg_v2c(gf_g(i, j), p);
        });
        const CCTK_REAL detg = calc_det(g);
        assert(detg > 0.0);
        const CCTK_REAL sqrt_detg = sqrt(detg);
        const vec<CCTK_REAL, 3> v_up{
            velx(p.I), vely(p.I), velz(p.I)};
        const vec<CCTK_REAL, 3> v_low = calc_contraction(g, v_up);
        const CCTK_REAL w_lorentz = calc_wlorentz(v_low, v_up);

        const CCTK_REAL R =
            abs_number(p.I) - mb * (R_eff_nue(p.I) - R_eff_nua(p.I));
        const CCTK_REAL Q =
            abs_energy(p.I) -
            (Q_eff_nue(p.I) + Q_eff_nua(p.I) + Q_eff_nux(p.I));
        const CCTK_REAL factor = alp_cc * sqrt_detg;
        const CCTK_REAL DYe_dot = factor * R;
        const CCTK_REAL tau_dot = factor * w_lorentz * Q;

        // Freeze each source-to-state ratio over the diagonal IMEX stage and
        // solve that linearized source implicitly. For a zero state or an
        // exactly singular linearization, fall back to its additive limit.
        DYe(p.I) =
            apply_semi_implicit_source(DYe(p.I), dt * DYe_dot);
        momx(p.I) += dt * factor * w_lorentz * v_low(0) * Q;
        momy(p.I) += dt * factor * w_lorentz * v_low(1) * Q;
        momz(p.I) += dt * factor * w_lorentz * v_low(2) * Q;
        tau(p.I) =
            apply_semi_implicit_source(tau(p.I), dt * tau_dot);

        assert(isfinite(DYe(p.I)));
        assert(isfinite(momx(p.I)));
        assert(isfinite(momy(p.I)));
        assert(isfinite(momz(p.I)));
        assert(isfinite(tau(p.I)));
      });
}

extern "C" void nuX_M0_CalcLuminosity(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_nuX_M0_CalcLuminosity;

  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{
      gxx, gxy, gxz, gyy, gyz, gzz};
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const CCTK_REAL alp_cc = calc_avg_v2c(alp, p);
        const smat<CCTK_REAL, 3> g([&](int i, int j) ARITH_INLINE {
          return calc_avg_v2c(gf_g(i, j), p);
        });
        const CCTK_REAL detg = calc_det(g);
        assert(detg > 0.0);
        const CCTK_REAL sqrt_detg = sqrt(detg);
        const vec<CCTK_REAL, 3> v_up{
            velx(p.I), vely(p.I), velz(p.I)};
        const CCTK_REAL w_lorentz =
            calc_wlorentz(calc_contraction(g, v_up), v_up);
        const CCTK_REAL factor =
            alp_cc * alp_cc * sqrt_detg * w_lorentz;

        luminosity_nue(p.I) = factor * Q_eff_nue(p.I);
        luminosity_nua(p.I) = factor * Q_eff_nua(p.I);
        luminosity_nux(p.I) = factor * Q_eff_nux(p.I);
      });
}

} // namespace nuX_M0
