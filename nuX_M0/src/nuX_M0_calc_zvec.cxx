#include <loop_device.hxx>
#include <mat.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"

#include "aster_utils.hxx"

namespace nuX_M0 {

using namespace Arith;
using namespace AsterUtils;
using namespace Loop;

extern "C" void nuX_M0_CalcZvec(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_CalcZvec;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});
  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{
      GF3D2<const CCTK_REAL>(layout_vc, gxx),
      GF3D2<const CCTK_REAL>(layout_vc, gxy),
      GF3D2<const CCTK_REAL>(layout_vc, gxz),
      GF3D2<const CCTK_REAL>(layout_vc, gyy),
      GF3D2<const CCTK_REAL>(layout_vc, gyz),
      GF3D2<const CCTK_REAL>(layout_vc, gzz)};

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
          return calc_avg_v2c(gf_g(i, j), p);
        });

        vec<CCTK_REAL, 3> v_up;
        v_up(0) = velx[ijk];
        v_up(1) = vely[ijk];
        v_up(2) = velz[ijk];
        const vec<CCTK_REAL, 3> v_low = calc_contraction(g_avg, v_up);
        const CCTK_REAL w_lorentz = calc_wlorentz(v_low, v_up);

        zvec_cartx[ijk] = w_lorentz * v_up(0);
        zvec_carty[ijk] = w_lorentz * v_up(1);
        zvec_cartz[ijk] = w_lorentz * v_up(2);
      });
}

} // namespace nuX_M0
