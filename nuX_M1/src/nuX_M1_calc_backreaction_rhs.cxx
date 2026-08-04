#include <cassert>
#include <cmath>
#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "avg_baryon_mass.hpp"

namespace nuX_M1 {

using namespace std;
using namespace Loop;
extern "C" void nuX_M1_CalcBackreactionRHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcBackreactionRHS;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_M1_CalcBackreactionRHS");

  if (ngroups != 1 || nspecies != 3)
    CCTK_ERROR("nuX_M1::backreact requires ngroups=1 and nspecies=3");

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});

  const CCTK_REAL mb = AverageBaryonMass(particle_mass);

  // Keep the equal-and-opposite matter source in the implicit partition.  In
  // particular, the first (zero-diagonal) stage is used by later rows of both
  // IMEX32L and IMEX42L and must be conservative as well.
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        const int groupspec = ngroups * nspecies;

        for (int ig = 0; ig < groupspec; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          const CCTK_REAL Ndot = rN_rhs[i4D];

          momxrhs[ijk] -= rFx_rhs[i4D];
          momyrhs[ijk] -= rFy_rhs[i4D];
          momzrhs[ijk] -= rFz_rhs[i4D];
          taurhs[ijk] -= rE_rhs[i4D];
          const CCTK_REAL DYe_dot_ig =
              -mb * ((ig == 0 ? Ndot : 0.0) - (ig == 1 ? Ndot : 0.0));
          DYe_rhs[ijk] += DYe_dot_ig;

          assert(isfinite(momxrhs[ijk]));
          assert(isfinite(momyrhs[ijk]));
          assert(isfinite(momzrhs[ijk]));
          assert(isfinite(taurhs[ijk]));
          assert(isfinite(DYe_rhs[ijk]));
        }
      });
}

} // namespace nuX_M1
