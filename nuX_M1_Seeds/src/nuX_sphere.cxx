#include <cmath>
#include <cassert>
#include <loop_device.hxx>
#include <mat.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_seed_utils.hxx"
#include "setup_eos.hxx"
#include "aster_utils.hxx"
#include "nuX_volume.hxx"

namespace nuX_M1_Seeds {

using namespace Loop;
using namespace EOSX;
using namespace AsterUtils;
using namespace nuX_M1_Seeds_volume;

// -----------------------------------------------------------------------------
// Main setup routine
// -----------------------------------------------------------------------------

extern "C" void nuX_M1_Seeds_SetupHydroTest_sphere(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_Seeds_SetupHydroTest_sphere;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_M1_Seeds_SetupHydroTest_sphere");

  auto eos_3p_ig = global_eos_3p_ig;
  if (not CCTK_EQUALS(evolution_eos, "IdealGas")) {
    CCTK_VERROR("Invalid evolution EOS type '%s'. Please, set "
                "EOSX::evolution_eos = \"IdealGas\" in your parameter file.",
                evolution_eos);
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  grid.loop_all_device<1, 1, 1>(grid.nghostzones, [=] CCTK_DEVICE(
                                                      const PointDesc &p) {
    const int ijk = layout_cc.linear(p.i, p.j, p.k);
    rho[ijk] = volume_f(1.0, p.x, p.y, p.z, p.dx, p.dy, p.dz);
    eps[ijk] = static_eps;
    velx[ijk] = 0.0;
    vely[ijk] = 0.0;
    velz[ijk] = 0.0;
    Ye[ijk] = static_ye;
    press[ijk] =
        eos_3p_ig->press_from_rho_eps_ye(rho[ijk], eps[ijk], Ye[ijk]);
  });

}

extern "C" void nuX_M1_Seeds_SetupNeutTest_sphere(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_Seeds_SetupNeutTest_sphere;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_M1_Seeds_SetupNeutTest_sphere");

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  int const ncomponents = radiation_components();
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        for (int ig = 0; ig < ncomponents; ++ig) {
          int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rE[i4D] = rN[i4D] = rFx[i4D] = rFy[i4D] = rFz[i4D] = 0.0;
        }
      });
}

} // namespace nuX_M1_Seeds
