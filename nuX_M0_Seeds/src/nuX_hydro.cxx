#include <cmath>

#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "setup_eos.hxx"

namespace nuX_M0_Seeds {

using namespace Loop;
using namespace EOSX;

enum class test_case { diff_limit_gaussian, diff_limit_square };

extern "C" void
nuX_M0_Seeds_SetupHydroTest_adv_velocity_jump(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_Seeds_SetupHydroTest_adv_velocity_jump;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_M0_Seeds_SetupHydroTest_adv_velocity_jump");

  auto eos_3p_ig = global_eos_3p_ig;
  if (not CCTK_EQUALS(evolution_eos, "IdealGas")) {
    CCTK_VERROR("Invalid evolution EOS type '%s'. Please, set "
                "EOSX::evolution_eos = \"IdealGas\" in your parameter file.",
                evolution_eos);
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});

  CCTK_REAL nx = test_nvec[0];
  CCTK_REAL ny = test_nvec[1];
  CCTK_REAL nz = test_nvec[2];
  CCTK_REAL n2 = nx * nx + ny * ny + nz * nz;

  if (n2 > 0) {
    CCTK_REAL nn = sqrt(n2);
    nx /= nn;
    ny /= nn;
    nz /= nn;
  } else {
    nx = 0.0;
    ny = 0.0;
    nz = 1.0;
  }

  grid.loop_all_device<1, 1, 1>(grid.nghostzones, [=] CCTK_DEVICE(
                                                      const PointDesc &p) {
    const int ijk = layout_cc.linear(p.i, p.j, p.k);
    CCTK_REAL const dotp3d = nx * p.x + ny * p.y + nz * p.z;
    if (dotp3d < 0.0) {
      velx[ijk] = static_velx;
      vely[ijk] = static_vely;
      velz[ijk] = static_velz;
    } else {
      velx[ijk] = -static_velx;
      vely[ijk] = -static_vely;
      velz[ijk] = -static_velz;
    }
    rho[ijk] = static_rho;
    eps[ijk] = static_eps;
    Ye[ijk] = static_ye;
    press[ijk] =
        eos_3p_ig->press_from_rho_eps_ye(rho[ijk], eps[ijk], Ye[ijk]);
  });
}

extern "C" void nuX_M0_Seeds_SetupHydroTest_diff_limit_test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_Seeds_SetupHydroTest_diff_limit_test;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_M0_Seeds_SetupHydroTest_diff_limit_test");

  auto eos_3p_ig = global_eos_3p_ig;
  if (not CCTK_EQUALS(evolution_eos, "IdealGas")) {
    CCTK_VERROR("Invalid evolution EOS type '%s'. Please, set "
                "EOSX::evolution_eos = \"IdealGas\" in your parameter file.",
                evolution_eos);
  }

  test_case tc;
  if CCTK_EQUALS (nuX_test_case, "diff_limit_gaussian") {
    tc = test_case::diff_limit_gaussian;
  } else if (CCTK_EQUALS(nuX_test_case, "diff_limit_square")) {
    tc = test_case::diff_limit_square;
  } else {
    CCTK_ERROR("Unknown value for parameter \"nuX_test_case\"");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(grid.nghostzones, [=] CCTK_DEVICE(
                                                      const PointDesc &p) {
    const int ijk = layout_cc.linear(p.i, p.j, p.k);
    velx[ijk] = static_velx;
    vely[ijk] = static_vely;
    velz[ijk] = static_velz;
    rho[ijk] = profile_hydro_density &&
                       tc == test_case::diff_limit_gaussian
                   ? fmax(rho_atmosphere,
                          static_rho * exp(-9.0 * p.z * p.z))
                   : static_rho;
    eps[ijk] = static_eps;
    Ye[ijk] = static_ye;
    press[ijk] =
        eos_3p_ig->press_from_rho_eps_ye(rho[ijk], eps[ijk], Ye[ijk]);
  });
}

} // namespace nuX_M0_Seeds
