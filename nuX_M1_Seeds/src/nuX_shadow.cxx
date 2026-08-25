#include <cmath>

#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_seed_utils.hxx"
#include "setup_eos.hxx"
#include "nuX_volume.hxx"

namespace nuX_M1_Seeds {

using namespace Loop;
using namespace EOSX;
using namespace nuX_M1_Seeds_volume;

// -----------------------------------------------------------------------------
// Main setup routine
// -----------------------------------------------------------------------------

extern "C" void nuX_M1_Seeds_SetupHydroTest_shadow(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_Seeds_SetupHydroTest_shadow;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_M1_Seeds_SetupHydroTest_shadow");

  auto eos_3p_ig = global_eos_3p_ig;
  if (not CCTK_EQUALS(evolution_eos, "IdealGas")) {
    CCTK_VERROR("Invalid evolution EOS type '%s'. Please, set "
                "EOSX::evolution_eos = \"IdealGas\" in your parameter file.",
                evolution_eos);
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        // Match THC shadow hydro setup: unit-radius density sphere at rest.
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

extern "C" void nuX_M1_Seeds_SetupNeutTest_shadow(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_Seeds_SetupNeutTest_shadow;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_M1_Seeds_SetupNeutTest_shadow");

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  int const ncomponents = radiation_components();

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        for (int ig = 0; ig < ncomponents; ++ig) {
          int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rE[i4D] = rN[i4D] = rFx[i4D] = rFy[i4D] = rFz[i4D] = 0.0;
        }
      });
}

extern "C" void nuX_M1_Seeds_ShadowBCs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_Seeds_ShadowBCs;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_M1_Seeds_ShadowBCs");

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  int const ncomponents = radiation_components();

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
    nx = 1.0;
    ny = 0.0;
    nz = 0.0;
  }

  // Initialize physical ghost zones to vacuum.
  grid.loop_bnd_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        for (int ig = 0; ig < ncomponents; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rE[i4D] = 0.0;
          rFx[i4D] = 0.0;
          rFy[i4D] = 0.0;
          rFz[i4D] = 0.0;
          rN[i4D] = 0.0;
        }
      });

  const auto set_inflow =
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        bool is_inflow_boundary = false;
        CCTK_REAL E_new = 0.0;
        CCTK_REAL Fx_new = 0.0;
        CCTK_REAL Fy_new = 0.0;
        CCTK_REAL Fz_new = 0.0;
        CCTK_REAL N_new = 0.0;

        // Set interior (BI) and ghost (NI) inflow states for reconstruction.
        // The sign selects the face; 0.5 filters small transverse components.
        if (nx > 0.5 && (p.BI[0] == -1 || p.NI[0] < 0)) {
          is_inflow_boundary = true;
          if (fabs(p.y) < beam_radius) {
            E_new = N_new = 1.0;
            Fx_new = 1.0;
          }
        }
        if (nx < -0.5 && (p.BI[0] == 1 || p.NI[0] > 0)) {
          is_inflow_boundary = true;
          if (fabs(p.y) < beam_radius) {
            E_new = N_new = 1.0;
            Fx_new = -1.0;
          }
        }
        if (ny > 0.5 && (p.BI[1] == -1 || p.NI[1] < 0)) {
          is_inflow_boundary = true;
          if (fabs(p.x) < beam_radius) {
            E_new = N_new = 1.0;
            Fx_new = 0.0;
            Fy_new = 1.0;
          }
        }
        if (ny < -0.5 && (p.BI[1] == 1 || p.NI[1] > 0)) {
          is_inflow_boundary = true;
          if (fabs(p.x) < beam_radius) {
            E_new = N_new = 1.0;
            Fx_new = 0.0;
            Fy_new = -1.0;
          }
        }
        if (!is_inflow_boundary)
          return;

        for (int ig = 0; ig < ncomponents; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rE[i4D] = E_new;
          rFx[i4D] = Fx_new;
          rFy[i4D] = Fy_new;
          rFz[i4D] = Fz_new;
          rN[i4D] = N_new;
        }
      };
  grid.loop_int_device<1, 1, 1>(grid.nghostzones, set_inflow);
  grid.loop_bnd_device<1, 1, 1>(grid.nghostzones, set_inflow);
}

} // namespace nuX_M1_Seeds
