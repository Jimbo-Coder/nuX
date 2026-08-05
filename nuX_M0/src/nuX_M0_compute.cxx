#include <AMReX_Gpu.H>

#include <array>
#include <cassert>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "nuX_M0_device.hxx"
#include "nuX_M0_kernel.hxx"
#include "nuX_M0_schedule.hxx"
#include "nuX_baryon_mass.hxx"

namespace nuX_M0 {

extern "C" void nuX_M0_Compute(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_Compute;
  DECLARE_CCTK_PARAMETERS;

  if (!CarpetX_AllLevelsSynchronized() ||
      !iteration_is_due(cctk_iteration, compute_every) ||
      load_device_scalar(nuX_M0_is_on) == 0) {
    return;
  }

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_M0_Compute");
  }

  CCTK_REAL const dt = cctk_time - load_device_scalar(nuX_M0_time);
  assert(dt >= 0.0);

  int const group_id = CCTK_GroupIndex("nuX_M0::nuX_M0_grid_vars");
  cGroupDynamicData group_data;
  int const ierr = CCTK_GroupDynamicData(cctkGH, group_id, &group_data);
  assert(!ierr);

  int const nrad_local = group_data.ash[0];
  int const nray_local = group_data.ash[1];
  int const npts = nrad_local * nray_local;

  CCTK_REAL const mb = nuX_Utils::AverageBaryonMass(particle_mass);
  SphericalGrid const grid = *M0Grid;

  amrex::ParallelFor(npts, [=] CCTK_DEVICE(int const idx) {
    nuX_M0_N_nue_old[idx] = nuX_M0_N_nue[idx];
    nuX_M0_N_nua_old[idx] = nuX_M0_N_nua[idx];
    nuX_M0_N_nux_old[idx] = nuX_M0_N_nux[idx];
    nuX_M0_E_nue_old[idx] = nuX_M0_E_nue[idx];
    nuX_M0_E_nua_old[idx] = nuX_M0_E_nua[idx];
    nuX_M0_E_nux_old[idx] = nuX_M0_E_nux[idx];
  });

  constexpr int nfluxes = 6;
  CCTK_REAL *const ray_fluxes = static_cast<CCTK_REAL *>(
      amrex::The_Arena()->alloc(nfluxes * nray_local * sizeof(CCTK_REAL)));
  CCTK_REAL *const ray_nue_num = ray_fluxes;
  CCTK_REAL *const ray_nua_num = ray_fluxes + nray_local;
  CCTK_REAL *const ray_nux_num = ray_fluxes + 2 * nray_local;
  CCTK_REAL *const ray_nue_ene = ray_fluxes + 3 * nray_local;
  CCTK_REAL *const ray_nua_ene = ray_fluxes + 4 * nray_local;
  CCTK_REAL *const ray_nux_ene = ray_fluxes + 5 * nray_local;

  amrex::ParallelFor(nray_local, [=] CCTK_DEVICE(int const iray) {
    int const offset = index(nrad_local, 0, iray);

    for (int irad = 0; irad < nrad_local; ++irad) {
      int const ij = offset + irad;
      nuX_M0_mask[ij] = 0;
      if (excision && nuX_M0_alp[ij] <= 0.0) {
        nuX_M0_mask[ij] = 1;
      }

      CCTK_REAL kr = 0.0;
      rad_null(grid, irad, iray, nuX_M0_alp[ij], nuX_M0_betax[ij],
               nuX_M0_betay[ij], nuX_M0_betaz[ij], nuX_M0_gxx[ij],
               nuX_M0_gxy[ij], nuX_M0_gxz[ij], nuX_M0_gyy[ij], nuX_M0_gyz[ij],
               nuX_M0_gzz[ij], nuX_M0_zvecx[ij], nuX_M0_zvecy[ij],
               nuX_M0_zvecz[ij], nuX_M0_mask[ij], nuX_M0_kt[ij], kr,
               nuX_M0_chi[ij], nuX_M0_sdetg[ij]);
      nuX_M0_flux_fac[ij] = nuX_M0_kt[ij] > 0.0 ? kr / nuX_M0_kt[ij] : 0.0;
    }

    for (int irad = 0; irad < nrad_local; ++irad) {
      int const ij = offset + irad;
      if (nuX_M0_mask[ij]) {
        nuX_M0_abs_nue[ij] = 0.0;
        nuX_M0_abs_nua[ij] = 0.0;
        nuX_M0_ndens_nue[ij] = 0.0;
        nuX_M0_ndens_nua[ij] = 0.0;
        nuX_M0_ndens_nux[ij] = 0.0;
        nuX_M0_eave_nue[ij] = 0.0;
        nuX_M0_eave_nua[ij] = 0.0;
        nuX_M0_eave_nux[ij] = 0.0;
        nuX_M0_abs_number[ij] = 0.0;
        nuX_M0_abs_energy[ij] = 0.0;
        nuX_M0_N_nue_old[ij] = 0.0;
        nuX_M0_N_nua_old[ij] = 0.0;
        nuX_M0_N_nux_old[ij] = 0.0;
        nuX_M0_E_nue_old[ij] = 0.0;
        nuX_M0_E_nua_old[ij] = 0.0;
        nuX_M0_E_nux_old[ij] = 0.0;
        nuX_M0_N_nue[ij] = 0.0;
        nuX_M0_N_nua[ij] = 0.0;
        nuX_M0_N_nux[ij] = 0.0;
        nuX_M0_E_nue[ij] = 0.0;
        nuX_M0_E_nua[ij] = 0.0;
        nuX_M0_E_nux[ij] = 0.0;
        continue;
      }

      CCTK_REAL abs_nue = nuX_M0_abs_0_nue[ij];
      CCTK_REAL abs_nua = nuX_M0_abs_0_nua[ij];
      if (use_reduced_opacity) {
        abs_nue *= ::exp(-nuX_M0_optd_0_nue[ij]);
        abs_nua *= ::exp(-nuX_M0_optd_0_nua[ij]);
      }
      nuX_M0_abs_nue[ij] = abs_nue;
      nuX_M0_abs_nua[ij] = abs_nua;
    }

    evol_density_slice(grid, dt, &nuX_M0_mask[offset], &nuX_M0_flux_fac[offset],
                       &nuX_M0_sdetg[offset], &nuX_M0_kt[offset],
                       &nuX_M0_R_nue[offset], &nuX_M0_abs_nue[offset],
                       &nuX_M0_N_nue_old[offset], &nuX_M0_N_nue[offset],
                       &nuX_M0_ndens_nue[offset]);
    evol_density_slice(grid, dt, &nuX_M0_mask[offset], &nuX_M0_flux_fac[offset],
                       &nuX_M0_sdetg[offset], &nuX_M0_kt[offset],
                       &nuX_M0_R_nua[offset], &nuX_M0_abs_nua[offset],
                       &nuX_M0_N_nua_old[offset], &nuX_M0_N_nua[offset],
                       &nuX_M0_ndens_nua[offset]);
    evol_density_slice(grid, dt, &nuX_M0_mask[offset], &nuX_M0_flux_fac[offset],
                       &nuX_M0_sdetg[offset], &nuX_M0_kt[offset],
                       &nuX_M0_R_nux[offset], nullptr,
                       &nuX_M0_N_nux_old[offset], &nuX_M0_N_nux[offset],
                       &nuX_M0_ndens_nux[offset]);

    evol_energy_slice(grid, dt, &nuX_M0_mask[offset], &nuX_M0_flux_fac[offset],
                      &nuX_M0_kt[offset], &nuX_M0_chi[offset],
                      &nuX_M0_R_nue[offset], &nuX_M0_Q_nue[offset],
                      &nuX_M0_ndens_nue[offset], &nuX_M0_E_nue_old[offset],
                      &nuX_M0_E_nue[offset], &nuX_M0_eave_nue[offset]);
    evol_energy_slice(grid, dt, &nuX_M0_mask[offset], &nuX_M0_flux_fac[offset],
                      &nuX_M0_kt[offset], &nuX_M0_chi[offset],
                      &nuX_M0_R_nua[offset], &nuX_M0_Q_nua[offset],
                      &nuX_M0_ndens_nua[offset], &nuX_M0_E_nua_old[offset],
                      &nuX_M0_E_nua[offset], &nuX_M0_eave_nua[offset]);
    evol_energy_slice(grid, dt, &nuX_M0_mask[offset], &nuX_M0_flux_fac[offset],
                      &nuX_M0_kt[offset], &nuX_M0_chi[offset],
                      &nuX_M0_R_nux[offset], &nuX_M0_Q_nux[offset],
                      &nuX_M0_ndens_nux[offset], &nuX_M0_E_nux_old[offset],
                      &nuX_M0_E_nux[offset], &nuX_M0_eave_nux[offset]);

    for (int irad = 0; irad < nrad_local; ++irad) {
      int const ij = offset + irad;
      nuX_M0_abs_number[ij] =
          mb * (nuX_M0_abs_nue[ij] * nuX_M0_ndens_nue[ij] -
                nuX_M0_abs_nua[ij] * nuX_M0_ndens_nua[ij]);
      nuX_M0_abs_energy[ij] =
          nuX_M0_abs_nue[ij] * nuX_M0_ndens_nue[ij] * nuX_M0_eave_nue[ij] +
          nuX_M0_abs_nua[ij] * nuX_M0_ndens_nua[ij] * nuX_M0_eave_nua[ij];
    }

    int const itheta = nuX_Utils::sph_grid::get_itheta(grid, iray);
    int const iphi = nuX_Utils::sph_grid::get_iphi(grid, iray);
    CCTK_REAL dS = grid.dtheta * grid.dphi;
    if (iphi == 0 || iphi == grid.nphi - 1)
      dS *= 0.5;
    if (itheta == 0 || itheta == grid.ntheta - 1)
      dS /= grid.nphi;

    int const iout = offset + (nrad_local - 1);
    CCTK_REAL const wt = nuX_M0_flux_fac[iout] * dS;
    ray_nue_num[iray] = nuX_M0_N_nue[iout] * wt;
    ray_nua_num[iray] = nuX_M0_N_nua[iout] * wt;
    ray_nux_num[iray] = nuX_M0_N_nux[iout] * wt;
    ray_nue_ene[iray] = nuX_M0_N_nue[iout] * nuX_M0_eave_nue[iout] * wt;
    ray_nua_ene[iray] = nuX_M0_N_nua[iout] * nuX_M0_eave_nua[iout] * wt;
    ray_nux_ene[iray] = nuX_M0_N_nux[iout] * nuX_M0_eave_nux[iout] * wt;
  });

  CCTK_REAL *const ray_fluxes_host = static_cast<CCTK_REAL *>(
      amrex::The_Pinned_Arena()->alloc(nfluxes * nray_local *
                                       sizeof(CCTK_REAL)));
  amrex::Gpu::copy(amrex::Gpu::deviceToHost, ray_fluxes,
                   ray_fluxes + nfluxes * nray_local, ray_fluxes_host);

  std::array<CCTK_REAL, nfluxes> fluxes{};
  for (int iray = 0; iray < nray_local; ++iray) {
    for (int iflux = 0; iflux < nfluxes; ++iflux)
      fluxes[iflux] += ray_fluxes_host[iflux * nray_local + iray];
  }

  store_device_scalar(nuX_M0_nue_num_flux, fluxes[0]);
  store_device_scalar(nuX_M0_nua_num_flux, fluxes[1]);
  store_device_scalar(nuX_M0_nux_num_flux, fluxes[2]);
  store_device_scalar(nuX_M0_nue_ene_flux, fluxes[3]);
  store_device_scalar(nuX_M0_nua_ene_flux, fluxes[4]);
  store_device_scalar(nuX_M0_nux_ene_flux, fluxes[5]);

  amrex::The_Pinned_Arena()->free(ray_fluxes_host);
  amrex::The_Arena()->free(ray_fluxes);

  store_device_scalar(nuX_M0_time, CCTK_REAL(cctk_time));
}

} // namespace nuX_M0
