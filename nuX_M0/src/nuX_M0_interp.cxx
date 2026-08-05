#include <AMReX_Gpu.H>

#include <array>
#include <cassert>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

#include "nuX_M0_device.hxx"
#include "nuX_M0_kernel.hxx"
#include "nuX_M0_schedule.hxx"

namespace nuX_M0 {

using namespace Loop;

extern "C" void nuX_M0_InterpToSph(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_InterpToSph;
  DECLARE_CCTK_PARAMETERS;

  if (!CarpetX_AllLevelsSynchronized() ||
      !iteration_is_due(cctk_iteration, compute_every) ||
      load_device_scalar(nuX_M0_is_on) == 0) {
    return;
  }

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_M0_InterpToSph");
  }

  MPI_Datatype mpi_real;
  if (sizeof(CCTK_REAL) == sizeof(float))
    mpi_real = MPI_FLOAT;
  else if (sizeof(CCTK_REAL) == sizeof(double))
    mpi_real = MPI_DOUBLE;
  else if (sizeof(CCTK_REAL) == sizeof(long double))
    mpi_real = MPI_LONG_DOUBLE;
  else
    CCTK_ERROR("Unsupported CCTK_REAL size for the M0 MPI broadcast");

  int const group_id = CCTK_GroupIndex("nuX_M0::nuX_M0_grid_vars");
  cGroupDynamicData group_data;
  int const group_ierr = CCTK_GroupDynamicData(cctkGH, group_id, &group_data);
  assert(!group_ierr);
  int const npts = group_data.ash[0] * group_data.ash[1];

  constexpr int ndim = 3;
  constexpr int nvars = 22;
  std::array<char const *, nvars> const variable_names = {
      "ADMBaseX::alp",
      "ADMBaseX::betax",
      "ADMBaseX::betay",
      "ADMBaseX::betaz",
      "ADMBaseX::gxx",
      "ADMBaseX::gxy",
      "ADMBaseX::gxz",
      "ADMBaseX::gyy",
      "ADMBaseX::gyz",
      "ADMBaseX::gzz",
      "HydroBaseX::rho",
      "nuX_M0::zvec_cartx",
      "nuX_M0::zvec_carty",
      "nuX_M0::zvec_cartz",
      "HydroBaseX::temperature",
      "HydroBaseX::Ye",
      "nuX_LeakageBase::optd_0_nue",
      "nuX_LeakageBase::optd_0_nua",
      "nuX_LeakageBase::optd_0_nux",
      "nuX_LeakageBase::optd_1_nue",
      "nuX_LeakageBase::optd_1_nua",
      "nuX_LeakageBase::optd_1_nux",
  };

  std::array<CCTK_REAL *, nvars> const ray_arrays = {
      nuX_M0_alp,         nuX_M0_betax,       nuX_M0_betay,
      nuX_M0_betaz,       nuX_M0_gxx,         nuX_M0_gxy,
      nuX_M0_gxz,         nuX_M0_gyy,         nuX_M0_gyz,
      nuX_M0_gzz,         nuX_M0_rho,         nuX_M0_zvecx,
      nuX_M0_zvecy,       nuX_M0_zvecz,       nuX_M0_temp,
      nuX_M0_Ye,          nuX_M0_optd_0_nue,  nuX_M0_optd_0_nua,
      nuX_M0_optd_0_nux,  nuX_M0_optd_1_nue,  nuX_M0_optd_1_nua,
      nuX_M0_optd_1_nux,
  };

  CCTK_REAL *const coords_host = static_cast<CCTK_REAL *>(
      amrex::The_Pinned_Arena()->alloc(ndim * npts * sizeof(CCTK_REAL)));
  CCTK_REAL *const values_host = static_cast<CCTK_REAL *>(
      amrex::The_Pinned_Arena()->alloc(nvars * npts * sizeof(CCTK_REAL)));

  int const rank = CCTK_MyProc(cctkGH);
  if (rank == 0) {
    std::array<CCTK_REAL const *, ndim> const coords_device = {
        nuX_M0_x, nuX_M0_y, nuX_M0_z};
    for (int d = 0; d < ndim; ++d) {
      amrex::Gpu::copy(amrex::Gpu::deviceToHost, coords_device[d],
                       coords_device[d] + npts, coords_host + d * npts);
    }
  }

  std::array<void const *, ndim> const interp_coords = {
      coords_host, coords_host + npts, coords_host + 2 * npts};
  std::array<CCTK_INT, nvars> variable_indices;
  std::array<CCTK_INT, nvars> output_types;
  std::array<void *, nvars> output_arrays;
  for (int ivar = 0; ivar < nvars; ++ivar) {
    variable_indices[ivar] = CCTK_VarIndex(variable_names[ivar]);
    if (variable_indices[ivar] < 0)
      CCTK_VERROR("Could not find interpolation variable %s",
                  variable_names[ivar]);
    output_types[ivar] = CCTK_VARIABLE_REAL;
    output_arrays[ivar] = values_host + ivar * npts;
  }

  int const operator_handle = CCTK_InterpHandle("CarpetX");
  if (operator_handle < 0)
    CCTK_VERROR("Could not obtain the CarpetX interpolator handle: %d",
                operator_handle);

  int interpolation_order_type;
  const void *const interpolation_order_p = CCTK_ParameterGet(
      "interpolation_order", "CarpetX", &interpolation_order_type);
  assert(interpolation_order_p);
  assert(interpolation_order_type == PARAMETER_INT);
  const CCTK_INT operator_order =
      *static_cast<const CCTK_INT *>(interpolation_order_p);

  int const options_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  if (options_handle < 0)
    CCTK_VERROR("Could not create the M0 interpolation options table: %d",
                options_handle);
  int const table_ierr =
      Util_TableSetInt(options_handle, operator_order, "order");
  if (table_ierr < 0)
    CCTK_VERROR("Could not set the M0 interpolation order: %d", table_ierr);

  int const interp_ierr = CCTK_InterpGridArrays(
      cctkGH, ndim, operator_handle, options_handle, 0,
      rank == 0 ? npts : 0, CCTK_VARIABLE_REAL, interp_coords.data(), nvars,
      variable_indices.data(), nvars, output_types.data(),
      output_arrays.data());
  int const destroy_ierr = Util_TableDestroy(options_handle);
  if (interp_ierr < 0)
    CCTK_VERROR("M0 ray-grid interpolation failed: %d", interp_ierr);
  if (destroy_ierr < 0)
    CCTK_VERROR("Could not destroy the M0 interpolation table: %d",
                destroy_ierr);

  int const bcast_ierr = MPI_Bcast(values_host, nvars * npts, mpi_real, 0,
                                   MPI_COMM_WORLD);
  if (bcast_ierr != MPI_SUCCESS)
    CCTK_VERROR("Could not broadcast M0 ray-grid data: %d", bcast_ierr);

  for (int ivar = 0; ivar < nvars; ++ivar) {
    CCTK_REAL const *const source = values_host + ivar * npts;
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, source, source + npts,
                     ray_arrays[ivar]);
  }

  amrex::The_Pinned_Arena()->free(values_host);
  amrex::The_Pinned_Arena()->free(coords_host);
}

extern "C" void nuX_M0_InterpToCart(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_InterpToCart;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});

  const interp_to_cart_action action = select_interp_to_cart_action(
      CarpetX_AllLevelsSynchronized(),
      load_device_scalar(nuX_M0_is_on) != 0, cctk_iteration, compute_every);
  if (action == interp_to_cart_action::skip)
    return;

  if (action == interp_to_cart_action::clear) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const int ijk = layout_cc.linear(p.i, p.j, p.k);
          abs_number[ijk] = 0.0;
          abs_energy[ijk] = 0.0;
        });
    return;
  }

  if (verbose && CCTK_MyProc(cctkGH) == 0) {
    CCTK_INFO("nuX_M0_InterpToCart");
  }

  SphericalGrid const sph_grid = *M0Grid;
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        CCTK_REAL r, theta, phi;
        nuX_Utils::sph_grid::coord_cart_to_sph(p.x, p.y, p.z, r, theta, phi);
        if (r > sph_grid.rmax)
          r = sph_grid.rmax;
        abs_number[ijk] =
            sample_spherical_linear(sph_grid, nuX_M0_abs_number, r, theta, phi);
        abs_energy[ijk] =
            sample_spherical_linear(sph_grid, nuX_M0_abs_energy, r, theta, phi);
      });
}

} // namespace nuX_M0
