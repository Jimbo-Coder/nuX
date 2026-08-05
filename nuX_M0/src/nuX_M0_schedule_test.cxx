#include <AMReX_Gpu.H>
#include <loop_device.hxx>

#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "nuX_M0_device.hxx"
#include "nuX_M0_schedule.hxx"

namespace nuX_M0 {

using namespace Loop;

extern "C" void nuX_M0_TestScheduleBeforeCart(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_TestScheduleBeforeCart;
  DECLARE_CCTK_PARAMETERS;

  if (!CarpetX_AllLevelsSynchronized())
    return;

  if (compute_every != 2)
    CCTK_ERROR("nuX_M0 schedule regression requires compute_every=2");
  if (!CCTK_Equals(neu_abs_type, "nuX_M0"))
    CCTK_ERROR("nuX_M0 schedule regression requires neu_abs_type=nuX_M0");
  if (wait_until_time > 0.0 || bns_sep_threshold > 0.0)
    CCTK_ERROR("nuX_M0 schedule regression requires unconditional activation");

  const bool due = iteration_is_due(cctk_iteration, compute_every);
  const CCTK_REAL last_transport_time =
      load_device_scalar(nuX_M0_time);
  const bool advanced_this_iteration = last_transport_time == cctk_time;
  if (due != advanced_this_iteration)
    CCTK_VERROR("M0 cadence regression at iteration %d: due=%d, "
                "transport_time=%.17g, cctk_time=%.17g",
                int(cctk_iteration), int(due), double(last_transport_time),
                double(cctk_time));

  if (cctk_iteration != 2)
    return;

  store_device_scalar(nuX_M0_is_on, CCTK_INT(0));

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        abs_number[ijk] = 1.0;
        abs_energy[ijk] = 2.0;
      });
}

extern "C" void nuX_M0_TestScheduleAfterCart(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M0_TestScheduleAfterCart;

  if (!CarpetX_AllLevelsSynchronized() || cctk_iteration != 2)
    return;

  const int npoints = cctk_ash[0] * cctk_ash[1] * cctk_ash[2];
  int *const failed =
      static_cast<int *>(amrex::The_Arena()->alloc(npoints * sizeof(int)));
  amrex::ParallelFor(npoints,
                     [=] CCTK_DEVICE(int i) { failed[i] = 0; });

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        failed[ijk] = abs_number[ijk] != 0.0 || abs_energy[ijk] != 0.0;
      });

  std::vector<int> failed_host(npoints);
  amrex::Gpu::copy(amrex::Gpu::deviceToHost, failed, failed + npoints,
                   failed_host.begin());
  amrex::Gpu::streamSynchronize();
  amrex::The_Arena()->free(failed);

  for (int i = 0; i < npoints; ++i) {
    if (failed_host[i])
      CCTK_VERROR("M0 shutdown regression left stale absorption at local "
                  "point %d",
                  i);
  }
}

} // namespace nuX_M0
