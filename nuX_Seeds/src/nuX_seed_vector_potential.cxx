#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"

namespace nuX_Seeds {

using namespace Loop;

extern "C" void nuX_Seeds_InitVectorPotential(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_InitVectorPotential;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_x(cctkGH, {1, 0, 0});
  const GF3D2layout layout_y(cctkGH, {0, 1, 0});
  const GF3D2layout layout_z(cctkGH, {0, 0, 1});

  grid.loop_all_device<1, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Avec_x[layout_x.linear(p.i, p.j, p.k)] = 0.0;
      });
  grid.loop_all_device<0, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Avec_y[layout_y.linear(p.i, p.j, p.k)] = 0.0;
      });
  grid.loop_all_device<0, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Avec_z[layout_z.linear(p.i, p.j, p.k)] = 0.0;
      });
}

} // namespace nuX_Seeds
