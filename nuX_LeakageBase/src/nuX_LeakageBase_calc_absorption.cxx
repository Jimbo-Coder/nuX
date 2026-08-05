#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_LeakageBase {

using namespace Loop;

extern "C" void nuX_LeakageBase_NoAbsorption(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_LeakageBase_NoAbsorption;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        abs_number[ijk] = 0.0;
        abs_energy[ijk] = 0.0;
      });
}

} // namespace nuX_LeakageBase
