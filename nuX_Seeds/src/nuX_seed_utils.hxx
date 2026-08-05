#ifndef NUX_SEED_UTILS_HXX
#define NUX_SEED_UTILS_HXX

#include "cctk.h"

namespace nuX_Seeds {

inline int radiation_components() {
  int const ncomponents = CCTK_NumVarsInGroup("nuX_M1::rE");
  if (ncomponents <= 0)
    CCTK_ERROR("Could not determine the number of nuX_M1 radiation fields");
  return ncomponents;
}

} // namespace nuX_Seeds

#endif
