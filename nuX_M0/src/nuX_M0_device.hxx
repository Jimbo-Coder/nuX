#ifndef NUX_M0_DEVICE_HXX
#define NUX_M0_DEVICE_HXX

#include <AMReX_Gpu.H>

namespace nuX_M0 {

template <typename T> inline T load_device_scalar(T const *const ptr) {
  T value;
  amrex::Gpu::copy(amrex::Gpu::deviceToHost, ptr, ptr + 1, &value);
  return value;
}

template <typename T>
inline void store_device_scalar(T *const ptr, T const value) {
  amrex::Gpu::copy(amrex::Gpu::hostToDevice, &value, &value + 1, ptr);
}

} // namespace nuX_M0

#endif
