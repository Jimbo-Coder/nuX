#ifndef NUX_M0_SCHEDULE_HXX
#define NUX_M0_SCHEDULE_HXX

namespace nuX_M0 {

constexpr bool iteration_is_due(const int iteration, const int every) {
  return iteration > 0 && every > 0 &&
         (iteration - 1) % every == 0;
}

enum class interp_to_cart_action { skip, clear, interpolate };

constexpr interp_to_cart_action select_interp_to_cart_action(
    const bool levels_synchronized, const bool is_on, const int iteration,
    const int every) {
  if (!levels_synchronized)
    return interp_to_cart_action::skip;
  if (!is_on)
    return interp_to_cart_action::clear;
  return iteration_is_due(iteration, every)
             ? interp_to_cart_action::interpolate
             : interp_to_cart_action::skip;
}

// CarpetX increments cctk_iteration before entering CCTK_EVOL. Preserve the
// THC cadence: compute_every=2 runs on evolution iterations 1, 3, 5, ... .
static_assert(iteration_is_due(1, 2));
static_assert(!iteration_is_due(2, 2));
static_assert(iteration_is_due(3, 2));
static_assert(!iteration_is_due(4, 2));
static_assert(iteration_is_due(5, 2));

// Turning M0 off must clear stale Cartesian absorption independently of the
// transport cadence. Unsynchronized level batches remain untouched.
static_assert(select_interp_to_cart_action(true, false, 2, 2) ==
              interp_to_cart_action::clear);
static_assert(select_interp_to_cart_action(true, false, 3, 2) ==
              interp_to_cart_action::clear);
static_assert(select_interp_to_cart_action(true, true, 2, 2) ==
              interp_to_cart_action::skip);
static_assert(select_interp_to_cart_action(true, true, 3, 2) ==
              interp_to_cart_action::interpolate);
static_assert(select_interp_to_cart_action(false, false, 2, 2) ==
              interp_to_cart_action::skip);

} // namespace nuX_M0

#endif // NUX_M0_SCHEDULE_HXX
