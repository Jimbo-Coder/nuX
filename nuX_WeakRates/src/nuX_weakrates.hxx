#ifndef NUX_WEAKRATES_HXX
#define NUX_WEAKRATES_HXX

#include <cmath>

#include <cctk.h>
#include <loop_device.hxx>

namespace nuX_WeakRates {

inline constexpr int nspecies = 3;
inline constexpr int weakrates_nue = 0;
inline constexpr int weakrates_anue = 1;
inline constexpr int weakrates_nux = 2;

struct EOSState {
  CCTK_REAL rho;
  CCTK_REAL temp;
  CCTK_REAL ye;
  CCTK_REAL particle_mass;
  CCTK_REAL mu_e;
  CCTK_REAL mu_p;
  CCTK_REAL mu_n;
  CCTK_REAL xa;
  CCTK_REAL xh;
  CCTK_REAL xn;
  CCTK_REAL xp;
  CCTK_REAL abar;
  CCTK_REAL zbar;
};

struct RateCoefficients {
  CCTK_REAL eta_0[nspecies];
  CCTK_REAL eta[nspecies];
  CCTK_REAL kappa_0_a[nspecies];
  CCTK_REAL kappa_a[nspecies];
  CCTK_REAL kappa_0_s[nspecies];
  CCTK_REAL kappa_s[nspecies];
};

struct EquilibriumDensities {
  CCTK_REAL number[nspecies];
  CCTK_REAL energy[nspecies];
};

namespace detail {

inline constexpr CCTK_REAL eta0 = 1.0e-3;
inline constexpr CCTK_REAL pi = 3.141592653589793238462643383279502884;
inline constexpr CCTK_REAL mev_to_erg = 1.60217733e-6;
inline constexpr CCTK_REAL clight = 2.99792458e10;
inline constexpr CCTK_REAL me_mev = 0.510998910;
inline constexpr CCTK_REAL me_erg = 8.187108692567103e-7;
inline constexpr CCTK_REAL sigma_0 = 1.76e-44;
inline constexpr CCTK_REAL alpha = 1.23;
inline constexpr CCTK_REAL qnp = 1.293333;
inline constexpr CCTK_REAL hc_mevcm = 1.23984172e-10;
inline constexpr CCTK_REAL cv = 0.5 + 2.0 * 0.23;
inline constexpr CCTK_REAL ca = 0.5;
inline constexpr CCTK_REAL gamma_0 = 5.565e-2;
inline constexpr CCTK_REAL fsc = 1.0 / 137.036;

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
safe_exp(const CCTK_REAL x) {
  return exp(fmin(CCTK_REAL(700), fmax(CCTK_REAL(-700), x)));
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
fermi2p(const CCTK_REAL eta) {
  const CCTK_REAL x = 1.0 / eta;
  return (1.0 / 3.0 + 3.2899 * x * x) /
         (1.0 - safe_exp(-1.8246 * eta));
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
fermi2m(const CCTK_REAL eta) {
  return 2.0 / (1.0 + 0.1092 * safe_exp(0.8908 * eta));
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
fermi2(const CCTK_REAL eta) {
  return eta > eta0 ? eta * eta * eta * fermi2p(eta)
                    : safe_exp(eta) * fermi2m(eta);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
fermi3p(const CCTK_REAL eta) {
  const CCTK_REAL x = 1.0 / eta;
  const CCTK_REAL x2 = x * x;
  return (0.25 + 4.9348 * x2 + 11.3644 * x2 * x2) /
         (1.0 + safe_exp(-1.9039 * eta));
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
fermi3m(const CCTK_REAL eta) {
  return 6.0 / (1.0 + 0.0559 * safe_exp(0.9069 * eta));
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
fermi3(const CCTK_REAL eta) {
  const CCTK_REAL eta2 = eta * eta;
  return eta > eta0 ? eta2 * eta2 * fermi3p(eta)
                    : safe_exp(eta) * fermi3m(eta);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
fermi4p(const CCTK_REAL eta) {
  const CCTK_REAL x = 1.0 / eta;
  const CCTK_REAL x2 = x * x;
  return (0.2 + 6.5797 * x2 + 45.4576 * x2 * x2) /
         (1.0 - safe_exp(-1.9484 * eta));
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
fermi4m(const CCTK_REAL eta) {
  return 24.0 / (1.0 + 0.0287 * safe_exp(0.9257 * eta));
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
fermi4(const CCTK_REAL eta) {
  const CCTK_REAL eta2 = eta * eta;
  return eta > eta0 ? eta2 * eta2 * eta * fermi4p(eta)
                    : safe_exp(eta) * fermi4m(eta);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
fermi5p(const CCTK_REAL eta) {
  const CCTK_REAL x = 1.0 / eta;
  const CCTK_REAL x2 = x * x;
  return (1.0 / 6.0 + 8.2247 * x2 + 113.6439 * x2 * x2 +
          236.5323 * x2 * x2 * x2) /
         (1.0 + safe_exp(-1.9727 * eta));
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
fermi5m(const CCTK_REAL eta) {
  return 120.0 / (1.0 + 0.0147 * safe_exp(0.9431 * eta));
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
fermi4o2(const CCTK_REAL eta) {
  return eta > eta0 ? eta * eta * fermi4p(eta) / fermi2p(eta)
                    : fermi4m(eta) / fermi2m(eta);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
fermi4o3(const CCTK_REAL eta) {
  return eta > eta0 ? eta * fermi4p(eta) / fermi3p(eta)
                    : fermi4m(eta) / fermi3m(eta);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
fermi5o3(const CCTK_REAL eta) {
  return eta > eta0 ? eta * eta * fermi5p(eta) / fermi3p(eta)
                    : fermi5m(eta) / fermi3m(eta);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
fermi5o4(const CCTK_REAL eta) {
  return eta > eta0 ? eta * fermi5p(eta) / fermi4p(eta)
                    : fermi5m(eta) / fermi4m(eta);
}

struct DerivedState {
  CCTK_REAL nb;
  CCTK_REAL eta_nue;
  CCTK_REAL eta_anue;
  CCTK_REAL eta_e;
  CCTK_REAL eta_np;
  CCTK_REAL eta_pn;
  CCTK_REAL xn;
  CCTK_REAL xp;
  CCTK_REAL xh;
  CCTK_REAL abar;
  CCTK_REAL zbar;
};

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline DerivedState
derive_state(const EOSState &eos, const CCTK_REAL low_density_threshold) {
  const CCTK_REAL mass_cgs =
      eos.particle_mass * mev_to_erg / (clight * clight);
  const CCTK_REAL nb = eos.rho / mass_cgs;
  const CCTK_REAL eta_nue = (eos.mu_p + eos.mu_e - eos.mu_n) / eos.temp;
  const CCTK_REAL eta_anue = -eta_nue;
  const CCTK_REAL eta_e = eos.mu_e / eos.temp;
  const CCTK_REAL eta_hat =
      (eos.mu_n - eos.mu_p - qnp) / eos.temp;
  const CCTK_REAL xn = fmax(CCTK_REAL(0), eos.xn);
  const CCTK_REAL xp = fmax(CCTK_REAL(0), eos.xp);
  const CCTK_REAL xh = fmax(CCTK_REAL(0), eos.xh);
  const CCTK_REAL abar = fmax(CCTK_REAL(0), eos.abar);
  const CCTK_REAL zbar = fmax(CCTK_REAL(0), eos.zbar);

  CCTK_REAL eta_np =
      nb * (xp - xn) / expm1(-eta_hat);
  CCTK_REAL eta_pn =
      nb * (xn - xp) / expm1(eta_hat);
  if (eos.rho < low_density_threshold) {
    eta_pn = nb * xp;
    eta_np = nb * xn;
  }

  return {nb,
          eta_nue,
          eta_anue,
          eta_e,
          fmax(CCTK_REAL(0), eta_np),
          fmax(CCTK_REAL(0), eta_pn),
          xn,
          xp,
          xh,
          abar,
          zbar};
}

} // namespace detail

class WeakRates {
public:
  bool include_beta;
  bool include_pair;
  bool include_plasmon;
  bool include_bremsstrahlung;
  bool include_elastic;
  CCTK_REAL beta_low_density_threshold;

  CCTK_HOST void init();

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline
  RateCoefficients compute_rates(const EOSState &eos) const {
    using namespace detail;

    RateCoefficients rates = {};
    const DerivedState state =
        derive_state(eos, beta_low_density_threshold);
    const CCTK_REAL temp = eos.temp;
    const CCTK_REAL temp2_erg =
        temp * temp * mev_to_erg * mev_to_erg;

    CCTK_REAL absorption_zeta[nspecies] = {};
    if (include_beta) {
      const CCTK_REAL abs_zeta =
          (1.0 + 3.0 * alpha * alpha) * 0.25 * sigma_0 /
          (me_erg * me_erg);
      absorption_zeta[weakrates_nue] =
          state.eta_np * abs_zeta /
          (1.0 + safe_exp(state.eta_e - fermi5o4(state.eta_nue)));
      absorption_zeta[weakrates_anue] =
          state.eta_pn * abs_zeta /
          (1.0 + safe_exp(-state.eta_e - fermi5o4(state.eta_anue)));
    }

    CCTK_REAL scattering_zeta = 0.0;
    if (include_elastic) {
      const CCTK_REAL scatter_p =
          state.nb * (1.0 + 5.0 * alpha * alpha) / 24.0 * sigma_0 /
          (me_erg * me_erg);
      const CCTK_REAL scatter_n =
          state.nb * (4.0 * (cv - 1.0) * (cv - 1.0) +
                      5.0 * alpha * alpha) /
          24.0 * sigma_0 / (me_erg * me_erg);
      const CCTK_REAL scatter_h =
          state.abar > 0.0
              ? state.nb * 0.0625 * sigma_0 / (me_erg * me_erg) *
                    state.abar *
                    (1.0 - state.zbar / state.abar) *
                    (1.0 - state.zbar / state.abar)
              : 0.0;
      scattering_zeta = state.xn * scatter_n + state.xp * scatter_p +
                        state.xh * scatter_h;
    }

    const CCTK_REAL eta[nspecies] = {state.eta_nue, state.eta_anue, 0.0};
    for (int isp = 0; isp < nspecies; ++isp) {
      rates.kappa_0_a[isp] =
          absorption_zeta[isp] * temp2_erg * fermi4o2(eta[isp]);
      rates.kappa_a[isp] =
          absorption_zeta[isp] * temp2_erg * fermi5o3(eta[isp]);
      rates.kappa_0_s[isp] =
          scattering_zeta * temp2_erg * fermi4o2(eta[isp]);
      rates.kappa_s[isp] =
          scattering_zeta * temp2_erg * fermi5o3(eta[isp]);
    }

    CCTK_REAL r_beta_nue = 0.0;
    CCTK_REAL r_beta_anue = 0.0;
    CCTK_REAL q_beta_nue = 0.0;
    CCTK_REAL q_beta_anue = 0.0;
    if (include_beta) {
      const CCTK_REAL beta =
          pi * clight * (1.0 + 3.0 * alpha * alpha) * sigma_0 /
          (hc_mevcm * hc_mevcm * hc_mevcm * me_mev * me_mev);
      const CCTK_REAL temp5 = temp * temp * temp * temp * temp;
      const CCTK_REAL block_nue =
          1.0 + safe_exp(state.eta_nue - fermi5o4(state.eta_e));
      const CCTK_REAL block_anue =
          1.0 + safe_exp(state.eta_anue - fermi5o4(-state.eta_e));
      r_beta_nue =
          beta * state.eta_pn * temp5 * fermi4(state.eta_e) / block_nue;
      r_beta_anue = beta * state.eta_np * temp5 *
                    fermi4(-state.eta_e) / block_anue;
      q_beta_nue = r_beta_nue * temp * fermi5o4(state.eta_e);
      q_beta_anue = r_beta_anue * temp * fermi5o4(-state.eta_e);
    }

    CCTK_REAL r_pair_nue = 0.0;
    CCTK_REAL r_pair_anue = 0.0;
    CCTK_REAL r_pair_nux = 0.0;
    CCTK_REAL q_pair_nue = 0.0;
    CCTK_REAL q_pair_anue = 0.0;
    CCTK_REAL q_pair_nux = 0.0;
    if (include_pair) {
      const CCTK_REAL temp4 = temp * temp * temp * temp;
      const CCTK_REAL enr_p = 8.0 * pi /
                              (hc_mevcm * hc_mevcm * hc_mevcm) * temp4 *
                              fermi3(state.eta_e);
      const CCTK_REAL enr_m = 8.0 * pi /
                              (hc_mevcm * hc_mevcm * hc_mevcm) * temp4 *
                              fermi3(-state.eta_e);
      const CCTK_REAL pair_const =
          sigma_0 * clight / (me_mev * me_mev) * enr_m * enr_p;
      const CCTK_REAL pair_energy =
          0.5 * temp *
          (fermi4o3(-state.eta_e) + fermi4o3(state.eta_e));
      const CCTK_REAL blocking_energy = 0.5 *
          (fermi4o3(state.eta_e) + fermi4o3(-state.eta_e));
      const CCTK_REAL block_nue =
          1.0 + safe_exp(state.eta_nue - blocking_energy);
      const CCTK_REAL block_anue =
          1.0 + safe_exp(state.eta_anue - blocking_energy);
      const CCTK_REAL block_nux = 1.0 + safe_exp(-blocking_energy);

      r_pair_nue = pair_const /
                   (36.0 * block_nue * block_anue) *
                   ((cv - ca) * (cv - ca) + (cv + ca) * (cv + ca));
      r_pair_anue = r_pair_nue;
      r_pair_nux = pair_const / (9.0 * block_nux * block_nux) *
                   ((cv - ca) * (cv - ca) +
                    (cv + ca - 2.0) * (cv + ca - 2.0));
      q_pair_nue = r_pair_nue * pair_energy;
      q_pair_anue = r_pair_anue * pair_energy;
      q_pair_nux = r_pair_nux * pair_energy;
    }

    CCTK_REAL r_plasmon_nue = 0.0;
    CCTK_REAL r_plasmon_anue = 0.0;
    CCTK_REAL r_plasmon_nux = 0.0;
    CCTK_REAL q_plasmon_nue = 0.0;
    CCTK_REAL q_plasmon_anue = 0.0;
    CCTK_REAL q_plasmon_nux = 0.0;
    if (include_plasmon) {
      const CCTK_REAL gamma =
          gamma_0 * sqrt((pi * pi + 3.0 * state.eta_e * state.eta_e) / 3.0);
      const CCTK_REAL blocking_energy =
          1.0 + 0.5 * gamma * gamma / (1.0 + gamma);
      const CCTK_REAL block_nue =
          1.0 + safe_exp(state.eta_nue - blocking_energy);
      const CCTK_REAL block_anue =
          1.0 + safe_exp(state.eta_anue - blocking_energy);
      const CCTK_REAL block_nux = 1.0 + safe_exp(-blocking_energy);
      const CCTK_REAL temp2 = temp * temp;
      const CCTK_REAL temp4 = temp2 * temp2;
      const CCTK_REAL temp8 = temp4 * temp4;
      const CCTK_REAL gamma2 = gamma * gamma;
      const CCTK_REAL gamma6 = gamma2 * gamma2 * gamma2;
      const CCTK_REAL gamma_const =
          pi * pi * pi * sigma_0 * clight * temp8 /
          (me_mev * me_mev * 3.0 * fsc *
           pow(hc_mevcm, CCTK_REAL(6))) *
          gamma6 * safe_exp(-gamma) * (1.0 + gamma);
      const CCTK_REAL plasmon_energy =
          0.5 * temp * (2.0 + gamma2 / (1.0 + gamma));

      const CCTK_REAL r_gamma =
          cv * cv * gamma_const / (block_nue * block_anue);
      r_plasmon_nue = r_gamma;
      r_plasmon_anue = r_gamma;
      r_plasmon_nux = (cv - 1.0) * (cv - 1.0) * 4.0 * gamma_const /
                      (block_nux * block_nux);
      q_plasmon_nue = r_plasmon_nue * plasmon_energy;
      q_plasmon_anue = r_plasmon_anue * plasmon_energy;
      q_plasmon_nux = r_plasmon_nux * plasmon_energy;
    }

    CCTK_REAL r_brem = 0.0;
    CCTK_REAL q_brem = 0.0;
    if (include_bremsstrahlung) {
      r_brem = 0.231 * (2.0778e2 / mev_to_erg) * 0.5 *
               (state.xn * state.xn + state.xp * state.xp +
                28.0 / 3.0 * state.xn * state.xp) *
               eos.rho * eos.rho * pow(temp, CCTK_REAL(4.5));
      q_brem = r_brem * temp / 0.231 * 0.504;
    }

    rates.eta_0[weakrates_nue] =
        r_beta_nue + r_pair_nue + r_plasmon_nue + r_brem;
    rates.eta_0[weakrates_anue] =
        r_beta_anue + r_pair_anue + r_plasmon_anue + r_brem;
    rates.eta_0[weakrates_nux] =
        r_pair_nux + r_plasmon_nux + 4.0 * r_brem;
    rates.eta[weakrates_nue] =
        q_beta_nue + q_pair_nue + q_plasmon_nue + q_brem;
    rates.eta[weakrates_anue] =
        q_beta_anue + q_pair_anue + q_plasmon_anue + q_brem;
    rates.eta[weakrates_nux] =
        q_pair_nux + q_plasmon_nux + 4.0 * q_brem;

    return rates;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline
  EquilibriumDensities equilibrium_densities(const EOSState &eos) const {
    using namespace detail;

    EquilibriumDensities densities = {};
    const CCTK_REAL eta_nue =
        (eos.mu_p + eos.mu_e - eos.mu_n) / eos.temp;
    const CCTK_REAL eta[nspecies] = {eta_nue, -eta_nue, 0.0};
    const CCTK_REAL temp2 = eos.temp * eos.temp;
    const CCTK_REAL temp3 = temp2 * eos.temp;
    const CCTK_REAL temp4 = temp3 * eos.temp;
    const CCTK_REAL prefactor = 4.0 * pi /
                                (hc_mevcm * hc_mevcm * hc_mevcm);
    for (int isp = 0; isp < nspecies; ++isp) {
      const CCTK_REAL multiplicity = isp == weakrates_nux ? 4.0 : 1.0;
      densities.number[isp] =
          multiplicity * prefactor * temp3 * fermi2(eta[isp]);
      densities.energy[isp] =
          multiplicity * prefactor * temp4 * fermi3(eta[isp]);
    }
    return densities;
  }
};

extern WeakRates *global_weakrates;

} // namespace nuX_WeakRates

#endif
