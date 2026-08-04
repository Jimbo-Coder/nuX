#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>
#include <float.h>

namespace nuX_FakeRates {

struct FakeOpacityCoefficients {
  CCTK_REAL eta_0[4];
  CCTK_REAL eta[4];
  CCTK_REAL kappa_0_a[4];
  CCTK_REAL kappa_a[4];
  CCTK_REAL kappa_s[4];
};

class FakeRatesDef {
public:
  CCTK_REAL kabs_nue;
  CCTK_REAL kabs_nua;
  CCTK_REAL kabs_nux;
  CCTK_REAL kabs_anux;

  CCTK_REAL kscat_nue;
  CCTK_REAL kscat_nua;
  CCTK_REAL kscat_nux;
  CCTK_REAL kscat_anux;

  CCTK_REAL et_nue;
  CCTK_REAL et_nua;
  CCTK_REAL et_nux;
  CCTK_REAL et_anux;

  CCTK_HOST void init();

  CCTK_DEVICE inline FakeOpacityCoefficients
  ComputeFakeOpacities(const CCTK_REAL rho) {

    FakeOpacityCoefficients coefficients = {0};
    coefficients.eta_0[0] = rho * et_nue;
    coefficients.eta_0[1] = rho * et_nua;
    coefficients.eta_0[2] = rho * et_nux;
    coefficients.eta_0[3] = rho * et_anux;

    coefficients.kappa_0_a[0] = rho * kabs_nue;
    coefficients.kappa_0_a[1] = rho * kabs_nua;
    coefficients.kappa_0_a[2] = rho * kabs_nux;
    coefficients.kappa_0_a[3] = rho * kabs_anux;

    coefficients.eta[0] = rho * et_nue;
    coefficients.eta[1] = rho * et_nua;
    coefficients.eta[2] = rho * et_nux;
    coefficients.eta[3] = rho * et_anux;

    coefficients.kappa_a[0] = rho * kabs_nue;
    coefficients.kappa_a[1] = rho * kabs_nua;
    coefficients.kappa_a[2] = rho * kabs_nux;
    coefficients.kappa_a[3] = rho * kabs_anux;

    coefficients.kappa_s[0] = rho * kscat_nue;
    coefficients.kappa_s[1] = rho * kscat_nua;
    coefficients.kappa_s[2] = rho * kscat_nux;
    coefficients.kappa_s[3] = rho * kscat_anux;
    return coefficients;
  }

  CCTK_DEVICE CCTK_HOST inline void
  FakeNeutrinoDens(CCTK_REAL rho, CCTK_REAL &num_nue, CCTK_REAL &num_nua,
                   CCTK_REAL &num_nux, CCTK_REAL &ene_nue, CCTK_REAL &ene_nua,
                   CCTK_REAL &ene_nux) {

    if (rho * kabs_nue > FLT_EPSILON * et_nue) {
      num_nue = et_nue / (rho * kabs_nue);
      ene_nue = et_nue / (rho * kabs_nue);
    } else {
      num_nue = 1.0;
      ene_nue = 1.0;
    }

    if (rho * kabs_nua > FLT_EPSILON * et_nua) {
      num_nua = et_nua / (rho * kabs_nua);
      ene_nua = et_nua / (rho * kabs_nua);
    } else {
      num_nua = 1.0;
      ene_nua = 1.0;
    }

    num_nux = 1.0;
    ene_nux = 1.0;
  }
};

extern FakeRatesDef *global_fakerates;

} // namespace nuX_FakeRates
