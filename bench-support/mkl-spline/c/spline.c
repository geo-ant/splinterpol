#include <mkl_df.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>

struct MklInterpolator_1d_f64 {
  // internal MKL task pointer
  DFTaskPtr task;
  // internal for the range [0,Ndata-1]
  double *range; // array len 2
  // internal for the boundary conditions (zero first derivs)
  double *bc; // array len 2
};

// free an interpolator
void MklInterpolator_1d_f64_free(struct MklInterpolator_1d_f64 *interp) {
  if (interp->task) {
    dfDeleteTask(&interp->task);
    interp->task = NULL;
  }
  if (interp->bc) {
    free(interp->bc);
    interp->bc = NULL;
  }
  if (interp->range) {
    free(interp->range);
    interp->range = NULL;
  }
}

// create a fitting task handle to fit a number of `ndata` data points of type
// float 64.
// Returns false on error.
//
// NOTE: this creates a uniform partition for 1D fitting.
bool MklInterpolator_1d_f64_init(
    struct MklInterpolator_1d_f64 *interp,
    // must be valid for the lifetime of the interpolator
    double const *const in_data, MKL_INT ndata,
    // output data buffer for the coefficients and
    // it's length. Must have at least length
    // 4*(ndata-1).
    double *out_coeffs, size_t const coeffs_len) {

  if (interp->task != NULL || interp->range != NULL || interp->bc != NULL) {
    return false;
  }

  // just to be safe
  if (ndata < 4) {
    return false;
  }

  // coefficients len, see:
  // https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2026-2/df-editppspline1d.html
  if (coeffs_len < ((size_t)DF_PP_CUBIC) * ((size_t)(ndata - 1))) {
    return false;
  }

  interp->range = (double *)malloc(2 * sizeof(double));
  if (interp->range == NULL) {
    return false;
  }
  interp->range[0] = 0.;
  interp->range[1] = (double)(ndata - 1);

  if (DF_STATUS_OK != dfdNewTask1D(&interp->task, ndata, interp->range,
                                   DF_UNIFORM_PARTITION, 1, in_data,
                                   DF_NO_HINT)) {
    MklInterpolator_1d_f64_free(interp);
    return false;
  }

  // setting the left and right boundary conditions to zero SHOULD produce
  // the same spline interpolation as whole symmetric, if I understand
  // correctly.
  MKL_INT bc_type = DF_BC_1ST_LEFT_DER | DF_BC_1ST_RIGHT_DER;

  interp->bc = (double *)malloc(2 * sizeof(double));
  if (interp->bc == NULL) {
    MklInterpolator_1d_f64_free(interp);
    return false;
  }
  interp->bc[0] = 0.;
  interp->bc[1] = 0.;

  if (DF_STATUS_OK !=
      dfdEditPPSpline1D(
          interp->task, DF_PP_CUBIC,
          // see
          // https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2026-2/df-editppspline1d.html
          // >  NOTE: To construct a natural cubic spline, pass these
          // > settings to the editor: DF_PP_CUBIC as the spline order,
          // > DF_PP_NATURAL as the spline type, and DF_BC_FREE_END
          // > as the boundary condition.
          // > To construct a cubic spline with other boundary conditions,
          // > pass these settings to the editor: DF_PP_CUBIC as the spline
          // order, > DF_PP_NATURAL as the spline type, and the required type of
          // > boundary condition.
          //
          // which means passing NATURAL here doesn't force natural
          // boundary conditions.
          DF_PP_NATURAL, bc_type, interp->bc,
          // No internal conditions. The docs in the table say ic has
          // to be non-null, but here:
          // https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2026-2/df-editppspline1d.html
          // on that page in the prose it says to pass NULL when
          // no boundary or internal conditions are given
          DF_NO_IC, NULL,
          // confusing signature, this doesn't actually write into the buffer
          // but accepts the pointer to it. The data must be writable for
          // the coefficient calculation to work later.
          out_coeffs, DF_NO_HINT)) {
    MklInterpolator_1d_f64_free(interp);
    return false;
  }

  return true;
}

// perform the coefficient calculation with the input and output arrays
// given at construction.
//
// NOTE: I think this should perform a recalculation, even if nothing changed.
// If the values seem spuriously small, I'll have to check my assumptions.
bool MklInterpolator_1d_f64_calculate_coefficients(
    struct MklInterpolator_1d_f64 interp) {
  return DF_STATUS_OK ==
         dfdConstruct1D(interp.task, DF_PP_SPLINE, DF_METHOD_STD);
}

// perform interpolation at multiple sites. No special ordering for the sites
// is assumed.
//
// NOTE: for non-random resampling, there are actually hints and ordering
// parameters that we can pass to speed up the computation potentially!
bool MklInterpolator_1d_f64_interp_multi_generic(
    struct MklInterpolator_1d_f64 interp, double const *const sites,
    MKL_INT nsites, double *const out_values) {

  // see: https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2026-2/df-interpolate1d-df-interpolateex1d.html
  MKL_INT dorder = 1; // request function values (derivative order zero)
  return DF_STATUS_OK == dfdInterpolate1D(interp.task, DF_INTERP, DF_METHOD_PP,
                                          nsites, sites, DF_NO_HINT, 1, &dorder,
                                          NULL, out_values, DF_NO_HINT, NULL);
}
