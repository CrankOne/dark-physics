# ifndef H_DPHMC_APRIME_WW_INTEGRATION_WRAPPERS_H
# define H_DPHMC_APRIME_WW_INTEGRATION_WRAPPERS_H

# include <stdlib.h>

# ifdef __cplusplus
extern "C" {
# endif

/**@brief A parameters set for integrating with GSL Gauss-Kronrod algoritm.
 *
 * @details This integration algorith is used in multiple places for
 * integrating 1D functions in this library. Iterative procedure over the GSL
 * function is applied to gain numerical estimation, gradually loosing the
 * requirement on the relative error in case of convergence failures. */
struct dphmc_IterativeQAGSParameters {
    double epsabs  /**< absolute error */
         , epsrel  /**< relative error */
         , epsrelIncFt  /**< relative error increment factor */
         ;
    size_t limit  /**< maximum number of subintervals (see GSL QAGS) */
         , nnodes  /**< GSL integration workspace nodes num */
         ;
};

/**An iterative wrapper around the GSL QAGS integration algorithm.*/
double
dphmc_QAGS_integrate_iteratively( const struct dphmc_IterativeQAGSParameters * qagsp
                                , double (*f)(double, void*), const void * ps
                                , double xLow, double xUp
                                , double * relErrPtr, double * absErrPtr
                                , void * qagspWS );

# ifdef __cplusplus
}
# endif

# ifdef __cplusplus
namespace dphmc {
typedef dphmc_IterativeQAGSParameters IterativeQAGSParameters;
}
# endif

# endif  /* H_DPHMC_APRIME_WW_APPROX_MAJORANT_H */

