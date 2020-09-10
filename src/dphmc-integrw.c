# include "dphmc-integrw.h"

# include <gsl/gsl_monte_vegas.h>
# include <gsl/gsl_integration.h>

# include "dphmc-inspection.h"

/** Performs repeating integration using Gauss-Kronrod method from GSL,
 * mitigating relative error requirement at each failed step.
 *
 * Temporary disables native GSL error handler to prevent abortion of the
 * execution and re-iterate the QAGS procedure with relative error mitigated
 * by factor `epsrelIncFt`.
 *
 * @param qagsp A fully-initialized instance of dphmc_IterativeQAGSParameters
 * @param f Integrand function
 * @param ps Parameters set for integrand function
 * @param xLow lower limit for integration
 * @param xUp upper limit of integration
 * @param relError pointer to `double` variable where the relative error is to
 *        be written
 * @param qagspWS pointer to allocated `gsl_integration_workspace`, valid for
 * amount of nodes specified in `qagsp` structure. May be NULL.
 * @returns integrap estimation with relative error written by `relErrPtr`
 * */
double
dphmc_QAGS_integrate_iteratively( const struct dphmc_IterativeQAGSParameters * qagsp
                                , double (*f)(double, void*), const void * ps
                                , double xLow, double xUp
                                , double * relErrPtr, double * absErrPtr
                                , void * qagspWS_
                                ) {
    double res = NAN;
    int rr;
    gsl_integration_workspace * intWS = ( qagspWS_ ? qagspWS_
                                        : gsl_integration_workspace_alloc( qagsp->nnodes ) );
    gsl_function F = {
            .function = f,
            .params = (void *) ps
        };

    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();

    for( *relErrPtr = qagsp->epsrel;
         *relErrPtr < 1.;
         *relErrPtr *= qagsp->epsrelIncFt ) {
        if( GSL_EROUND == (rr = gsl_integration_qags( &F
                                                    , xLow, xUp
                                                    , qagsp->epsabs, *relErrPtr
                                                    , qagsp->limit, intWS
                                                    , &res, absErrPtr)) ) {
            DPhMC_msg( "Integration: at gsl_integration_qags(...) needs to be"
                    " more tolerant [%e, %e]; changing epsrel %e -> %e."
                    , xLow, xUp, *relErrPtr, (*relErrPtr) * qagsp->epsrelIncFt );
            continue;
        } else if( !rr ) {
            break;  /* integrated correctly */
        } else {
            DPhMC_msge( "Lower integration limit: %e (value %e)", xLow, f(xLow, ps) );
            DPhMC_msge( "Upper integration limit: %e (value %e)", xUp, f(xUp, ps) );
            DPhMC_msge( "Current relative error: %e", *relErrPtr );
            DPhMC_error( DPHMC_THIRD_PARTY_ERROR
                       , "Got a gsl_integration_qags(...) error %d: \"%s\""
                       , rr
                       , (gsl_strerror(rr) ? gsl_strerror(rr) : "<unknown>" )
                       );
        }
    }
    gsl_set_error_handler(old_handler);
    if( !isfinite(res) ) {
        DPhMC_error( DPHMC_INTEGRATION_ERROR
                   , "dphmc_aprime_chi() integration evaluated to infinite"
                   " value: %e."
                   , res );
    }
    if( ! qagspWS_) {
        gsl_integration_workspace_free( intWS );
    }
    return res;
}



