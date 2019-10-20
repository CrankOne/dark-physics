# include "dphmc-config.h"

# include <assert.h>
# include <stdlib.h>
# include <stdio.h>
# include <stdint.h>
# include <math.h>
# include <unistd.h>
# include <time.h> 
# include <sys/types.h>
# include <sys/stat.h>
# include <string.h>
# include <strings.h>

# include <gsl/gsl_math.h>
# include <gsl/gsl_multimin.h>
# include <gsl/gsl_monte_vegas.h>
# include <gsl/gsl_integration.h>
# include <gsl/gsl_rng.h>

# include "models/ww-approx.h"

# define me     _DPhMC_CST_electronMass_GeV
# define alpha  _DPhMC_CST_fineStructureConstant
# define muP    _DPhMC_CST_muProton
# define mp     _DPhMC_CST_protonMass_GeV

/* This type is local C-structure, hidden for other units. */
struct Workspace {
    struct dphmc_APrimeWSParameters parameterSet;
    double chi, tRange[2], chiRelError, chiAbsError;
    gsl_integration_workspace * gslIntWS;

    double xRange[2],
           thetaMax,
           integralEstVal ;  /**< according to (A16) */

    /** Reentrant VEGAS integration workspace, used for numerical integration. */
    gsl_monte_vegas_state * gslVegasSt;
};

struct Workspace * /* blank workspace ctr */
_static_allocate_workspace() {
    struct Workspace * p = malloc( sizeof(struct Workspace) );
    bzero( p, sizeof(struct Workspace) );
    p->tRange[0] = p->tRange[1] = NAN;
    p->chi = NAN;
    p->gslIntWS = NULL;
    p->gslVegasSt = NULL;
    return p;
}

void /* workspace dtr */
_static_free_workspace( struct Workspace * w ) {
    if( w ) {
        if( w->gslIntWS ) {
            gsl_integration_workspace_free(w->gslIntWS);
        }
        if( w->gslVegasSt ) {
            gsl_monte_vegas_free( w->gslVegasSt );
        }
        free( w );
    }
}

static double
_static_calc_analytic_integral_estimation( void * wsPtr ) {
    struct Workspace * ws = (struct Workspace *) wsPtr;
    const struct dphmc_APrimeWSParameters * ps = &(ws->parameterSet);
    const double AM = ps->massA_GeV
               , E0 = ps->EBeam_GeV
               ;
    double c1 = me/ps->massA_GeV
         , c2 = ps->massA_GeV / ps->EBeam_GeV  /*ws->xRange[0]*/
         , upLim = c1 > c2 ? c1 : c2;
    upLim *= upLim;
    /*upLim = c2;  // NOTE: NA64 tooks it like that (no ^2, no max, just m_A/E_0)*/
    /* See (A15) and footnote 2 on page 2 of arXiv 1209.6083 */
    /* Note: I tried to perform similar integration procedure and got factor 6
     * instead of 3 in denominator after dx integration, so predicted output
     * should be ~4 times lesser than (A15) predicts. */
    return /*(8./3)*/ /*(2./3)*/ (4./3)
          * pow((alpha*ps->epsilon)/AM, 2) * alpha * sqrt( 1 - pow(AM/E0, 2) )
          * ws->chi
          * log( 1 / upLim )
          * (ps->factor)
          * _DPhMC_CST_pbarn_per_GeV
          ;
}

void
dphmc_print_aprime_ws_parameters( FILE * f
                                , const struct dphmc_APrimeWSParameters * p ) {
    fprintf( f, "{Z=%d; A=%f; m_{A'}=%e; E_0=%e; \\eps=%e; factor=%e; "
              , p->Z, p->A, p->massA_GeV, p->EBeam_GeV, p->epsilon, p->factor );
    fprintf( f, "epsabs=%e; epsrel=%e; epsrelIncFt=%e; limit=%zu; nnodes=%zu}"
              , p->epsabs, p->epsrel, p->epsrelIncFt, p->limit, p->nnodes );
}

/**The physical parameters should be provided by pointer on corresponding
 * const structure. This data is then copied to the "workspace" object that
 * stores also some intermediate caches needed to compute the cross-section
 * values.
 *
 * @param n nodes for GSL integration workspace (see GSL QAGS);
 * @returns -1 on error, 0 on success */
int
dphmc_init_aprime_cs_workspace( const struct dphmc_APrimeWSParameters * newPs
                              , void ** wsPtr_ ) {
    /* Local enum, storing evaluation instructions */
    const int doRecalculation             = 0x1
            , doGSLWorkspaceAllocation    = 0x2 | 0x1
            , doParametersReassignment    = 0x4 | 0x1
            ;
    int routineFlags = 0x0;
    struct Workspace ** wsPtr = (struct Workspace **) wsPtr_
                   ,  * ws = NULL;
    struct dphmc_APrimeWSParameters * ps = NULL;
    /* decide, what to do: initialize workspace and parameter set pointers */
    if( *wsPtr ) {  /* parameter set already allocated */
        if( memcmp( &((*wsPtr)->parameterSet)
                  , newPs
                  , sizeof(struct dphmc_APrimeWSParameters)
                  ) ) {  /* parameter set has changed */
            if( (*wsPtr)->parameterSet.nnodes ) {
                /* GSL nodes changed -- need to re-allocate GSL workspace) */
                gsl_integration_workspace_free( (*wsPtr)->gslIntWS );
                (*wsPtr)->gslIntWS = NULL;
                routineFlags |= doGSLWorkspaceAllocation;
            }
            routineFlags |= doRecalculation;
            routineFlags |= doParametersReassignment;
        }
    } else {
        /* parameters were not allocated */
        *wsPtr = _static_allocate_workspace();
        routineFlags |= doRecalculation;
        routineFlags |= doParametersReassignment;
    }

    ws = *wsPtr;
    ps = &(ws->parameterSet);

    if( routineFlags & doParametersReassignment ) {
        memcpy( ps
              , newPs
              , sizeof(struct dphmc_APrimeWSParameters)
              );
    }

    if( routineFlags & doGSLWorkspaceAllocation ) {
        ws->gslIntWS = gsl_integration_workspace_alloc( ps->nnodes );
    }

    if( routineFlags & doRecalculation ) {
        /*DPhMC_msg( "Recalculation of chi form-factor invoked for E0=%e.", ps->EBeam_GeV );*/
        ws->tRange[0] = pow(ps->massA_GeV, 4)/(4*(ps->EBeam_GeV)*(ps->EBeam_GeV));
        /* TODO: ^^^ pow(U/(2*E0*(1-x)), 2); for excellent approximation */
        ws->tRange[1] = ps->massA_GeV*ps->massA_GeV;
        ws->chi = dphmc_aprime_chi( &(ws->chiAbsError), &(ws->chiRelError), ws );
        if(!(isfinite(ws->chi) &&
             isfinite(ws->tRange[0]) &&
             isfinite(ws->tRange[1])) ) {
            DPhMC_error( DPHMC_BAD_VALUE
                       , "Bad value of t_min=%e/t_max=%e/chi=%e."
                       , ws->tRange[0], ws->tRange[1], ws->chi );
        }
        /* Set the theta integration limits.
         * theta_up -- according to formula (9) of \cite{Bjorken} */
        {
            double val1 = sqrt( ps->massA_GeV*_DPhMC_CST_electronMass_GeV )/ps->EBeam_GeV,
                   val2 = pow( ps->massA_GeV/ps->EBeam_GeV, 1.5 );
            ws->thetaMax = (val1 > val2 ? val1 : val2);
        }
        /* Set the x integration limits.
         *  x_low -- mass threshold
         *  x_up -- given by singularity conditions */
        ws->xRange[0] = ps->massA_GeV / ps->EBeam_GeV;
        {   /* from condition (A16) */
            double c1 = me/ps->massA_GeV,
                   c2 = ws->xRange[0];
            ws->xRange[1] = c1 > c2 ? c1 : c2;
            ws->xRange[1] = 1 - ws->xRange[1]*ws->xRange[1];
        }
        if( ws->xRange[0] > 1 || ws->xRange[0] < 0. ) {
            DPhMC_msgw( "x_min > 1 || x_min < 0 for the A' integration workspace %p"
                         "(E_0 = %e GeV, m_A'= %e GeV, x_{min} = %e)!",
                         ws, ps->EBeam_GeV, ps->massA_GeV, ws->xRange[0] );
        }
        if( ws->xRange[1] > 1 || ws->xRange[1] < 0. ) {
            DPhMC_msgw( "x_max > 1 || x_max < 0 for the A' integration workspace %p"
                         "(E_0 = %e GeV, m_A'= %e GeV, x_{max} = %e)!",
                         ws, ps->EBeam_GeV, ps->massA_GeV, ws->xRange[1] );
        }
        if( ws->thetaMax > M_PI ) {
            DPhMC_msge(
                         "theta_max = %e >= pi for the A' integration workspace %p"
                         " (E_0 = %e GeV, m_A'= %e GeV, x_{min} = %e, x_{max} = %e)."
                         " Correcting it to pi."
                       , ws->thetaMax
                       , ws, ps->EBeam_GeV, ps->massA_GeV, ws->xRange[0], ws->xRange[1] );
            ws->thetaMax = M_PI;
        }
        if( ! (ws->xRange[0] < ws->xRange[1]) ) {
            DPhMC_msge(
                         "x_min >= x_max for the A' integration workspace %p"
                         " (E_0 = %e GeV, m_A'= %e GeV, x_{min} = %e, x_{max} = %e)!"
                       , ws, ps->EBeam_GeV, ps->massA_GeV, ws->xRange[0], ws->xRange[1] );
        }
        ws->integralEstVal = _static_calc_analytic_integral_estimation( ws );
    }
    assert( ps->A > 0 );
    assert( ps->Z > 0 );
    return 0;
}

void
dphmc_free_aprime_cs_workspace( void * ws_ ) {
    if( !ws_ ) return;
    _static_free_workspace( (struct Workspace *) ws_ );
}

/** For given nuclei parameters, calculates integrand function as a summation
 * of elastic/inelastic components.
 * see \ref{JDBjorken} (A18,19). */
static double
chi_integrand(double t, void * ws_) {
    struct Workspace * ws = (struct Workspace *) ws_;
    assert( ws );
    assert( ws->parameterSet.A > 0 );
    assert( ws->parameterSet.Z > 0 );
    const double A = ws->parameterSet.A,
                 Z = ws->parameterSet.Z,
                 d = 0.164/pow(A, 2./3),
                ap = 773/(pow(Z, 2./3)*me),
                 a = 111/(pow(Z, 1./3)*me),
                a2 = a*a,
                t2 = t*t,
               ap2 = ap*ap,
              G2el = pow( a2*t*Z/( (1 + a2*t)*(1 + t/d) ), 2),
              G2in = Z*pow( ( (ap2*t)/(1+ap2*t) )*( 1+t*(muP*muP-1)/(4*mp*mp) ) /
                                    pow(1 + t/.71, 4)
                  , 2 ),
                G2 = G2el + G2in,
               res = G2*(t - ws->tRange[0])/t2
                   ;
    if( !isfinite(res) ) {
        DPhMC_msge( "A=%e, Z=%e, d=%e, a'=%e, G2el=%e, G2in=%e\n", A, Z, d, ap, G2el, G2in );
        DPhMC_error( DPHMC_BAD_VALUE
                   , "chi_integrand() evaluated to bad value: %e."
                   , res );
    }
    return res;
}


/**For given nuclei parameters, calculates integrated photons flux
 * involving GSL integration procedure.
 *
 * see \ref{JDBjorken} (A17). */
double
dphmc_aprime_chi( double * absError
                , double * relError
                , void * ws_ ) {
    assert(ws_);
    struct Workspace * ws = (struct Workspace *) ws_;
    assert( ws->parameterSet.A > 0 );
    assert( ws->parameterSet.Z > 0 );
    const struct dphmc_APrimeWSParameters * ps = &(ws->parameterSet);

    double chi_ = NAN; {
        gsl_function F = {
                .function = chi_integrand,
                .params = ws_
            };
        int rr;

        gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
        for( *relError = ps->epsrel;
             *relError < 1.;
             *relError *= ps->epsrelIncFt ) {
            if( GSL_EROUND == (rr = gsl_integration_qags(
                    &F,
                    ws->tRange[0],
                    ws->tRange[1],
                    ps->epsabs,
                    *relError,
                    ps->limit,
                    ws->gslIntWS,
                    &chi_,
                    absError)) ) {
                DPhMC_msg( "chi parameter integration: at"
                        " gsl_integration_qags(...) needs to be more tolerant"
                        " [%e, %e]; changing epsrel %e -> %e.",
                                ws->tRange[0],
                                ws->tRange[1],
                                *relError,
                                (*relError) * ps->epsrelIncFt);
                continue;
            } else if( !rr ) {
                break;  /* integrated correctly */
            } else {
                DPhMC_msge( "Lower integration limit: %e", ws->tRange[0] );
                DPhMC_msge( "Upper integration limit: %e", ws->tRange[1] );
                DPhMC_msge( "Current relative error: %e", *relError );
                DPhMC_error( DPHMC_THIRD_PARTY_ERROR
                           , "Got a gsl_integration_qags(...) error %d: \"%s\""
                           , rr
                           , (gsl_strerror(rr) ? gsl_strerror(rr) : "<unknown>" )
                           );
            }
        }
        gsl_set_error_handler(old_handler);
    }
    if( !isfinite(chi_) ) {
        DPhMC_error( DPHMC_INTEGRATION_ERROR
                   , "dphmc_aprime_chi() integration evaluated to infinite"
                   " value: %e."
                   , chi_ );
    }

    return chi_;
}

/**Represents volatile part of function \f$d \d sigma / d x\f$
 * corrsponding to non-normalized PDF, that may be of use in generators that
 * does not require angular distribution. */
double
dphmc_aprime_unnormed_pdf_a12( double x
                             , double theta
                             , void * ws ) {
    assert(ws);
    const struct dphmc_APrimeWSParameters * ps
            = &(((const struct Workspace *) ws)->parameterSet);
    const double AM = ps->massA_GeV
               , E0 = ps->EBeam_GeV
               , am2 = AM*AM
               , U = pow(E0*theta, 2)*x + am2*(1 - x)/x + me*me*x
               , U2 = U*U
               , factor1 = (1 - x + x*x/2)/U2
               , factor2 = pow( (1 - x)*AM/U2, 2)
               , factor3 = am2 - (U*x)/(1 - x)
               , res = (factor1 + factor2*factor3)*x
               ;
    assert(isfinite(res));
    return res;
}

/**Calculates sigma according to (A12) formula \cite{JDBjorken}, returning
 * result of numerical calculation of the following value:
 * \f[
 * \frac{1}{E_{0}^2 x} \frac{ d \sigma_{ 2 \rightarrow 3 } }{ d x d \cos{\theta_{A'}} }
 * \f]
 *
 * This relation includes calculation of constant term that is not needed for
 * simulation purposes. For the event generation purposes, se the "unnormed
 * PDF" function dphmc_aprime_unnormed_pdf_a12() instead. */
double
dphmc_aprime_cross_section_a12( double x
                              , double theta
                              , void * ws_ ) {
    assert(ws_);
    const struct Workspace * ws = (const struct Workspace *) ws_;
    const struct dphmc_APrimeWSParameters * ps = &(ws->parameterSet);
    const double AM = ps->massA_GeV
               , E0 = ps->EBeam_GeV
               , betaAPrime = sqrt( 1 - pow(AM/E0, 2) )
               , factor4 = 8*alpha*pow(2*alpha*ps->epsilon, 2)*betaAPrime*(ws->chi)
               , res = factor4 * dphmc_aprime_unnormed_pdf_a12(x, theta, ws_)
                               * (E0 * E0 * x) * _DPhMC_CST_pbarn_per_GeV * sin(theta)
               ;
    return res * (ps->factor);
}

/**Represents function \f$d \d sigma / d x\f$
 * corrsponding to dfferential by emitted A' energy cross section, that may be
 * useful in generators that does not require angular distribution. */
double
dphmc_aprime_cross_section_a14( double x
                              , void * ws_ ) {
    assert(ws_);
    const struct Workspace * ws = (const struct Workspace *) ws_;
    const struct dphmc_APrimeWSParameters * ps = &(ws->parameterSet);
    const double AM = ps->massA_GeV
               , E0 = ps->EBeam_GeV
               , betaAPrime = sqrt( 1 - pow(AM/E0, 2) )
               , factor4 = 4*alpha*pow(alpha*ps->epsilon, 2)*betaAPrime*(ws->chi)
               , res = factor4 * ( 1 - x + x*x/3. ) / ( AM*AM*(1 - x)/x + me*me*x )
                               * _DPhMC_CST_pbarn_per_GeV
               ;
    return res * (ps->factor);
}

double
dphmc_aprime_lower_cut_x( void * wsPtr ) {
    struct Workspace * ws = (struct Workspace *) wsPtr;
    return ws->xRange[0];
}

double
dphmc_aprime_upper_cut_x( void * wsPtr ) {
    struct Workspace * ws = (struct Workspace *) wsPtr;
    return ws->xRange[1];
}

double
dphmc_aprime_upper_cut_theta( void * wsPtr ) {
    struct Workspace * ws = (struct Workspace *) wsPtr;
    return ws->thetaMax;
}

double
dphmc_aprime_ww_fast_integral_estimation( void * wsPtr ) {
    struct Workspace * ws = (struct Workspace *) wsPtr;
    return ws->integralEstVal;
}

static double
_static_ww_aprime_cs_wrapper( double * x, size_t dim, void * params ) {
    assert( 2 == dim );
    return dphmc_aprime_cross_section_a12( x[0], x[1], params );
}

/**Performs Gauss-Kronrod integration procedure of
 * dphmc_aprime_cross_section_a14(). Faster then dphmc_aprime_ww_full_numeric_2(),
 * but yields a less-precise value.
 * @todo: is still inefficient for MC generator. Consider to move it to
 * testing */
double
dphmc_aprime_ww_full_numeric_1( void * ws_
                              , double epsabs
                              , double epsrel
                              , double epsrelIncFt
                              , size_t limit, size_t nnodes
                              , double * relError ) {
    assert(ws_);
    struct Workspace * ws = (struct Workspace *) ws_;

    double sigma = NAN, absError = 0;
    gsl_integration_workspace * gslIntWSSigma 
        = gsl_integration_workspace_alloc( nnodes );
    {
        gsl_function F = {
                .function = dphmc_aprime_cross_section_a14,
                .params = ws_
            };
        int rr;

        gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
        for( *relError = epsrel;
             *relError < 1.;
             *relError *= epsrelIncFt ) {
            if( GSL_EROUND == (rr = gsl_integration_qags(
                    &F,
                    ws->xRange[0], ws->xRange[1],
                    epsabs,
                    *relError,
                    limit,
                    gslIntWSSigma,
                    &sigma,
                    &absError)) ) {
                DPhMC_msg( "chi parameter integration: at"
                        " gsl_integration_qags(...) needs to be more tolerant"
                        " [%e, %e]; changing epsrel %e -> %e.",
                                ws->tRange[0], ws->tRange[1],
                                *relError,
                                (*relError) * epsrelIncFt);
                continue;
            } else if( !rr ) {
                break;  /* integrated correctly */
            } else {
                DPhMC_msge( "Lower integration limit ... %e", ws->xRange[0] );
                DPhMC_msge( "Upper integration limit ... %e", ws->xRange[1] );
                DPhMC_msge( "Current rel error ... %e\n", *relError );
                DPhMC_error( DPHMC_THIRD_PARTY_ERROR
                           , "Got a gsl_integration_qags(...) error %d: \"%s\""
                           , rr
                           , (gsl_strerror(rr) ? gsl_strerror(rr) : "<unknown>" )
                           );
            }
        }
        gsl_set_error_handler(old_handler);
    }
    gsl_integration_workspace_free( gslIntWSSigma );
    if( !isfinite(sigma) ) {
        DPhMC_error( DPHMC_INTEGRATION_ERROR
                   , " %s integration evaluated to infinite"
                   " value: %e.", __FUNCTION__
                   , sigma );
    }

    return sigma;
}

/**Performs VEGAS two-folded integration of
 * dphmc_aprime_unnormed_pdf_a12().
 * @todo: may be precise, but is extemely inefficient. Move it to testing? */
double
dphmc_aprime_ww_full_numeric_2( void * ws_
                              , size_t nCalls
                              , double * absError
                              , double * chi2 ) {
    struct Workspace * ws = (struct Workspace *) ws_;
    if( ! ws->gslVegasSt ) {
        fflush( stdout );
        ws->gslVegasSt = gsl_monte_vegas_alloc(2);
        int rc = gsl_monte_vegas_init( ws->gslVegasSt );
        assert( ! rc );
    }
    # if 1
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc( T );
    # endif
    gsl_monte_function F = { &_static_ww_aprime_cs_wrapper, 2, ws_ };
    double ranges[] = { ws->xRange[0], 0
                      , ws->xRange[1], ws->thetaMax
                      };
    double result, absError_;
    gsl_monte_vegas_integrate( &F
                             , ranges
                             , ranges + 2
                             , 2
                             , nCalls
                             , r
                             , ws->gslVegasSt
                             , &result
                             , &absError_ );
    if( absError ) *absError = absError_;
    if( chi2 ) *chi2 = gsl_monte_vegas_chisq(ws->gslVegasSt);
    return result;
}

/**Used mainly by accept-reject (direct von Neumann) sampling method */
double dphmc_aprime_ww_upper_limit( void * ws ) {
    # warning "Here be dragons"
}
