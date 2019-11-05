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
# include <gsl/gsl_deriv.h>
# include <gsl/gsl_roots.h>

# include "models/ww-approx.h"

# define me     _DPhMC_CONST_electronMass_GeV
# define alpha  _DPhMC_CONST_fineStructureConstant
# define muP    _DPhMC_CONST_muProton
# define mp     _DPhMC_CONST_protonMass_GeV


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
 * be written
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
            DPhMC_msge( "Lower integration limit: %e", xLow );
            DPhMC_msge( "Upper integration limit: %e", xUp );
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

struct dphmc_ChiCache {
    double chi  /* constant approximation of \chi */
         , tRange[2]  /*< virtuality range to integrate over, for \chi */
         , chiRelError  /*< relative error of \chi */
         , chiAbsError  /*< absolute error of \chi XXX? */
         ;
    struct dphmc_APrimeFormFactor ff;
};

/* Caches for CS-calculus procedures, local for this object file.
 * This type is local C-structure, hidden for other units. */
struct dphmc_APrimeWWCaches {
    /* Physical parameters used for this caches */
    struct dphmc_APrimePhysParameters physPars;

    /* These values are derived from `physPars` and remains constant */
    double xRange[2]  /*< x-range, for integration */
         , thetaMax  /*< maximum theta value, for integration */
         , integralEstVal  /*< according to (A16) */
         , cstFact  /*< constant term of the (A12): 8*\alpha^3*\eps^2\geta_{A'} */
         ;

    /* GSL QAGS parameters used to integrate \chi */
    struct dphmc_IterativeQAGSParameters chiIntPars;
    /* Stores result of last calculation of \chi. This data will be set
     * unconditionally, even if \chi(\theta) calculus is enabled. */
    struct dphmc_ChiCache constChiCache;
    /* Reentrant QAGS integration workspace used for calculation of \chi. If
     * set to NULL, the re-calculation is disabled */
    gsl_integration_workspace * gslChiIntWS;

    /* Reentrant VEGAS integration workspace, used for numerical integration. */
    gsl_monte_vegas_state * gslVegasSt;
};

void
dphmc_print_aprime_ws_parameters( FILE * f
                                , const struct dphmc_APrimePhysParameters * p ) {
    fprintf( f, "{Z=%d; A=%f; m_{A'}=%e; E_0=%e; \\eps=%e;"
                " factor=%e; maxThetaFactor=%e}"
              , p->Z, p->A, p->massA_GeV, p->EBeam_GeV, p->epsilon
              , p->factor, p->maxThetaFactor );
}

/**Certain basic check will be performed to assure the approximation is appliable.
 *
 * @returns 0 if values are acceptable.
 * @return -1 if \f$A'\f$ mass is too small
 * @return -2 \f$E_0\f$ is comparable to \f$m_{A'}\f$
 * */
int
dphmc_aprime_ww_check_phys_parameters_are_valid( const struct dphmc_APrimePhysParameters * ps ) {
    if( ps->massA_GeV <= 4*me ) {
        return -1;  /* A' mass is too small */
    }
    if( ps->massA_GeV > ps->EBeam_GeV/10 ) {
        return -2;  /* A' mass is of the same order as E0 */
    }
    /* ... TODO: more */
    return 0;
}

/**The physical parameters should be provided by pointer on corresponding
 * const structure. This data is then copied to the "workspace" object that
 * stores also some intermediate caches needed to compute the cross-section
 * values.
 *
 * The following bits in `flags` attribute of the structure have meaning:
 *  0x1 -- takes into accunt \f$\chi\f$ dependency of \f$\theta\f$ by
 *         forcing numerical integration of \f$\chi\f$. May dramatically affect
 *         the performance, the effect is expected to be negligible
 *  0x2 -- Sets the upper \f$\theta_{max}\f$ cut-off angle to \f$\pi\f$. The
 *         \f$\theta_{max}\f$ cut-off angle is expected to be proportional
 *         to \f$ \max(\frac{\sqrt{m_{A'}m_e}}{E_0}, (m_{A'}/E_0)^{3/2}), with
 *         coefficiont \f$~ 1\f$. Enabling this flag turns integration on for
 *         the full \f$[0,\pi]\f$ range of \f$\theta\f$. One may control the
 *         coefficient of proprotionality with `maxThetaFactor` attribute with
 *         this flag being disabled, or set it to \$\pi\$ by enabling this
 *         flag.
 *  0x4 -- makes routine to omit the physics parameters validity check. May
 *         cause incorrect results, floating point errors, but may be of use
 *         for some sloppy approximations / preliminary estimations for
 *         inheriting models, etc.
 *
 * @param ps physics parameters; copied into caches
 * @param chiIntPars \f$\chi\f$ integration parameters, copied into caches
 * */
struct dphmc_APrimeWWCaches *
dphmc_aprime_new( const struct dphmc_APrimePhysParameters * ps
                , const struct dphmc_IterativeQAGSParameters * chiIntPars
                , double (*  elastic_f)( double, struct dphmc_APrimePhysParameters * )
                , double (*inelastic_f)( double, struct dphmc_APrimePhysParameters * )
                , int flags ) {
    struct dphmc_APrimeWWCaches * caches = malloc(sizeof(struct dphmc_APrimeWWCaches));
    /* Copy the parameters */
    memcpy( &(caches->physPars)
          , ps
          , sizeof(struct dphmc_APrimePhysParameters) );
    if( (!(0x4 & flags))
      && dphmc_aprime_ww_check_phys_parameters_are_valid(&(caches->physPars)) ) {
        DPhMC_msge( "Invalid physical parameters. Refuse to allocate workspace." );
        free( caches );
        return NULL;
    }
    /* Compute some values derived directly from phys. pars */
    {
        /* lower x bound is determined from common sense -- mass threshold */
        caches->xRange[0] = caches->physPars.massA_GeV / caches->physPars.EBeam_GeV;
        /* upper x bound is defined by divergence of the cross section, -- see
         * condition (A16) in \cite Bjorken */
        {
            double c1 = me/ps->massA_GeV
                 , c2 = caches->xRange[0];
            caches->xRange[1] = c1 > c2 ? c1 : c2;
            caches->xRange[1] = 1 - caches->xRange[1] /*caches->xRange[1]*/;  // TODO: squared?
            assert( caches->xRange[1] < 1. );
            assert( caches->xRange[1] > caches->xRange[0] );
        }
        if( 0x2 & flags ) {
            /* flags forces the theta upper limit to be set strictly to \pi. */
            caches->thetaMax = M_PI;
        } else {
            double val1 = sqrt( caches->physPars.massA_GeV * me )/caches->physPars.EBeam_GeV
                 , val2 = pow( caches->physPars.massA_GeV / caches->physPars.EBeam_GeV, 1.5 );
            caches->thetaMax = caches->physPars.maxThetaFactor * (val1 > val2 ? val1 : val2);
            if( caches->thetaMax > M_PI ) {
                /* since user can define their own proportionality factor for
                 * theta max, this check is not redundant, though may indicate
                 * user's mistake; todo: consider put warning here? */
                caches->thetaMax = M_PI;
            }
        }
        /* pre-computing this factor may save some CPU time for future as it
         * is used at many places here. Note the absence of \chi here --
         * depending on user's flag we, contrary to original article, may not
         * neglect the dependance of \theta */
        {
            const double AM = caches->physPars.massA_GeV
                       , E0 = caches->physPars.EBeam_GeV
                       , betaAPrime = sqrt( 1 - pow(AM/E0, 2) )
                       ;
            caches->cstFact = 4*alpha*pow(alpha * caches->physPars.epsilon, 2)
                            * betaAPrime
                            * _DPhMC_CONST_pbarn_per_GeV
                            * caches->physPars.factor
                            ;
        }
        /* Set the full sigma approximation given by formula (A15) to NAN. Once
         * requested by `dphmc_aprime_ww_fast_integral_estimation()` will be
         * computed and cached for further usage */
        caches->integralEstVal = nan("lazy");
    }
    /* Copy \chi parameters and allocate the integration workspace for further
     * usage */
    memcpy( &(caches->chiIntPars)
          , chiIntPars
          , sizeof(struct dphmc_APrimeWWCaches) );
    caches->gslChiIntWS = gsl_integration_workspace_alloc( caches->chiIntPars.nnodes );
    /* Compute the constant chi approximation. Even if \chi dependency of
     * \theta has to be taken into account, there is some approximations that
     * rely on this value */
    {
        const double ma = caches->physPars.massA_GeV
                   , E0 = caches->physPars.EBeam_GeV
                   ;
        struct dphmc_ChiCache * chiCache = &(caches->constChiCache);
        /* The t-range has dependency of theta, however it is supposed
         * that following values give a good constant approximation
         * for \chi term */
        chiCache->tRange[0] = pow(ma, 4)
                            / (4*E0*E0);
        chiCache->tRange[1] = ma*ma;
        /* todo: consider this structure to be reentrant */
        caches->constChiCache.ff.ps = & caches->physPars;
        caches->constChiCache.ff.elastic_f   =   elastic_f;
        caches->constChiCache.ff.inelastic_f = inelastic_f;
        /* if \theta dependence for \chi is enabled, allocate reentrant
         * integration workspace object, otherwise set it to NULL. In latter
         * case, the dphmc_QAGS_integrate_iteratively() routine will use own
         * temporary workspace. */
        if( ! (0x1 & flags) ) {  /* chi-recalc disabled */
            caches->gslChiIntWS = NULL;
        } else {  /* chi-recalc enabled */
            caches->gslChiIntWS
                = gsl_integration_workspace_alloc( caches->chiIntPars.nnodes );
        }
        chiCache->chi = dphmc_QAGS_integrate_iteratively( &(caches->chiIntPars)
                                , dphmc_aprime_gsl_wrapper_chi_integrand
                                , &(caches->constChiCache.ff)
                                , chiCache->tRange[0], chiCache->tRange[1]
                                , &(chiCache->chiRelError), &(chiCache->chiAbsError)
                                , caches->gslChiIntWS );
    }
    /* Initialize GSL VEGAS integration workspace for full \sigma with NULL */
    caches->gslVegasSt = NULL;
    return caches;
    # if 0
    struct Workspace ** wsPtr = (struct Workspace **) wsPtr_
                   ,  * ws = NULL;
    struct dphmc_APrimeWSParameters * ps = NULL;
    /* Flags steering (re-)initialization of the procedure */
    const int doRecalculation             = 0x1  /* re-calc chi */
            , doGSLWorkspaceAllocation    = 0x2 | 0x1  /* re-allocate chi-int ws */
            , doParametersReassignment    = 0x4 | 0x1  /* re-assign user params */
            ;
    int routineFlags = 0x0;

    if( *wsPtr ) {  /* workspace is already allocated */
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
        if( ps->flags & 0x2 ) {
            ws->thetaMax = M_PI;
        } else {
            double val1 = sqrt( ps->massA_GeV*_DPhMC_CONST_electronMass_GeV )/ps->EBeam_GeV
                 , val2 = pow( ps->massA_GeV/ps->EBeam_GeV, 1.5 );
            ws->thetaMax = ps->maxThetaFactor * (val1 > val2 ? val1 : val2);
            if( ws->thetaMax > M_PI ) {
                ws->thetaMax = M_PI;
            }
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
    # endif
}

void
dphmc_aprime_delete( struct dphmc_APrimeWWCaches * caches ) {
    if( !caches ) return;
    if( caches->gslChiIntWS ) {
        /* Free GSL integration workspace for \chi */
        gsl_integration_workspace_free( caches->gslChiIntWS );
    }
    if( caches->gslVegasSt ) {
        /* Free GSL integration workspace for full \sigma */
        gsl_monte_vegas_free( caches->gslVegasSt );
    }
    free(caches);
}

/** For given nuclei parameters, calculates integrand function as a sum
 * of elastic and, optionally, inelastic components (when `inelastic` callback
 * is not NULL). See \cite Bjorken (A18,19) and \cite KimTsaiWW (3.20):
 * \f[
 *  \chi' = \frac{t - t_{min}}{t^2} [G_2^{el}(t) + [G_2^{inel}(t) ],
 * \f]
 * where elastic and inelastic components of electric form factor \f$G_2(t)\f$
 * are determined by the callbacks within the `ff` structure.
 *
 * @param t the virtuality variable
 * @param ff_ the instance of `dphmc_APrimeFormFactor`
 * @returns the \f$\chi'(t)\f$ value
 * @note: tMin = pow(ps->massA_GeV, 4)/(4*(ps->EBeam_GeV)*(ps->EBeam_GeV)) */
double
dphmc_aprime_gsl_wrapper_chi_integrand( double t, void * ff_ ) {
    struct dphmc_APrimeFormFactor * ff = (struct dphmc_APrimeFormFactor *) ff_;
    assert(ff->elastic_f);
    double res = ff->elastic_f( t, ff->ps );
    if( ff->inelastic_f ) {
        res += ff->inelastic_f( t, ff->ps );
    }
    return res;
}

/** Computes result of the following expression:
 * \f[
 * G_{2,el})(t) = (\frac{ a^2 t }{ 1 + a^2 t })^2 (\frac{1}{1 + t/d})^2 Z^2,
 * \f]
 * where $a = 111 Z^{-1/3}$, $d = 0.164 GeV^2 A^{-2/3}.$ */
double
dphmc_aprime_form_factor_elastic( double t
                                , struct dphmc_APrimePhysParameters * ps ) {
    assert( ps );
    const double  A = ps->A
               ,  Z = ps->Z
               ,  d = 0.164/pow(A, 2./3)
               ,  a = 111/(pow(Z, 1./3)*me)
               , a2 = a*a
               ;
    return pow( a2*t*Z/( (1 + a2*t)*(1 + t/d) ), 2);
}

/**Computes result of following expression:
 * \f[
 * G_{2,inel})(t) = (\frac{a^2 t}{1 + a^2 t})^2 (\frac{1 + (\mu_p^2 - 1)/(4 m_p^2)}
 *      { (1 + t / (0.71 GeV^2))^4 })^2 Z,
 * \f]
 * where \f$ a = 773 Z^{-2/3} / m_e \f$, \f$m_p\f$ is the proton mass
 * and \f$\mu_{p}\f$ is the proton magnetic moment.
 *
 * For most of the cases the contribution of inelastic part is expected to be
 * negligible.
 * */
double
dphmc_aprime_form_factor_inelastic( double t
                                  , struct dphmc_APrimePhysParameters * ps ) {
    assert( ps );
    const double   Z = ps->Z
               ,  ap = 773/(pow(Z, 2./3)*me)
               , ap2 = ap*ap
               ;
    return Z*pow( ( (ap2*t)/(1+ap2*t) )*( 1+t*(muP*muP-1)/(4*mp*mp) ) /
                                    pow(1 + t/.71, 4)
                  , 2 );
}

/**If caches were initialized for taking into account the \f$\theta\f$
 * dependence, calculates the corresponding \f$\chi\f$ value, otherwise
 * returns pre-computed constant approximation.
 *
 * @todo the \f$t_{max}\f$ is taken to be \f$m_{A'}\f$ that leads to wrong
 * order of equivalence b/w vituality values. Has to be changed after more
 * careful kinematics analysis.
 * */
double
dphmc_aprime_ww_photon_flux_for( double x
                               , double theta
                               , void * caches_ ) {
    assert( caches_ );
    struct dphmc_APrimeWWCaches * caches
        = (struct dphmc_APrimeWWCaches *) caches_;
    if( ! (caches->gslChiIntWS) ) {
        /* no ws allocd -- return const approx */
        return caches->constChiCache.chi;
    }
    /* compute chi for given \theta */
    const double E0 = caches->physPars.EBeam_GeV
               , ma = caches->physPars.massA_GeV
               , U = E0*E0 * theta*theta * x + ma*ma * (1 - x)/x + me*me*x
               , tLowSq = U/(2*E0*(1-x))
               ;
    return dphmc_QAGS_integrate_iteratively( &(caches->chiIntPars)
                                , dphmc_aprime_gsl_wrapper_chi_integrand
                                , &(caches->constChiCache.ff)
                                , tLowSq*tLowSq, caches->constChiCache.tRange[1]
                                , &(caches->constChiCache.chiRelError), &(caches->constChiCache.chiAbsError)
                                , caches->gslChiIntWS );
}

/**Calculates sigma according to (A12) formula \cite{JDBjorken}, returning
 * result of numerical calculation of the following value:
 * \f[
 * \frac{ d \sigma_{ 2 \rightarrow 3 } }{ d x d \cos{\theta_{A'}} } = \sin(x)
 * x (\frac{1-x+x^2/2}{U^2} + \frac{(1-x)^2 m_{A'}^2}{U^4} ( m_{A'}^2 - \frac{U x}{1-x} )),
 * \f]
 * where
 * \f[
 * U = E_0^2 \theta_{A'}^2 x + m_{A'}^2 \frac{1-x}{x} + m_e^2 x.
 * \f]
 *
 * Note the absence of \f$ \frac{1}{E_{0}^2 x} \f$ contrary to original
 * expression (moved to left).
 *
 * This relation includes calculation of constant term that is not needed for
 * simulation purposes. For the event generation purposes, consider the
 * "unnormed PDF" function dphmc_aprime_unnormed_pdf_a12() instead.
 *
 * @todo Add constant term in LaTeX formula here.
 * */
double
dphmc_aprime_cross_section_a12( double x
                              , double theta
                              , void * caches_ ) {
    assert(caches_);
    struct dphmc_APrimeWWCaches * caches
        = (struct dphmc_APrimeWWCaches *) caches_;
    const double chi = dphmc_aprime_ww_photon_flux_for( x, theta, caches_ );
    const struct dphmc_APrimePhysParameters * ps = &(caches->physPars);
    const double AM = ps->massA_GeV
               , E0 = ps->EBeam_GeV
               , am2 = AM*AM
               , U = pow(E0*theta, 2)*x + am2*(1 - x)/x + me*me*x
               , U2 = U*U
               , factor1 = (1 - x + x*x/2)/U2
               , factor2 = pow( (1 - x)*AM/U2, 2)
               , factor3 = am2 - (U*x)/(1 - x)
               , factor4 = 2 * caches->cstFact * chi * (E0 * E0) * x * sin(theta)
               , res = factor4 * (factor1 + factor2*factor3)
               ;
    assert(isfinite(res));
    return res;
}

/**Represents function \f$d \d sigma / d x\f$
 * corrsponding to dfferential by emitted A' energy cross section, that may be
 * useful in generators that does not require angular distribution.
 * @note this relation is derived with number of neglections, e.g. does not
 * take into account form factor dependency, imprecise \f$~O(m_{e})\f$, etc. */
double
dphmc_aprime_cross_section_a14( double x
                              , void * caches_ ) {
    assert(caches_);
    const struct dphmc_APrimeWWCaches * caches
        = (const struct dphmc_APrimeWWCaches *) caches_;
    const double chi = caches->constChiCache.chi;
    const struct dphmc_APrimePhysParameters * ps = &(caches->physPars);
    const double AM = ps->massA_GeV
               , factor4 = caches->cstFact * chi
               , res = factor4 * ( 1 - x + x*x/3. ) / ( AM*AM*(1 - x)/x + me*me*x )
               ;
    assert(isfinite(res));
    return res;
}


/** Getter for pre-computed \f$x_{low}\f$ (lower integration limit for relative
 * energy of emitted A'). */
double
dphmc_aprime_lower_cut_x( const void * caches_ ) {
    struct dphmc_APrimeWWCaches * caches
        = (struct dphmc_APrimeWWCaches *) caches_;
    return caches->xRange[0];
}

/** Getter for pre-computed \f$x_{up}\f$ (upper integration limit for relative
 * energy of emitted A'). */
double
dphmc_aprime_upper_cut_x( const void * caches_ ) {
    struct dphmc_APrimeWWCaches * caches
        = (struct dphmc_APrimeWWCaches *) caches_;
    return caches->xRange[1];
}

/** Getter for pre-computed \f$\theta_{max}\f$ (upper integration limit for
 * polar of emitted A'). */
double
dphmc_aprime_upper_cut_theta( const void * caches_ ) {
    struct dphmc_APrimeWWCaches * caches
        = (struct dphmc_APrimeWWCaches *) caches_;
    return caches->thetaMax;
}



static double
dphmc_aprime_ww_approx_sigma_a15( const struct dphmc_APrimeWWCaches * caches_ ) {
    const struct dphmc_APrimeWWCaches * caches
                = (struct dphmc_APrimeWWCaches *)(caches_);
    double c1 = me/caches->physPars.massA_GeV
         , c2 = caches->physPars.massA_GeV / caches->physPars.EBeam_GeV  /*ws->xRange[0]*/
         , upLim = c1 > c2 ? c1 : c2;
    upLim *= upLim;
    /*upLim = c2;  // NOTE: NA64 tooks it like that (no ^2, no max, just m_A/E_0)*/
    /* See (A15) and footnote 2 on page 2 of arXiv 1209.6083 */
    /* Note: I tried to perform similar integration procedure and got factor 6
     * instead of 3 in denominator after dx integration, so predicted output
     * should be ~4 times lesser than (A15) predicts. */
    return /*(8./3)*/ /*(2./3)*/ (1./3)
          * caches->cstFact
          * caches->constChiCache.chi
          * log( 1 / upLim )
          ;
}

double
dphmc_aprime_ww_fast_integral_estimation( void * caches_ ) {
    struct dphmc_APrimeWWCaches * caches
                = (struct dphmc_APrimeWWCaches *)(caches_);
    if( ! isnan(caches->integralEstVal) ) {
        return caches->integralEstVal;
    }
    /*assert( isnan(caches->integralEstVal, "lazy") );*/
    return caches->integralEstVal = dphmc_aprime_ww_approx_sigma_a15( caches_ );
}

static double
_static_ww_aprime_cs_wrapper( double * x, size_t dim, void * params ) {
    assert( 2 == dim );
    return dphmc_aprime_cross_section_a12( x[0], x[1], params );
}

/**Applies Gauss-Kronrod integration procedure to
 * `dphmc_aprime_cross_section_a14()`. Faster then
 * `dphmc_aprime_ww_full_numeric_2()`, but yields a less-precise value.
 * @note Same assumptions as for `dphmc_aprime_cross_section_a14()` take place.
 * @todo is still inefficient for MC generator. Consider to move it to
 * testing library */
double
dphmc_aprime_ww_full_numeric_1( const struct dphmc_IterativeQAGSParameters * qagsp
                              , void * caches_
                              , double * relErrPtr, double * absErrPtr ) {
    assert(caches_);
    struct dphmc_APrimeWWCaches * caches
        = (struct dphmc_APrimeWWCaches *) caches_;
    return dphmc_QAGS_integrate_iteratively( qagsp
                                           , dphmc_aprime_cross_section_a14, caches_
                                           , caches->xRange[0], caches->xRange[1]
                                           , relErrPtr, absErrPtr
                                           , NULL );
}

/**Performs VEGAS two-folded integration of
 * dphmc_aprime_unnormed_pdf_a12().
 * @todo move out random generator
 * @todo may be precise, but is extemely inefficient. Move it to testing?
 * @todo may be reasonable to wrap some adaptive loop around, until some
 * required chisq threshold is reached.
 */
double
dphmc_aprime_ww_full_numeric_2( void * caches_
                              , size_t nCalls
                              , double * absError
                              , double * chi2 ) {
    assert(caches_);
    struct dphmc_APrimeWWCaches * caches
        = (struct dphmc_APrimeWWCaches *) caches_;

    if( ! caches->gslVegasSt ) {
        caches->gslVegasSt = gsl_monte_vegas_alloc(2);
        int rc = gsl_monte_vegas_init( caches->gslVegasSt );
        assert( GSL_SUCCESS == rc );
    }
    # if 1
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc( T );
    # endif
    gsl_monte_function F = { _static_ww_aprime_cs_wrapper, 2, caches_ };
    double ranges[] = { caches->xRange[0], 0
                      , caches->xRange[1], caches->thetaMax
                      };
    double result, absError_;
    gsl_monte_vegas_integrate( &F
                             , ranges
                             , ranges + 2
                             , 2
                             , nCalls
                             , r
                             , caches->gslVegasSt
                             , &result
                             , &absError_ );
    if( absError ) *absError = absError_;
    if( chi2 ) *chi2 = gsl_monte_vegas_chisq(caches->gslVegasSt);
    return result;
}

struct dphmc_WSWithX {
    void * caches;
    double x;
};

/* Internal GSL wrapper for dphmc_aprime_cross_section_a12 */
static double
_dphmc_cs_theta_for_deriv( double theta, void * wsWithX_ ) {
    struct dphmc_WSWithX * wsWithX = (struct dphmc_WSWithX *) wsWithX_;
    return dphmc_aprime_cross_section_a12( wsWithX->x, theta, wsWithX->caches );
}

double
dphmc_aprime_ww_n_theta_function( double x, double theta
                                , void * caches
                                , double initStep
                                , double * abserrPtr ) {
    struct dphmc_WSWithX wsWithX = { caches, x };
    gsl_function F;
    F.function = _dphmc_cs_theta_for_deriv;
    F.params = &wsWithX;
    double result;
    int rc = gsl_deriv_central( &F, theta, initStep, &result, abserrPtr);
    assert(!rc);
    return result;
}



struct dphmc_ThetaMaxCache {
    double ma, ma2
         , E0
         , x, x2, xm1
         , a, b, d
         , e, chi
         ;
};

void *
dphmc_aprime_ww_new_sa_theta_max_cache( double x, const void * caches_ ) {
    assert(caches_);
    const struct dphmc_APrimePhysParameters * ps
            = &(((const struct dphmc_APrimeWWCaches *) caches_)->physPars);
    struct dphmc_ThetaMaxCache * r = malloc( sizeof(struct dphmc_ThetaMaxCache) );

    r->ma = ps->massA_GeV;
    r->ma2 = r->ma * r->ma;
    r->E0 = ps->EBeam_GeV;
    r->x = x;
    r->x2 = x*x;
    r->xm1 = 1 - x;
    r->a = r->xm1 + r->x2/2;
    r->b = r->xm1 * r->xm1 * r->ma2 * r->ma2;
    r->d = r->xm1 * r->ma2 * r->x;
    r->e = ps->epsilon;
    r->chi = ((const struct dphmc_APrimeWWCaches *) caches_)->constChiCache.chi;

    return r;
}

void
dphmc_aprime_ww_free_sa_theta_max_cache( void * cache_ ) {
    if( !cache_ ) return;
    free( (struct dphmc_ThetaMaxCache *) cache_ );
}

/**
 * \f[
 * \frac{1}{E_{0}^2 x} \frac{ d \sigma_{ 2 \rightarrow 3 } }{ d x d \cos{\theta_{A'}} }
 * \f]
 * Function is expected to be zero in proximity of maximum cross-section value. */
double
dphmc_aprime_ww_d2_theta_function( double theta
                                 , void * cache_ ) {
    assert(cache_);
    const struct dphmc_ThetaMaxCache * c
        = (const struct dphmc_ThetaMaxCache *) cache_;
    const double U = pow(c->E0 * theta, 2) * c->x + c->ma2 * c->xm1 / c->x + me*me* c->x
               , U2 = U*U
               , E02 = (c->E0 * c->E0)
               //, numerator = c->a * U2 - 2 * c->d * U + 3 * c->b
               //, denominator = c->a * U2 - c->d * U + c->b
               , cotTheta = 1./tan(theta)
               //, resFactor = cotTheta + 2 * E02 * c->x * (theta / U) * numerator / denominator
               /* following might be redundant for maximum look-up as
                * corresponds to theta-independent term */
               , betaAPrime = sqrt( 1 - pow(c->ma/c->E0, 2) )
               , f = c->b/U2 - c->d / U + c->a
               , C = 8 * alpha * pow( alpha*c->e, 2) * betaAPrime * (c->chi) * E02 * c->x
               , resFactor = cotTheta + 2 * E02 * c->x * (theta / U)*( (c->d - 2*c->b/U)/(U*f) - 2 )
               ;
    //return cotTheta - 2 * (c->E0 * c->E0) * c->x * (theta / U) * numerator / denominator;
    return C*sin(theta)*resFactor*f/U2 * _DPhMC_CONST_pbarn_per_GeV;
}

# define DPHMC_APRIME_WW_CS_D2THETA_ITERATION_MAX 100
# define DPHMC_APRIME_WW_CS_D2THETA_BRACKET_LIMIT 1e-6

/** Uses Brent-Dekker (GSL implem) root fnding algorithm over the analytical
 * derivative formula to locate the maximum emission cross section
 * \f$theta_{max}\f$ for given \f$x\f$ (iteratively).
 *
 * @param x \f$x\f$ for which the look-up procedure is performed
 * @param ws_ general A' WW approx cache object allocated
 *              with dphmc_init_aprime_cs_workspace()
 * @param cache_ (can be NULL) -- reentrant cache object
 *              allocated with dphmc_aprime_ww_free_sa_theta_max_cache()
 * @returns Emission angle \f$\theta\f$ if succeed, `nan("1")` if diverged,
 * `nan("2")` if did not converged during
 * `DPHMC_APRIME_WW_CS_D2THETA_ITERATION_MAX` iterations.
 * */
double dphmc_aprime_ww_theta_max_semianalytic( double x
                                             , const void * caches_
                                             , void * thetaMaxCache_ ) {
    assert(caches_);
    const struct dphmc_APrimeWWCaches * caches
        = ((const struct dphmc_APrimeWWCaches *) caches_);
    //const struct dphmc_APrimeWSParameters * ps = &(ws->parameterSet);  // XXX
    /* Set up the GSL wrapper around 2-derivative function */
    gsl_function F;
    F.function = dphmc_aprime_ww_d2_theta_function;
    if( ! thetaMaxCache_ ) {
        F.params = dphmc_aprime_ww_new_sa_theta_max_cache( x, caches_ );
    } else {
        F.params = thetaMaxCache_;
    }
    /* Set up the solver */
    double thetaLowEst = 1e-12
         , thetaUpEst = caches->thetaMax
         ;

    # if 0
    printf( "init: %e -> %e ; %e -> %e\n"
          , thetaLowEst
          , F.function( thetaLowEst, F.params )
          , thetaUpEst
          , F.function( thetaUpEst,  F.params ) );  // XXX
    # endif
    # ifndef NDEBUG
    /* Check that points straddle. We expect that if all the mathematical
     * routines were implemented correctly, the x/theta limits will always
     * contain maximum emission angle */
    if( F.function( thetaLowEst, F.params )
       *F.function( thetaUpEst,  F.params ) > 0. ) {
        DPhMC_msge( "Semi-analytic theta_max look-up failed:"
                    " border points do not straddle at x=%e.", x );
        return nan("3");
    }
    # endif

    gsl_root_fsolver * s = gsl_root_fsolver_alloc( gsl_root_fsolver_brent );
    gsl_root_fsolver_set( s, &F, thetaLowEst, thetaUpEst );
    /* Perform the iterative look-up for function zero */
    int status;
    size_t nIter = 0;
    double assumedMax;
    for( status = GSL_CONTINUE
       ; GSL_CONTINUE == status
       ; ++nIter ) {
        status = gsl_root_fsolver_iterate( s );
        # ifndef NDEBUG
        /* We expect that if all the mathematical routines were implemented
         * correctly, the derivative will always be a numeric value, so this
         * check is redundant for non-debug builds */
        if( GSL_SUCCESS != status ) {
            /* GSL_EBADFUNC means function returned nan or inf */
            DPhMC_msge( "Semi-analytic theta_max look-up failed:"
                        " Brent-Dekker iteration returned %d GSL exit code."
                      , status );
            assumedMax = nan("4");
            break;
        }
        # endif
        # if 0
        printf( "%3zu : %e -> %e ; %e -> %e\n", nIter
              , thetaLowEst
              , F.function( thetaLowEst, F.params )
              , thetaUpEst
              , F.function( thetaUpEst,  F.params ) );  // XXX
        # endif
        assumedMax = gsl_root_fsolver_root( s );
        thetaLowEst = gsl_root_fsolver_x_lower( s );
        thetaUpEst = gsl_root_fsolver_x_upper( s );

        status = gsl_root_test_interval( thetaLowEst
                                       , thetaUpEst
                                       , 0
                                       , DPHMC_APRIME_WW_CS_D2THETA_BRACKET_LIMIT );
        if( GSL_CONTINUE != status && GSL_SUCCESS != status ) {
            /* diverged */
            DPhMC_msge( "Semi-analytic theta_max look-up failed:"
                        " Brent-Dekker algo diverged on interval [%e, %e]."
                      , thetaLowEst, thetaUpEst );
            assumedMax = nan("1");
            break;
        }
        if( nIter > DPHMC_APRIME_WW_CS_D2THETA_ITERATION_MAX ) {
            /* out of iterations exit */
            DPhMC_msge( "Semi-analytic theta_max look-up failed:"
                        " Brent-Dekker algo failed to converge during"
                        " %d iterations."
                      , DPHMC_APRIME_WW_CS_D2THETA_ITERATION_MAX );
            assumedMax = nan("2");
            break;
        }
    }
    if( ! thetaMaxCache_ ) {
        dphmc_aprime_ww_free_sa_theta_max_cache( F.params );
    }
    gsl_root_fsolver_free( s );
    return assumedMax;
}

double
dphmc_aprime_ww_integrated_a14( double x, const void * caches_ ) {
    const struct dphmc_APrimeWWCaches * caches
        = ((const struct dphmc_APrimeWWCaches *) caches_);
    const struct dphmc_APrimePhysParameters * ps
        = &(caches->physPars);
    const double ma = ps->massA_GeV
               , ma2 = ma*ma
               , me2 = me*me
               , ma2mme2 = ma2 - me2
               , x2 = x*x
               , denom = sqrt( ma2 - 4*me2 )
               , fact = x2*me*me + 2*x*me2*ma2mme2
                      + log( (x - 1)*ma2 -x2*me2 )*ma2mme2*ma2mme2
                      + 2 * atanh( (ma2 - 2*x*me2)/( ma*denom ) )
                          * ma * (ma2*ma2 - 4*ma2*me2 + 3*me2*me2)
                          / denom
               ;
    return caches->cstFact*fact/(2*me2*me2*me2);  // TODO check that cstFact really contains everything (no)
}

# if 0
double
dphmc_aprime_ww_fast_integral_estimation( void * ws ) {
}

/**Used mainly by accept-reject (direct von Neumann) sampling method */
double
dphmc_aprime_ww_upper_limit( void * caches ) {
    //# warning "Here be dragons"  // TODO
    return 0.;
}
# endif

