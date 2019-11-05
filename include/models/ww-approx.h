# ifndef H_DPHMC_APRIME_WW_APPROX_CROSS_SECTION_H
# define H_DPHMC_APRIME_WW_APPROX_CROSS_SECTION_H

/**@file ww-approx.h
 * @brief Pure C computational routines based on code originally written
 * by D. Kirpichnikov.
 *
 * Header contains declarations of main computational routines combined
 * in order to provide probability density function (PDF) for production
 * reaction of hypothetical «dark photons» according to \cite Bjorken.
 *
 * Referencies:
 *  - JDBjorken // http://arxiv.org/pdf/0906.0580v1.pdf
 */

# include "dphmc-config.h"
# include "dphmc-inspection.h"

# include <stdlib.h>
# include <stdint.h>

/*
 * Constants
 * Since C compiler performs some run-time arithmetics, it is better to keep
 * these constants defined as preprocessor macros.
 */
/**@brief Mass of electron in GeVs
 * @details Shall be close to `CLHEP::electron_mass_c2 / CLHEP::GeV` */
# define _DPhMC_CONST_electronMass_GeV        5.10998950e-4
/**@brief Fine structure constant,
 * @details Shall be close to `CLHEP::fine_structure_const` */
# define _DPhMC_CONST_fineStructureConstant   7.2973525693e-3
/** Proton magnetic mu */
# define _DPhMC_CONST_muProton                2.79284734462
/** Proton mass in GeV
 * @details Shall be close to `CLHEP::proton_mass_c2 / CLHEP::GeV` */
# define _DPhMC_CONST_protonMass_GeV          9.38272013e-1
/**@brief picobarn in GeV translation unit
 * @details Given by \f$ GeV^{-2} = 0.3894 mb = 0.3894 * 1e9 pb \f$
 * */
# define _DPhMC_CONST_pbarn_per_GeV           3.894e+8

# ifdef __cplusplus
extern "C" {
# endif

/**@brief A parameters set structure for integrating with GSL Gauss-Kronrod algo.
 * @details This integration algorith is common for integrating 1D functions
 * in library. Iterative procedure over the GSL function is applied to gain
 * numerical estimation, gradually loosing the requirement on the relative
 * error. */
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

/**@brief Common numerical Cross-Section calculus parameterization.
 * @details Represents a set of common parameters for calculating cross-section
 * (full and differential) values of Weisacker-Williams approximation of
 * A' production.
 * */
struct dphmc_APrimePhysParameters {
    uint16_t Z;  /**< Target Z*/
    double A  /**< Target A*/
         , massA_GeV  /**< A' (hypothetical) mass, GeV */
         , EBeam_GeV  /**< Incident beam energy, GeV */
         , epsilon  /**< Mixing constant for A' particle */
         , factor  /**< (dev) Multiplication factor for cross-section */
         , maxThetaFactor  /**< \f$\theta_{max}\f$ proportionality coefficient */
         ;
};

/** Prints the human-readable form of the workspace parameters to given
 * stream. */
void
dphmc_print_aprime_ws_parameters( FILE *
                                , const struct dphmc_APrimePhysParameters * );

/** Checks the input physics parameters for validity to WW approximation
 * approach. */
int
dphmc_aprime_ww_check_phys_parameters_are_valid( const struct dphmc_APrimePhysParameters * );

/**@brief Common internal caches of all the procedures */
struct dphmc_APrimeWWCaches;

/** Initializes supplementary parameter set for cross-section calculation. */
struct dphmc_APrimeWWCaches *
dphmc_aprime_new( const struct dphmc_APrimePhysParameters * newPs
                , const struct dphmc_IterativeQAGSParameters * chiIntParameters
                , double (*  elastic_f)( double, struct dphmc_APrimePhysParameters * )
                , double (*inelastic_f)( double, struct dphmc_APrimePhysParameters * )
                , int flags );

/** Frees supplementary parameters set for chi integration. */
void dphmc_aprime_delete( struct dphmc_APrimeWWCaches * );



/** Represents parameters set for form factor calculus within WW approx. */
struct dphmc_APrimeFormFactor {
    /** Physics parameters used for FF calculus */
    struct dphmc_APrimePhysParameters * ps;
    /** Elastic part of the form factor, integrated over virtuality \f$t\f$. */
    double (*  elastic_f)( double, struct dphmc_APrimePhysParameters * ps );
    /** Inelastic part of the form factor, integrated over virtuality \f$t\f$ */
    double (*inelastic_f)( double, struct dphmc_APrimePhysParameters * ps );
};

/**Calculates elastic part of electric component of the form factor given
 * by (A18) \cite Bjorken. */
double
dphmc_aprime_form_factor_elastic( double t
                                , struct dphmc_APrimePhysParameters * ps );

/**Inelastic part of electric component of the form factor given
 * by (A19) \cite Bjorken. */
double
dphmc_aprime_form_factor_inelastic( double t
                                  , struct dphmc_APrimePhysParameters * ps );

/**GSL-compatible wrapper for \f$\chi\f$ integrand function. */
double
dphmc_aprime_gsl_wrapper_chi_integrand( double t, void * ff_ );

/** ... */
double
dphmc_aprime_ww_photon_flux_for( double x
                               , double theta
                               , void * caches_ );


/** Calculates cross-section according to \cite{JDBjorken} (A12). */
double dphmc_aprime_cross_section_a12( double x
                                     , double theta
                                     , void * ws );

/** Full implementation of (A14) of \cite Bjorken. */
double dphmc_aprime_cross_section_a14( double x
                                     , void * ws );



/** Returns \f$x_{min}\f$ cut (from common sense: \f$m_{A'}/E_0\f$) */
double dphmc_aprime_upper_cut_x( const void * caches );

/** Returns \f$x_{max}\f$ cut (according to formula on p. 4 \cite{Bjorken}) */
double dphmc_aprime_lower_cut_x( const void * caches );

/** Returns \f$\theta_{A', max}\f$ cut (according to formula (9) \cite{Bjorken}) $ */
double dphmc_aprime_upper_cut_theta( const void * caches );



/** Returns integral estimation obtained according to (A16) formula
 * \cite{Bjorken}. */
double dphmc_aprime_ww_fast_integral_estimation( void * ws );

/**Numerical full cross-section estimate w.r.t. WW A' approximation (by
 * x only). */
double
dphmc_aprime_ww_full_numeric_1( const struct dphmc_IterativeQAGSParameters * qagsp
                              , void * caches_
                              , double * relErrPtr, double * absErrPtr );

/**Numerical full cross-section estimate w.r.t. WW A' approximation. */
double dphmc_aprime_ww_full_numeric_2( void * ws, size_t nCalls
                                     , double * absError, double * chi2 );

/**Returns estimation of an upper bound of the differential
 * cross-section given by dphmc_aprime_cross_section_a12(). */
double dphmc_aprime_ww_upper_limit( void * ws );


/**Numerical estimation of differential cross-section differentiated
 * by \f$\theta\f$. */
double dphmc_aprime_ww_n_theta_function( double x, double theta
                                       , void * ws
                                       , double initStep
                                       , double * abserrPtr );


// TODO: in a dedicated group (semianalytic estimation of maximal theta)

/** Allocates cache for semi-analytic estimation of \f$\theta_{max}\f$. */
void * dphmc_aprime_ww_new_sa_theta_max_cache( double x, const void * caches );

/** Second order derivative of cross-section expresion */
double dphmc_aprime_ww_d2_theta_function( double theta
                                        , void * cache_ );

/** Frees cache for semi-analytic estimation of \f$\theta_{max}\f$. */
void dphmc_aprime_ww_free_sa_theta_max_cache( void * cache );

/** Semi-analytical estimation of the maximum emission angle for
 * given \f$x\f$. */
double dphmc_aprime_ww_theta_max_semianalytic( double x
                                             , const void * ws
                                             , void * cache_ );

/** ... */
double dphmc_aprime_ww_integrated_a14( double x, const void * caches_ );

# ifdef __cplusplus
}
# endif

# endif  /* H_DPHMC_APRIME_WW_APPROX_CROSS_SECTION_H */

