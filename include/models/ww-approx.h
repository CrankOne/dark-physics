# ifndef H_APRIME_CROSS_SECTIONS_H
# define H_APRIME_CROSS_SECTIONS_H

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
# define _DPhMC_CST_electronMass_GeV        5.10998950e-4
/**@brief Fine structure constant,
 * @details Shall be close to `CLHEP::fine_structure_const` */
# define _DPhMC_CST_fineStructureConstant   7.2973525693e-3
/** Proton magnetic mu */
# define _DPhMC_CST_muProton                2.79284734462
/** Proton mass in GeV
 * @details Shall be close to `CLHEP::proton_mass_c2 / CLHEP::GeV` */
# define _DPhMC_CST_protonMass_GeV          9.38272013e-1
/**@brief picobarn in GeV translation unit
 * @details Given by \f$ GeV^{-2} = 0.3894 mb = 0.3894 * 1e9 pb \f$
 * */
# define _DPhMC_CST_pbarn_per_GeV           3.894e+8

# ifdef __cplusplus
extern "C" {
# endif

/** Numerical Cross-Section calculus parameterization. */
struct dphmc_APrimeWSParameters {
    uint16_t Z;  /**< Target Z*/
    double A  /**< Target A*/
         , massA_GeV  /**< A' (hypothetical) mass, GeV */
         , EBeam_GeV  /**< Incident beam energy, GeV */
         , epsilon  /**< Mixing constant for A' particle */
         , factor  /**< Multiplication factor for cross-section */
         ;

    double epsabs  /**< epsabs Gauss-Kronrod absolute error (see GSL QAGS) */
         , epsrel  /**< epsabs Gauss-Kronrod relative error (see GSL QAGS) */
         , epsrelIncFt  /**< epsrel Gauss-Kronrod relative error inc factor */
         ;
    size_t limit  /**< maximum number of subintervals (see GSL QAGS) */
         , nnodes  /**< GSL integration workspace nodes num */
         ;
};

/** Prints the human-readable form of the workspace parameters to given
 * stream. */
void
dphmc_print_aprime_ws_parameters( FILE *
                                , const struct dphmc_APrimeWSParameters * );

/** Initializes supplementary parameter set for cross-section calculation. */
int dphmc_init_aprime_cs_workspace( const struct dphmc_APrimeWSParameters *
                                  , void ** wsPtr );

/** Frees supplementary parameters set for chi integration. */
void dphmc_free_aprime_cs_workspace( void * );

/**Integrated photons flux for A' CS calculations. */
double dphmc_aprime_chi( double * absError,
            double * relError,
            void * ws );

/** Calculates cross-section according to \cite{JDBjorken} (A12). */
double dphmc_aprime_cross_section_a12( double x
                                     , double theta
                                     , void * ws );

/** Volatile part of (A12) of \cite Bjorken. */
double dphmc_aprime_unnormed_pdf_a12( double x
                                    , double theta
                                    , void * ws );

/** Full implementation of (A12) of \cite Bjorken. */
double dphmc_aprime_cross_section_a14( double x
                                     , void * ws );

/** Returns \f$x_{min}\f$ cut (from common sense: \f$m_{A'}/E_0\f$) */
double dphmc_aprime_upper_cut_x( void * ws );

/** Returns \f$x_{max}\f$ cut (according to formula on p. 4 \cite{Bjorken}) */
double dphmc_aprime_lower_cut_x( void * ws );

/** Returns \f$\theta_{A', max}\f$ cut (according to formula (9) \cite{Bjorken}) $ */
double dphmc_aprime_upper_cut_theta( void * ws );

/** Returns integral estimation obtained according to (A16) formula \cite{Bjorken} */
double dphmc_aprime_ww_fast_integral_estimation( void * ws );

/**Numerical full cross-section estimate w.r.t. WW A' approximation (by
 * x only). */
double
dphmc_aprime_ww_full_numeric_1( void * ws_
                              , double epsabs
                              , double epsrel
                              , double epsrelIncFt
                              , size_t limit, size_t nnodes
                              , double * relError );

/**Numerical full cross-section estimate w.r.t. WW A' approximation. */
double dphmc_aprime_ww_full_numeric_2( void * ws, size_t nCalls
                                     , double * absError, double * chi2 );

/**Returns estimation of an upper bound of the differential
 * cross-section given by dphmc_aprime_cross_section_a12(). */
double dphmc_aprime_ww_upper_limit( void * ws );

# ifdef __cplusplus
}
# endif

# endif  /* H_APRIME_CROSS_SECTIONS_H */

