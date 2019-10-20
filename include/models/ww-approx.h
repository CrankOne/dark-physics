# ifndef H_APRIME_CROSS_SECTIONS_H
# define H_APRIME_CROSS_SECTIONS_H

/**\file aprime-cs.h
 * \brief Pure C computational routines based on code originally written
 *        by D. Kirpichnikov.
 *
 * Header contains declarations of main computational routines combined
 * in order to provide probability density function (PDF) for production
 * reaction of hypothetical «dark photons» according to \cite{JDBjorken}.
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
/** Mass of electron in GeVs
 * CLHEP::electron_mass_c2 / CLHEP::GeV */
# define _DPhMC_CST_electronMass_GeV        5.10998950e-4
/** Fine structure constant
 * CLHEP::fine_structure_const */
# define _DPhMC_CST_fineStructureConstant   7.2973525693e-3
/** Proton magnetic mu */
# define _DPhMC_CST_muProton                2.79284734462
/** Proton mass in GeV
 * CLHEP::proton_mass_c2 / CLHEP::GeV */
# define _DPhMC_CST_protonMass_GeV          9.38272013e-1
/** picobarn in GeV translation unit
 * GeV^(-2) = 0.3894 mb = 0.3894 * 1e9 pb
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

/** Initializes supplementary parameter set for cross-section calculation.
 *
 * The physical parameters should be provided by pointer on corresponding
 * const structure. This data is then copied to the "workspace" object that
 * stores also some intermediate caches needed to compute the cross-section
 * values.
 *
 * @param n nodes for GSL integration workspace (see GSL QAGS);
 * @returns -1 on error, 0 on success
 * */
int dphmc_init_aprime_cs_workspace( const struct dphmc_APrimeWSParameters *
                                  , void ** wsPtr );

/** Frees supplementary parameters set for chi integration. */
void dphmc_free_aprime_cs_workspace( void * );

/**@brief Integrated photons flux for A' CS calculations.
 *
 * For given nuclei parameters, calculates integrated photons flux
 * involving GSL integration procedure. On error, raises standard Goo
 * exception with third-party error code.
 *
 * see \ref{JDBjorken} (A17). */
double dphmc_aprime_chi( double * absError,
            double * relError,
            void * ws );

/**@brief Calculates cross-section according to \cite{JDBjorken} (A12).
 * @details Calculates sigma according to (A12) formula \cite{JDBjorken}, returning
 * result of numerical calculation of the following value:
 * \f[
 * \frac{1}{E_{0}^2 x} \frac{ d \sigma_{ 2 \rightarrow 3 } }{ d x d \cos{\theta_{A'}} }
 * \f]
 *
 * This relation includes calculation of constant term that is not needed for
 * simulation purposes. For the event generation purposes, se the "unnormed
 * PDF" function dphmc_aprime_unnormed_pdf_a12() instead.
 * */
double dphmc_aprime_cross_section_a12( double x
                                     , double theta
                                     , void * ws );

/**@brief TODO
 * @details TODO */
double dphmc_aprime_unnormed_pdf_a12( double x
                                    , double theta
                                    , void * ws );

/**@brief TODO
 * @details TODO */
double dphmc_aprime_cross_section_a14( double x
                                     , void * ws );

/** Returns $x_{min}$ cut (from common sense: $m_{A'}/E_0$) */
double dphmc_aprime_upper_cut_x( void * ws );

/** Returns $x_{max}$ cut (according to formula on p. 4 \cite{Bjorken}) */
double dphmc_aprime_lower_cut_x( void * ws );

/** Returns $\theta_{A', max} cut (according to formula (9) \cite{Bjorken}) $ */
double dphmc_aprime_upper_cut_theta( void * ws );

/** Returns integral estimation obtained according to (A16) formula \cite{Bjorken} */
double dphmc_aprime_ww_fast_integral_estimation( void * ws );

double
dphmc_aprime_ww_full_numeric_1( void * ws_
                              , double epsabs
                              , double epsrel
                              , double epsrelIncFt
                              , size_t limit, size_t nnodes
                              , double * relError );

/**@brief Numerical full cross-section estimate w.r.t. WW A' approximation.
 * @details Performs VEGAS integration internally.
 * @todo: may be precise, but inefficient. Move it to testing? */
double dphmc_aprime_ww_full_numeric_2( void * ws, size_t nCalls
                                     , double * absError, double * chi2 );

# ifdef __cplusplus
}
# endif

# endif  /* H_APRIME_CROSS_SECTIONS_H */

