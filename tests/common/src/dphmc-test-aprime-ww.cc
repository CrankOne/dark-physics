# include "dphmc-test-aprime-ww.hh"

# include <stdexcept>
# include <cassert>
# include <cmath>
# include <tuple>

namespace dphmc {
namespace test {

double
chi_bjorken( double * x, double * p ) {
    const double EBeam = p[0],
                 APMass = x[0];
    // These parameters are somewhat comon ones. They are hardcoded here since
    // for \chi tests we usually want a reproducible tests.
    aprime::PhysParameters aps = {
            /* Z ............... */ (uint16_t) p[2], //p[2] 74,
            /* A ............... */ p[1],  //183.84,
            /* mass A', GeV .... */ APMass,
            /* E beam, GeV ..... */ EBeam,
            /* mixing factor ... */ 1e-4,
            /* re-norm factor .. */ 1,
            /* theta cut-off c . */ 100
    };
    IterativeQAGSParameters qagsp = {
            /* epsabs .......... */ 1e-12,
            /* epsrel .......... */ 1e-12,
            /* epsrelIncFt ..... */ 1.1,
            /* limit ........... */ (size_t) 1e3,
            /* nnodes .......... */ (size_t) 1e3
    };
    aprime::WWCaches * caches = dphmc_aprime_new( &aps, &qagsp
                                                , dphmc_aprime_form_factor_elastic
                                                , dphmc_aprime_form_factor_inelastic
                                                , 0x0 );
    if( dphmc_error || !caches ) {
        if(caches)
            dphmc_aprime_delete( caches );
        throw std::runtime_error( "Failed to initialize CS-workspace for A'." );
    }
    double absErr, relErr;
    double chiRes = dphmc_aprime_ww_photon_flux_for( std::nan("unused")
                                                   , std::nan("unused")
                                                   , caches );
    dphmc_aprime_delete( caches );
    return chiRes/(p[2]*p[2]);
}



std::vector<DCrossSectionPoint>
aprime_dcs_ww( const aprime::PhysParameters & apcsPs
             , const IterativeQAGSParameters & qagsp
             , size_t xNPts, size_t tNPts
             , double * xRange, double * thetaRange ) {
    assert( xNPts > 1 );
    assert( tNPts > 1 );
    aprime::WWCaches * caches = dphmc_aprime_new( &apcsPs, &qagsp
                                        , dphmc_aprime_form_factor_elastic
                                        , dphmc_aprime_form_factor_inelastic
                                        , 0x0 );

    double & xLow = xRange[0],      & xUp = xRange[1]
         , & tLow = thetaRange[0],  & tUp = thetaRange[1]
         ;

    xLow = dphmc_aprime_lower_cut_x(caches);
    xUp = dphmc_aprime_upper_cut_x(caches);
    tUp = dphmc_aprime_upper_cut_theta(caches);
    tLow = 1e-9;  // todo: in case of Born approximation, has to be mitigated

    assert( xLow < xUp );
    assert( tLow < tUp );
    std::vector<DCrossSectionPoint> pts;

    for( ExpValues evX( xNPts, 1 - xUp, 1 - xLow )
       ; evX.good()
       ; ++evX ) {
        double x = 1 - evX.value();
        void * saThetaCache = dphmc_aprime_ww_new_sa_theta_max_cache( x, caches );
        double absErr;
        for( ExpValues evTheta( tNPts, tLow, tUp )
           ; evTheta.good()
           ; ++evTheta ) {
            pts.push_back( { evX.value()
                           , evTheta.value()
                           , dphmc_aprime_cross_section_a12( x, evTheta.value(), caches )
                           , dphmc_aprime_ww_d2_theta_function( evTheta.value(), saThetaCache )
                           , dphmc_aprime_ww_n_theta_function( x, evTheta.value(), caches, 1e-9, &absErr )
                           , absErr
                           } );
            assert( pts.size() <= xNPts*tNPts );
        }
        // Free temporary caches valid for certain x only
        dphmc_aprime_ww_free_sa_theta_max_cache( saThetaCache );
    }
    dphmc_aprime_delete( caches );
    return pts;
}



std::vector< std::tuple<double, double, double> >
aprime_dcs_max( const aprime::PhysParameters & apcsPs
              , const IterativeQAGSParameters & qagsp
              , size_t xNPts ) {
    assert( xNPts > 1 );
    aprime::WWCaches * caches = dphmc_aprime_new( &apcsPs, &qagsp
                                        , dphmc_aprime_form_factor_elastic
                                        , dphmc_aprime_form_factor_inelastic
                                        , 0x0 );

    double xLow = dphmc_aprime_lower_cut_x(caches)
         , xUp = dphmc_aprime_upper_cut_x(caches);
    assert( xLow < xUp );

    std::vector< std::tuple<double, double, double> > pts;
    const double tX1    = pow( ( 1 - xLow)/(1 - xUp), 1./(xNPts-1) );
    for( double x_ = (1 - xUp); x_ <= (1 - xLow)*(1 + 1e-3); x_ *= tX1 ) {
        double x = 1 - x_;
        void * saThetaCache = dphmc_aprime_ww_new_sa_theta_max_cache( x, caches );
        double thetaMax = dphmc_aprime_ww_theta_max_semianalytic( x, caches, saThetaCache );
        dphmc_aprime_ww_free_sa_theta_max_cache( saThetaCache );
        if( std::isfinite( thetaMax ) ) {
            pts.push_back( std::tuple<double, double, double>( x_, thetaMax
                        , dphmc_aprime_cross_section_a12( x, thetaMax, caches ) ) );
        }
    }
    dphmc_aprime_delete( caches );
    return pts;
}

}  // namespace ::dphmc::test
}  // namespace ddphmc

