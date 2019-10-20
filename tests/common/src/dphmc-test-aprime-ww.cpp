# include "dphmc-test-aprime-ww.hpp"

# include <stdexcept>
# include <cassert>
# include <cmath>

namespace dphmc {
namespace test {

double
chi_bjorken( double * x, double * p ) {
    const double EBeam = p[0],
                 APMass = x[0];
    void * ws = NULL;
    // These parameters are somewhat comon ones. They are hardcoded here since
    // for \chi tests we usually want a reproducible tests.
    dphmc_APrimeWSParameters aps = {
            /* Z ............... */ (uint16_t) p[2], //p[2] 74,
            /* A ............... */ p[1],  //183.84,
            /* mass A', GeV .... */ APMass,
            /* E beam, GeV ..... */ EBeam,
            /* mixing factor ... */ 1e-4,
            /* re-norm factor .. */ 1,
            /* epsabs .......... */ 1e-12,
            /* epsrel .......... */ 1e-12,
            /* epsrelIncFt ..... */ 1.1,
            /* limit ........... */ (size_t) 1e3,
            /* nnodes .......... */ (size_t) 1e3
    };
    dphmc_init_aprime_cs_workspace( &aps, &ws );
    if( dphmc_error ) {
        dphmc_free_aprime_cs_workspace( ws );
        throw std::runtime_error( "Failed to initialize CS-workspace for A'." );
    }
    double absErr, relErr;
    double chiRes = dphmc_aprime_chi( &absErr, &relErr, ws );
    dphmc_free_aprime_cs_workspace( ws );
    return chiRes/(p[2]*p[2]);
}

std::vector<DCrossSectionPoint>
aprime_dcs_ww( const dphmc_APrimeWSParameters & apcsPs
             , size_t xNPts, size_t tNPts
             , double * xRange, double * thetaRange ) {
    assert( xNPts > 1 );
    assert( tNPts > 1 );
    void * ws = NULL;
    dphmc_init_aprime_cs_workspace( &apcsPs, &ws );

    double & xLow = xRange[0],      & xUp = xRange[1]
         , & tLow = thetaRange[0],  & tUp = thetaRange[1]
         ;

    xLow = dphmc_aprime_lower_cut_x(ws);
    xUp = dphmc_aprime_upper_cut_x(ws);
    tUp = dphmc_aprime_upper_cut_theta(ws);
    tLow = 1e-9;  // todo: in case of Born approximation, has to be mitigated

    assert( xLow < xUp );
    assert( tLow < tUp );
    std::vector<DCrossSectionPoint> pts;
    # if 0  // evenly-spaced
    const double xStep = (xUp - xLow)/(xNPts - 1)
               , tStep = (tUp - tLow)/(tNPts - 1);
    assert( xStep > 0 );
    assert( tStep > 0 );
    for( double x = xLow; x <= xUp; x += xStep ) {
        for( double theta = tLow; theta <= tUp; theta += tStep ) {
            assert( pts.size() <= xNPts*tNPts );
            pts.push_back( {x, theta, dphmc_aprime_cross_section_a12( x, theta, ws )} );
        }
    }
    # else  // logarithmically-spaced
    const double tX     = pow( ( 1 - xLow)/(1 - xUp), 1./(xNPts-1) )
               , tTheta = pow( tUp/tLow, 1./(tNPts-1) )
               ;
    for( double x_ = (1 - xUp); x_ <= (1 - xLow)*(1 + 1e-3); x_ *= tX ) {
        for( double theta = tLow; theta <= tUp*(1 + 1e-3); theta *= tTheta ) {
            double x = 1 - x_;
            pts.push_back( {1 - x, theta, dphmc_aprime_cross_section_a12( x, theta, ws )} );
            assert( pts.size() <= xNPts*tNPts );
        }
    }
    # endif
    dphmc_free_aprime_cs_workspace( ws );
    return pts;
}

}  // namespace ::dphmc::test
}  // namespace ddphmc

