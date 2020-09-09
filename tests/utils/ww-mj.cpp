# include "dphmc-test-fixture.hh"
# include "dphmc-test-app.hh"

# include <iostream>
# include <cmath>
# include <cstring>

# define me     _DPhMC_CONST_electronMass_GeV

void
usage_info( std::ostream & os, const std::string & appName ) {
    os << "Usage:" << std::endl
       << "    $ " << appName << " [1|2] <outFile>" << std::endl
       << "Application samples the WW cross section values for certain physics"
       << " and writes the resulting 100-bin histogram into a file."
       << std::endl;
}


/** Double accept-reject scheme is applied here to generate the A' emission
 * parameters */
int
dphmc_aprime_ww_sample_emission( struct dphmc_URandomState * uRandom
                               , double E0
                               , double ma
                               , double * xPtr
                               , double * thetaPtr ) {
    int rc;
    # define URANDOM(v) DPHMC_RIF( uRandom->urandom_f(uRandom->statePtr, &v) );
    // pre-calculate common sub-expressions
    const double E02 = E0*E0  // initiating particle energy squared, E_0^2
               , ma2 = ma*ma  // squared A' mass, m_{A'}^2
               , E02ma2 = (E0/ma)*(E0/ma)  // squared ratio of E_0 and m_{A'}
               , ma4 = ma2*ma2  // m_{A'}^4
               , mema2 = (me/ma)*(me/ma)  // \frac{me^2}{ma^2}
               , pi2 = M_PI*M_PI
               ;
    {
        // Sample X according to 2nd majorant (M_{2,x}(x))
        double x
             , mxx2
             ;
        do {
            double mj1, mj2, conjX;
            // Generate a pair of random numbers
            URANDOM(x);
            URANDOM(conjX);
            // Map one uniform random variable from pair to x wrt M_{2,x}(x)
            // using reverse transform method:
            //  x = \frac{1}{2 E_0} \frac{1-(m_e/m_{A'})^{2 u}}{1-(m_e/m_{A'})^{2}}
            x = (1. - pow(mema2, x)) / (1. - mema2);
            mxx2 = (1-x)/(x*x);
            // Get the desired value M_{1,x}(x) ...
            mj1 = (pi2/(2*ma4*x))/( (mema2 + mxx2)
                                  * (mema2 + E02ma2 * pi2 + mxx2));
            // ... and current M_{2,x}(x) value.
            mj2 = (1/(2 * E02 * ma2 * x * x))/ (mxx2 + mema2) ;
            if( conjX*mj2 > mj1 ) {
                // continue until got one is above the desired
                continue;
            }
        } while(0);
        // Sample theta using inverse function of M_{2}(x, \theta) with x
        // taken as a parameter
        // TODO
        // Get the desired value V(x, \theta)
        // TODO
        // ... and current M_{1}(x, \theta)
        // TODO
    } // while
}

# if 0
class MajorantTestWrite : public {
    void fill_hst() {
        size_t nSamples = 1e5;
        size_t histogram[100];

        size_t nBins = sizeof(histogram)/sizeof(*histogram);
        double x;
        memset(histogram, 0, sizeof(histogram));
        for( int i = 0; i < nSamples; ++i ) {
            dphmc_aprime_ww_mj1_sample_x( &_rgs, _cachesPtr, &x );
            ASSERT_LE( x, 1 );
            ASSERT_GE( x, 0 );
            if( x < 1. ) {
                ++histogram[int(x*nBins)];
            } else {
                ++histogram[int((nBins-1))];
            }
        }
        // ...
    }
};
# endif

int
main(int argc, char * argv[]) {

    return EXIT_SUCCESS;
}
