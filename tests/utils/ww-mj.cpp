#include "dphmc-test-fixture.hh"
#include "dphmc-hst.hh"
#include "dphmc-integrw.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <cstring>
#include <unistd.h>
#include <cstdlib>
#include <fstream>
#include <getopt.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#define me _DPhMC_CONST_electronMass_GeV
# define URANDOM(v) DPHMC_RIF( uRandom->urandom_f(uRandom->statePtr, &v) );

class Sampler {
public:
    const double E0, ma;
//protected:
    double E02  ///< initiating particle energy squared, \f$E_0^2\f$
         , ma2  ///< squared A' mass, \f$m_{A'}^2\f$
         , E02ma2  ///< squared ratio of \f$E_0\f$ and \f$m_{A'}\f$
         , ma4  ///< \f$\m_{A'}^4\F$
         , mema2  ///< \f$\frac{me^2}{ma^2}\f$
         , pi2  ///< \f$\pi^2\f$
         , mj2Norm  ///< norm of \f$M_{2,x}(x)\f$
         ;
public:
    Sampler( double E0_, double ma_ ) : E0(E0_), ma(ma_) {
        E02 = E0*E0;
        ma2 = ma*ma;
        E02ma2 = (E0/ma)*(E0/ma);
        ma4 = ma2*ma2;
        mema2 = (me/ma)*(me/ma);
        pi2 = M_PI*M_PI;
        mj2Norm = log(1/mema2)/(2*E02*(1-mema2)*ma2);
    }

    /** This optimistic estimation for a random \f$x\f$ is derived by integrating
     * the majorant function:
     *
     * \f[
     * M_{x,2} (x) = (m_e^2 x + m_{A'}^2 (1-x) )^(-1),
     * \f]
     *
     * which yields the following reverse function for x:
     *
     * \f[
     * x = \frac{1}{2 E_0} \frac{1-(m_e/m_{A'})^{2 u}}{1-(m_e/m_{A'})^{2}}
     * \f]
     *
     * Used to derive X in two-step modified von Neumann method on A' with WW.
     * */
    double mj2_rev_x( double u ) const {
        return ( (1. - pow(mema2, u))
               / (1. - mema2)
               );
    }

    /** The \f$M_{2,x}(x) > \f$\int\limits_{0}{\pi} M_{1,x} d \theta\f$ for
     * $x \in [0:1]$ being thus a "majorant for integrated majorant". */
    double mj2( double x ) const {
        return (1/(2 * E02 * ma2 * x * x))
              / ( (1-x)/(x*x) + mema2)
              ;
    }

    /** Returns value of the "first" x-majorant (2D majorant integrated over theta):
     *
     * \f[
     * M_{x,1} = \frac{\pi^2 }{2 m_{A'}^4 x} /
     *      (( \frac{m_e^2}{m_{A'}}^2 + \frac{1-x}{x^2} )
     *       ( (\frac{m_e^2}{m_{A'}}^2 + \frac{E_0^2}{m_{A'}^2} \pi^2) + \frac{1-x}{x^2} ))
     * \f]
     *
     * */
    double mj1(double x) const {
        double mxx2 = (1-x)/(x*x)
             , factor1 = pi2/(2*ma2*ma2*x)
             , factor2 = (mema2 + mxx2)
             , factor3 = (mema2 + E02ma2 * pi2 + mxx2)
             ;
        return factor1/(factor2*factor3);
    }

    /** Uses \f$M_{x,2}(u)\f$ (defined by `dphmc_aprime_ww_mj2_rev_x()`)
     * to sample \f$\tilde{x}\f$ ("optimistic" \f$x\f$) according to \f$M_{x,1}(u)\f$
     * (defined by `dphmc_aprime_ww_mj1_x()`).
     *
     * Returned \f$\tilde{x}\f$ has to be used further as a parameter to first-order
     * majorant \f$M_{x,theta,1}(x, \theta)\f$ to sample A' kinematics.
     *
     * \warning Since \f$M_{x,1}(u)\f$ and \f$M_{x,2}(u)\f$ are considered to be
     * very close, no limit on iteration is done at this routine. */
    double sample_x_mj1( struct dphmc_URandomState * uRandom ) {
        int rc;
        // Sample X according to 2nd majorant (M_{2,x}(x))
        double x
             , mxx2
             ;
        double m1, m2, conjX;
        do {
            // Generate a pair of random numbers
            URANDOM(x);
            URANDOM(conjX);
            // Map one uniform random variable from pair to x wrt M_{2,x}(x)
            // using reverse transform method:
            //  x = \frac{1}{2 E_0} \frac{1-(m_e/m_{A'})^{2 u}}{1-(m_e/m_{A'})^{2}}
            x = mj2_rev_x(x);
            // Get the desired value M_{1,x}(x) ...
            m1 = mj1(x);
            // ... and current M_{2,x}(x) value.
            m2 = mj2(x);
            if( conjX*m2 <= m1 ) {
                return x;
            }
        } while(true);
    }

    #if 0
    double sample_theta_mj1( double x, struct dphmc_URandomState * uRandom ) {
        double u;
        URANDOM(u);
    }
    #endif
};

#if 0
/** Double accept-reject scheme is applied here to generate the A' emission
 * parameters */
int
dphmc_aprime_ww_sample_emission( struct dphmc_URandomState * uRandom
                               , double E0
                               , double ma
                               , double * xPtr
                               , double * thetaPtr ) {
    int rc;
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
            mj2 = (1/(2*E02*ma2*x*x))/(mxx2 + mema2);
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
    return 0;
}
#endif

class MajorantTestingApp : public dphmc::test::Fixture {
public:
    std::string oDirName;
    int N, nSamples;

    MajorantTestingApp();

    virtual void init() override;
    void run();
    virtual void clear() override;

private:
    /// Integrand wrapper for sampler
    static double _mj1fx_integrand(double x, void * sampler_) {
        return reinterpret_cast<Sampler*>(sampler_)->mj1(x);
    }
};

MajorantTestingApp::MajorantTestingApp() : N(100), nSamples(1e6) {
    parameters().define( "E0", "Incident energy (GeV)", 8e1 );
    parameters().define( "ma", "Mass of A' particle (GeV)", 0.04 );
}

void
MajorantTestingApp::init() {
    dphmc::test::Fixture::init();
}

void
MajorantTestingApp::run() {
    assert( !oDirName.empty() );

    dphmc::test::Histogram<1> hst(N, 0, 1);
    double x, theta;
    Sampler sampler( parameters()["E0"]
                   , parameters()["ma"]
                   );
    for( int i = 0; i < nSamples; ++i ) {
        hst.fill( sampler.sample_x_mj1(&_rgs) );
    }


    std::ofstream mjxf( oDirName + "/mjx.dat" )
                // ...
                ;

    dphmc::IterativeQAGSParameters mj1xInt = { 1e-8, 1e-8, 1.1, 1000, 1000 };
    gsl_integration_workspace * mj1xIntWSPtr
        = gsl_integration_workspace_alloc( mj1xInt.nnodes );

    double mj1xNInt, mj1xNIntAbsErr, mj1xNIntRelErr, mj1NInt;
    mj1NInt = dphmc_QAGS_integrate_iteratively( &mj1xInt
                                , _mj1fx_integrand
                                , &sampler
                                , 0, 1
                                , &mj1xNIntRelErr, &mj1xNIntAbsErr
                                , mj1xIntWSPtr );

    mjxf << "# x-sampling of M_{1,x} by M_{2,x};" << std::endl
         << "# E_0 = " << sampler.E0 << " GeV, m_a = "
            << sampler.ma << " GeV, n=" << nSamples << ";" << std::endl
         << "# underflow = " << hst.underflow()
            << ", overflow = " << hst.overflow() << ";" << std::endl
         << "# Integral value of Mj2 = " << sampler.mj2Norm << std::endl
         << "# x, M_{1,x}, M_{2,x}, w, \\integrate_{x_i}^{x_{i+1}}(M_{1,x}), relErr, absErr ;" << std::endl
         ;
    double mj2SumCheck = 0.;
    double xChi2 = 0.;
    const int nBins = hst.axis().nBins;
    for( int i = 0; i < nBins; ++i ) {
        double x = ((double) i)/nBins
             , nSimRatio = ((double) hst.bins()[i])/hst.sum()  // binned x-samples
             //, nSimRatio = (rand()/((double) RAND_MAX))/hst.sum()  // to assure chi^2 broken
             , deviation  // deviation to compute \chi^2
             ;
        mj1xNInt = dphmc_QAGS_integrate_iteratively( &mj1xInt
                                , _mj1fx_integrand
                                , &sampler
                                , x, x + 1./nBins
                                , &mj1xNIntRelErr, &mj1xNIntAbsErr
                                , mj1xIntWSPtr ) / mj1NInt;
        deviation = nSimRatio - mj1xNInt;
        xChi2 += deviation*deviation/mj1xNInt;
        mjxf << std::scientific
            << std::setw(14) << x << " "
            << std::setw(14) << sampler.mj1(x)/(nBins*sampler.mj2Norm) << " "
            << std::setw(14) << sampler.mj2(x)/(nBins*sampler.mj2Norm) << " "
            << std::setw(14) << nSimRatio << " "
            << std::setw(14) << mj1xNInt << " "
            << std::setw(14) << mj1xNIntRelErr << " "
            << std::setw(14) << mj1xNIntAbsErr << " "
            << std::endl;
    }
    xChi2 *= hst.sum();
    // To accept hypothesis of the values being distributed wrt given law
    // the calculated chi^2 must be < ch_{crit}^2
    std::cout << "chi^2/ndf for x: " << xChi2
              << " / " << (nBins-1) << " = " << xChi2/(nBins-1)
              << "; \\chi_{crit}(\\alpha=.05,ndf=" << nBins - 1
              << ") = " << gsl_cdf_chisq_Qinv(0.05, nBins - 1)
              << ", \\chi_{crit}(\\alpha=.01,ndf=" << nBins - 1
              << ") = " << gsl_cdf_chisq_Qinv(0.01, nBins - 1)
              << std::endl;
    gsl_integration_workspace_free( mj1xIntWSPtr );
    // ...
}

void
MajorantTestingApp::clear() {
    dphmc::test::Fixture::clear();
}

static void
usage_info( std::ostream & os
          , const std::string & appName
          , MajorantTestingApp & app ) {
    os << "Usage:" << std::endl
       << "    $ " << appName << " OPTIONS <output-dir>" << std::endl
       << "Where OPTIONS are:" << std::endl
       << " -N -- number of points (bins) on the plots. Default is 100."
          << std::endl
       << " -n,--n-samples -- number of samples to generate. Default is 1e6"
          << std::endl
       << "Application samples the WW cross section values for certain physics"
       << " and writes the resulting histograms and curves as dedicated data"
          " files into provided the directory. Files:" << std::endl
       << " ./mj2x.dat " << std::endl
       << " ./mj1x.dat " << std::endl
       ;
    os << "Parameters"
          " available for modification by `-D` option:" << std::endl;
    app.parameters().print_references(os);
    os << std::endl;
}

static int
configure_app( int argc, char * argv[]
             , MajorantTestingApp & app ) {

    static struct option longOpts[] = {
        { "help",               no_argument,        NULL, 'h' },
        { "set",                optional_argument,  NULL, 'D' },
        { "dir",                required_argument,  NULL, 'o' },
        { "n-samples",          required_argument,  NULL, 'n' },
        // ...
        { NULL, 0x0, NULL, 0x0 },
    };

    int c, optIdx;
    bool hadError = false
       , listRegistered = false;
    // Iterate over command line options to set fields in appCfg
    while((c = getopt_long( argc, argv, "hD:N:"
                          , longOpts, &optIdx )) != -1) {
        switch(c) {
            case '0' :
                std::cerr << "Unable to parse option \""
                          << longOpts[optIdx].name << "\""  << std::endl;
                break;
            case 'd' :
                if( 0 != app.parse_parameter_setting(optarg) ) {
                    std::cerr << "Unable to parse \"" << optarg << "\"" << std::endl;
                    return -1;
                }
                break;
            case 'N' :
                app.N = atoi(optarg);
                break;
            case 'n' :
                app.nSamples = atoi(optarg);
                break;
            case 'o':
                app.oDirName = optarg;
                break;
            // ... more runtime options may be added here
            case 'h' :
                usage_info( std::cout, argv[0], app );
                return 1;
            case ':' :
                std::cerr << "Option -" << optopt << " requires an argument."
                          << std::endl;
                hadError = true;
                break;
            case '?' :
                std::cerr << "Unrecognized option '-" << optopt << "'."
                          << std::endl;
                hadError = true;
        }
    }
    if( hadError )
        return -1;
    if( optind >= argc ) {
        std::cerr << "Error: expected single positional argument: output dir path."
                  << std::endl;
        return -1;
    }
    app.oDirName = argv[optind];
    return 0;
}




int
main(int argc, char * argv[]) {
    MajorantTestingApp app;
    int rc = configure_app( argc, argv, app );
    if( rc < 0 ) {
        std::cerr << "An error occured during command line arguments parsing."
                  << std::endl;
        usage_info( std::cerr, argv[0], app );
        return EXIT_FAILURE;
    } else if( rc > 0 ) {
        return EXIT_SUCCESS;
    }
    app.init();
    app.run();
    app.clear();
    return EXIT_SUCCESS;
}
