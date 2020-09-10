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
    double sample_x_mj1( struct dphmc_URandomState * uRandom ) const {
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

    /**For given \f$x\f$ returns random \f$theta\f$ distributed wrt
     * \f$M_{1}(x, \theta)\f$:
     * 
     * \f[
     * \theta = \pi \sqrt{\frac{ u \times ((1-x)/x^2 + m_{e}^2/m_{A'}^2) }
     *          { (1-x)/x^2 + m_{e}^2/m_{A'}^2 + E_0^2/m_{A'} \pi^2 (u-1) }}
     * \f]
     * */
    double sample_theta_mj1( double x, struct dphmc_URandomState * uRandom ) const {
        int rc;
        const double mxx2 = (1-x)/(x*x);
        double u;
        URANDOM(u);
        return M_PI*sqrt( u*(mxx2+mema2) / (mxx2 + mema2 + E02ma2*pi2*(1-u)) );
    }

    /** Returns value of \f$M_{1}(x, \theta)\f$ */
    double mj1(double x, double theta) const {
        const double mxx2 = (1-x)/(x*x);
        return x*theta / (E02ma2 * theta * theta + mxx2 + mema2 );
    }

    double reference_f( double x, double theta ) const {
        const double xm1 = 1 - x
                   , U = E02*theta*theta*x + ma2*xm1/x + me*me*x
                   , U2 = U*U
                   ;
        return ((1 - x + x*x/2)/U2 + ma2*(
                        xm1*xm1*ma2 - xm1*U*x
                    )/(U2*U2))*x*sin(theta);
    }

    /** Samples point and returns reference value at this point */
    double sample_x_theta( struct dphmc_URandomState * uRandom
                         , double & xRef
                         , double & thetaRef ) const {
        int rc;
        // Sample X according to 2nd majorant (M_{2,x}(x))
        double x, theta, probe, reference;
        do {
            // generate a candidate (triplet)
            x = sample_x_mj1(uRandom);
            theta = sample_theta_mj1(x, uRandom);
            URANDOM(probe);
            // test triplet
            reference = reference_f(x, theta);
            if( true /*probe*mj1(x, theta) <= reference*/ ) {  // TODO
                // within a target function -- accept:
                xRef = x;
                thetaRef = theta;
                return reference;
            }
        } while(true);
    }
};

class MajorantTestingApp : public dphmc::test::Fixture {
public:
    std::string mj2File  ///< when non-empty, the mj2 .dat will be generated
              , mj1File  ///< when non-empty, the mj1 .dat will be generated
              ;
    int nx, ny, nSamples;

    MajorantTestingApp();

    virtual void init() override;
    void run();
    void run_mj1(const Sampler &);
    void run_mj2(const Sampler &);
    virtual void clear() override;

private:
    /// Integrand wrapper for sampler
    static double _mj1fx_integrand(double x, void * sampler_) {
        return reinterpret_cast<Sampler*>(sampler_)->mj1(x);
    }
};

MajorantTestingApp::MajorantTestingApp() : nx(100), ny(100), nSamples(1e6) {
    parameters().define( "E0", "Incident energy (GeV)", 8e1 );
    parameters().define( "ma", "Mass of A' particle (GeV)", 0.04 );
}

void
MajorantTestingApp::init() {
    dphmc::test::Fixture::init();
}

void
MajorantTestingApp::run() {
    Sampler sampler( parameters()["E0"]
                   , parameters()["ma"]
                   );
    if( !mj1File.empty() ) run_mj1(sampler);
    if( !mj2File.empty() ) run_mj2(sampler);
}

void
MajorantTestingApp::run_mj2(const Sampler & sampler) {
    dphmc::test::Histogram<1> hst(nx, 0, 1);
    double x;
    for( int i = 0; i < nSamples; ++i ) {
        hst.fill( sampler.sample_x_mj1(&_rgs) );
    }

    std::ofstream mjxf( mj2File );

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
}

void
MajorantTestingApp::run_mj1(const Sampler & sampler) {
    dphmc::test::Histogram<2> hst( nx, 0, 1
                                 , ny, 0, M_PI
                                 );
    double x, theta;
    for( int i = 0; i < nSamples; ++i ) {
        sampler.sample_x_theta(&_rgs, x, theta);
        //std::cout << x << ", " << theta << std::endl;  // XXX
        hst.fill( x, theta );
    }

    const size_t n1Bins = hst.axis1().nBins
               , n2Bins = hst.axis2().nBins
               ;

    std::ostream & mjf = std::cout;  //std::ofstream(mj1File);  // TODO

    for( size_t i = 0; i < n1Bins; ++i ) {
        for( size_t j = 0; j < n2Bins; ++j ) {
            mjf << hst(i,j) << " ";
        }
        mjf << std::endl;
    }
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
       << "    $ " << appName << " OPTIONS <output-file>" << std::endl
       << "Where OPTIONS are:" << std::endl
       << " -N -- number of points (bins) on the plots. Default is 100."
          " Prefix number with 'x'/'t' to indicate dimension: x or theta."
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
    std::string oFile;
    bool eval1D = false;
    // ^^^ NOTE: 1D is for Majorant #2, 2D is for Majorant 1D... kind of tricky

    // Iterate over command line options to set fields in appCfg
    while((c = getopt_long( argc, argv, "h2D:N:n:"
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
                if( isdigit(optarg[0]) ) {
                    app.nx = atoi(optarg);
                } else if( 'x' == optarg[0] ) {
                    app.nx = atoi(optarg+1);
                } else if( 't' == optarg[0] ) {
                    app.ny = atoi(optarg+1);
                } else {
                    std::cerr << "Bad first character at cmd-line token: \""
                              << optarg << "\"."
                              << std::endl;
                    return -1;
                }
                break;
            case 'n' :
                app.nSamples = atoi(optarg);
                break;
            case '2':
                eval1D = true;
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
        std::cerr << "Error: expected single positional argument: output file path."
                  << std::endl;
        return -1;
    }
    if( eval1D )
        app.mj2File = argv[optind];
    else
        app.mj1File = argv[optind];
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
