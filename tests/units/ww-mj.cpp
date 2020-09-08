# include "models/ww-approx.h"
# include "dphmc-rnd.h"

# include <gtest/gtest.h>

class MajorantTest : public ::testing::Test {
protected:
    struct dphmc_IterativeQAGSParameters _integrPars;
    struct dphmc_APrimePhysParameters _phPars;
    struct dphmc_URandomState _rgs;
    struct dphmc_APrimeWWCaches * _cachesPtr;
public:
    MajorantTest() : _cachesPtr(nullptr) {
        // Physical parameters for tests
        // WARNING: hardcoded static data in UTs were calculated using these
        // parameters!
        _phPars.Z = 74;
        _phPars.A = 183.84;
        _phPars.massA_GeV = .1;
        _phPars.EBeam_GeV = 80;
        _phPars.epsilon = 1e-4;
        _phPars.factor = 1;
        _phPars.maxThetaFactor = 1;
        // Some default integration workspace parameters for \chi calculus
        _integrPars.epsabs = 1e-12;
        _integrPars.epsrel = 1e-12;
        _integrPars.epsrelIncFt = 1.1;
        _integrPars.limit = 1e3;
        _integrPars.nnodes = 1e3;
        // Initialize random number generator
        dphmc_rnd_gen_gsl_init( &_rgs
                              , gsl_rng_ranlux
                              , 1337  // seed
                              );
    }
    virtual void SetUp() override {
        _cachesPtr = dphmc_aprime_new( &_phPars
                                     , &_integrPars
                                     , dphmc_aprime_form_factor_elastic
                                     , dphmc_aprime_form_factor_inelastic
                                     , 0x0 );
    }
    virtual void TearDown() override {
        // clear parametric caches
        dphmc_aprime_delete( _cachesPtr );
        _cachesPtr = nullptr;
        // clear the random number generator
        dphmc_rnd_gen_gsl_free( &_rgs );
    }
};

// Assure the x sampling of M_{x, 1} produces resonable results
TEST_F( MajorantTest, x_chiSq_test ) {
    size_t nSamples = 1e7;
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

    for( int i = 0; i < nBins; ++i ) {
        std::cout << i << " " << histogram[i]/((double) nSamples) << std::endl;
    }
}

#if 0
TEST_F( MajorantTest, majorant2samples_majorant1 ) {
    // ...
}
#endif
