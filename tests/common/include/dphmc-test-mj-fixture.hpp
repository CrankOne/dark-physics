# ifndef H_DPHMC_APRIME_WW_TEST_MJ_FICTURE_H
# define H_DPHMC_APRIME_WW_TEST_MJ_FICTURE_H

# include "models/ww-approx.h"
# include "dphmc-rnd.h"

class MajorantTestCommon {
protected:
    struct dphmc_IterativeQAGSParameters _integrPars;
    struct dphmc_APrimePhysParameters _phPars;
    struct dphmc_URandomState _rgs;
    struct dphmc_APrimeWWCaches * _cachesPtr;
public:
    MajorantTestCommon() : _cachesPtr(nullptr) {
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
    virtual void init_caches() {
        _cachesPtr = dphmc_aprime_new( &_phPars
                                     , &_integrPars
                                     , dphmc_aprime_form_factor_elastic
                                     , dphmc_aprime_form_factor_inelastic
                                     , 0x0 );
    }
    virtual void free_caches() {
        // clear parametric caches
        dphmc_aprime_delete( _cachesPtr );
        _cachesPtr = nullptr;
        // clear the random number generator
        dphmc_rnd_gen_gsl_free( &_rgs );
    }
};

# endif
