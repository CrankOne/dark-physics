# ifndef H_DPHMC_APRIME_WW_TEST_MJ_FICTURE_H
# define H_DPHMC_APRIME_WW_TEST_MJ_FICTURE_H

# include "models/ww-approx.h"
# include "dphmc-rnd.h"
# include "dphmc-test-app.hh"

namespace dphmc {
namespace test {

class Fixture {
protected:
    ParametersSet _ps;
    URandomState _rgs;
public:
    /// Cosntructs the basic fixture instance with random number generator
    Fixture();

    /// Initializes random number generator
    virtual void init();
    /// Frees fixture's resources
    virtual void clear();

    /// Returns parameters instance
    ParametersSet & parameters() { return _ps; }
    /// Returns parameters instance (ro)
    const ParametersSet & parameters() const { return _ps; }

    /// Returns generator instance
    URandomState & rng() { return _rgs; }
    /// Returns generator instance (cosnt)
    const URandomState & rng() const { return _rgs; }

    /// Used by applications; parses parameter setting of the form <name=val>
    int parse_parameter_setting(const std::string &);
};

class APrimePhysFixture : public Fixture {
protected:
    struct dphmc_APrimePhysParameters _phPars;
public:
    APrimePhysFixture();

    /// Initializes physics parameters struct
    virtual void init() override;
};

class APrimeChiIntegrFixture : public APrimePhysFixture {
protected:
    struct dphmc_IterativeQAGSParameters _integrPars;
    struct dphmc_APrimeWWCaches * _cachesPtr;
public:
    APrimeChiIntegrFixture();

    /// Initializes chi integrating parameters struct
    virtual void init() override;
    /// Deletes cache
    virtual void clear() override;
};

}
}

# endif
