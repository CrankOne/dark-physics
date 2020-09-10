# include "dphmc-test-fixture.hh"

namespace dphmc {
namespace test {

Fixture::Fixture() {
    define_rng_parameters(_ps);
}

void
Fixture::init() {
    configure_rng( _ps, _rgs );
}

void
Fixture::clear() {
    dphmc_rnd_gen_gsl_free( &_rgs );
}

int
Fixture::parse_parameter_setting(const std::string & strv) {
    size_t eqPos = strv.find('=');
    if( std::string::npos == eqPos ) return -1;  // no '=' sign
    std::string parName = strv.substr(0, eqPos);
    if(parName.empty()) return -2;  // var name is empty
    std::string parValStr = strv.substr( eqPos+1, std::string::npos );
    if(parValStr.empty()) return -3;  // value string is empty
    char * endptr = nullptr;
    double val = strtod(parValStr.c_str(), &endptr);
    if( parValStr.c_str() == endptr
     || (endptr && *endptr != '\0') ) return -4;  // conversion not full or failed
    try {
        _ps[parName] = val;
    } catch(std::runtime_error & e) {
        return -5;  // failed to set value
    }
    return 0;
}

APrimePhysFixture::APrimePhysFixture() {
    define_aprime_particle_parameters(_ps);
}

void
APrimePhysFixture::init() {
    Fixture::init();
    configure_phys_parameters( parameters(), _phPars );
}


APrimeChiIntegrFixture::APrimeChiIntegrFixture() {
    define_media_parameters(_ps);
    define_cross_section_parameters(_ps);
}

void
APrimeChiIntegrFixture::init() {
    APrimePhysFixture::init();
    configure_chiint_parameters( parameters(), _integrPars );
    _cachesPtr = dphmc_aprime_new( &_phPars
                                 , &_integrPars
                                 , dphmc_aprime_form_factor_elastic  // TODO: parameterize
                                 , dphmc_aprime_form_factor_inelastic  // TODO: parameterize
                                 , 0x0  // TODO: parameterize
                                 );
    // ...
}

void
APrimeChiIntegrFixture::clear() {
    APrimePhysFixture::clear();
    if( _cachesPtr ) {
        // clear parametric caches
        dphmc_aprime_delete( _cachesPtr );
        _cachesPtr = nullptr;
    }
}

}
}

