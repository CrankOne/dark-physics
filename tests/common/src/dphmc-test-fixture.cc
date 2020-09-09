# include "dphmc-test-fixture.hh"

namespace dphmc {
namespace test {

Fixture::Fixture() {
    define_rng_parameters(_ps);
}

void
Fixture::read_parameters() {
    configure_rng( _ps, _rgs );
}

void
Fixture::free_resources() {
    dphmc_rnd_gen_gsl_free( &_rgs );
}


APrimePhysFixture::APrimePhysFixture() {
    define_aprime_particle_parameters(_ps);
}

void
APrimePhysFixture::read_parameters() {
    Fixture::read_parameters();
    configure_phys_parameters( parameters(), _phPars );
}


APrimeChiIntegrFixture::APrimeChiIntegrFixture() {
    define_media_parameters(_ps);
    define_cross_section_parameters(_ps);
}

void
APrimeChiIntegrFixture::read_parameters() {
    APrimePhysFixture::read_parameters();
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
APrimeChiIntegrFixture::free_resources() {
    APrimePhysFixture::free_resources();
    if( _cachesPtr ) {
        // clear parametric caches
        dphmc_aprime_delete( _cachesPtr );
        _cachesPtr = nullptr;
    }
}

}
}

