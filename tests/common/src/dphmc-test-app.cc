#include "dphmc-test-app.hh"

namespace dphmc {
namespace test {

void
ParametersSet::define( const std::string & name
                     , const std::string & description
                     , double value ) {
    auto ir = emplace( name, std::make_pair(value, description) );
    if( !ir.second ) {
        throw std::runtime_error( "Duplicating parameter name definition." );
    }
}

void
ParametersSet::print_references( std::ostream & os ) const {
    for( auto parPair : *this ) {
        os << "<" << parPair.first << "="
           << parPair.second.first << "> -- "
           << parPair.second.second << std::endl; 
    }
}

double &
ParametersSet::operator[](const std::string & parName) {
    auto it = find(parName);
    if( it == end() ) {
        throw std::runtime_error( "Unknown parameter name." );
    }
    return it->second.first;
}

double
ParametersSet::operator[](const std::string & parName) const {
    auto it = find(parName);
    if( it == end() ) {
        throw std::runtime_error( "Unknown parameter name." );
    }
    return it->second.first;
}

void
define_aprime_particle_parameters( ParametersSet & ps ) {
    ps.define( "ma"
             , "for A' physics: mass of the A' particle, GeV"
             , .004 );
    ps.define( "E0"
             , "for A' physics: energy of the incident particle, GeV"
             , 80 );
}

void
define_media_parameters( ParametersSet & ps ) {
    ps.define( "Z"
             , "Z number of the target media"
             , 74 );
    ps.define( "A"
             , "A number of the target nuclei"
             , 183.84 );
}

void
define_cross_section_parameters( ParametersSet & ps ) {
    ps.define( "epsilon"
             , "mixing constant"
             , 1e-4 );
    ps.define( "chiEpsAbs"
             , "effective photon flux N-integration absolute error"
             , 1e-12 );
    ps.define( "chiEpsRel"
             , "effective photon flux N-integration relative error"
             , 1e-12 );
    ps.define( "chiEpsRelInc"
             , "effective photon flux integration increment factor (must be >1)"
             , 1.1 );
    ps.define( "chiLimit"
             , "Limiting number of iterations for numerical"
               " integration of effective photon flux"
             , 1e3 );
    ps.define( "chiNNodes"
             , "Limiting number of iterations for numerical"
               " integration ???"
             , 1e3 );
}

void
configure_phys_parameters( const ParametersSet & ps
                         , aprime::PhysParameters & pp ) {
    pp.Z = ps["Z"];
    pp.A = ps["A"];
    pp.massA_GeV = ps["ma"];
    pp.EBeam_GeV = ps["E0"];
    pp.epsilon = ps["epsilon"];
    // pp.factor =
    // pp.maxThetaFactor =
}

void
configure_chiint_parameters( const ParametersSet & ps
                           , IterativeQAGSParameters & ip ) {
    ip.epsabs = ps["chiEpsAbs"];
    ip.epsrel = ps["chiEpsRel"];
    ip.epsrelIncFt = ps["chiEpsRelInc"];
    ip.limit = ps["chiLimit"];
    ip.nnodes = ps["chiNNodes"];
}


}
}

