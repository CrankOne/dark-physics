# include "models/ww-approx.hh"

# include <G4ParticleDefinition.hh>
# include <G4Element.hh>
# include <G4Electron.hh>

# include <cassert>
# include <limits>

namespace dphmc {

bool
APrimeWWApproximation::GeneratorKey::operator==(
        const APrimeWWApproximation::GeneratorKey & other ) const {
    return materialZ == other.materialZ
        && energyRange == other.energyRange
        && pType == other.pType
        ;
}

::std::size_t
APrimeWWApproximation::HashGeneratorKey::operator()(
        const APrimeWWApproximation::GeneratorKey & k) const noexcept {
    return  std::hash<uint16_t>{}(k.materialZ)
         ^ (std::hash<uint16_t>{}(k.energyRange) << 1) >> 1
         ^ (std::hash<G4int>{}(k.pType) << 1)
         ;
}

/** Constructs the model for given projectile energies. If the number of energy
 * ranges bounds supplied via `incidentEs` is <2, throws a `std::runtime_error`
 * exception.
 *
 * @param lowestZ Lowest element charge number. Lighter elements won't be considered
 * @param incidentEs Energy ranges bounds (may be unsorted, must be terminated with zero)
 * */
APrimeWWApproximation::APrimeWWApproximation( const G4double lowestZ
                                            , G4double * incidentEs ) : fLowestZ(lowestZ) {
    assert( incidentEs );
    while( *incidentEs ) {
        fRanges.insert( *(incidentEs++) );
    }
    if( fRanges.size() < 2 ) {
        throw std::runtime_error( "Need at least two energy range bounds." );
    }
}

APrimeWWApproximation::GeneratorsIndex::iterator
APrimeWWApproximation::RetrieveGeneratorFor( G4double incidentParticleEnergy
                                           , const G4ParticleDefinition * pDef
                                           , const G4Element * g4Element ) {
    if( *fRanges.begin() > incidentParticleEnergy
     || *fRanges.rbegin() < incidentParticleEnergy
     || g4Element->GetZ() < fLowestZ ) {
        // return end iterator signaling that we are out of ranges
        return fGenerators.end();
    }
    assert( ! fRanges.empty() );  // must be filled to the time
    auto genRangeIt = fRanges.lower_bound(incidentParticleEnergy);
    GeneratorKey key { g4Element->GetZ()
                     , std::distance( fRanges.begin()
                                    , genRangeIt )
                     , pDef->GetPDGEncoding()
                     };
    auto genIt = fGenerators.find( key );
    if( fGenerators.end() == genIt ) {
        // Instantiate and register new generator instance
        //auto range = fRanges.equal_range( incidentParticleEnergy );
        auto ub = fRanges.upper_bound( incidentParticleEnergy )
           , lb = ub;
        assert( ub != fRanges.begin() );  // discriminated by cond >= E_min above
        assert( ub != fRanges.end() );  // discriminated by cond <= E_max above
        --lb;
        assert( *lb < *ub );
        auto ir = fGenerators.emplace( key, NewGenerator( *lb, *ub
                                                        , pDef, g4Element ));
        assert( ir.second );
        genIt = ir.first;
    }
    return genIt;
}

/** Chooses appropriate generator instance and forwards the call to its
 * AbstractGenerator::GetFullCrossSection() */
G4double
APrimeWWApproximation::GetFullCrossSection( G4double incidentParticleEnergy
                                          , const G4ParticleDefinition * pDef
                                          , const G4Element * g4Element ) {
    auto genIt = RetrieveGeneratorFor( incidentParticleEnergy, pDef, g4Element );
    if( fGenerators.end() == genIt ) {
        // out of range of model's validity
        return std::numeric_limits<G4double>::has_quiet_NaN ? std::nan("1") : 0.;
    }
    return genIt->second->GetFullCrossSection();
}

/** Chooses appropriate generator instance and forwards the call to its
 * AbstractGenerator::ShootKinematics() */
G4double
APrimeWWApproximation::GenerateOn( G4double incidentParticleEnergy
                                 , const G4ParticleDefinition * pDef
                                 , const G4Element * g4Element
                                 , G4double * energy
                                 , G4double * theta
                                 , G4double * phi ) {
    auto genIt = RetrieveGeneratorFor( incidentParticleEnergy, pDef, g4Element );
    if( fGenerators.end() == genIt ) {
        // out of range of model's validity
        return std::numeric_limits<G4double>::has_quiet_NaN ? std::nan("1") : 0.;
    }
    *energy = incidentParticleEnergy;
    return genIt->second->ShootKinematics( incidentParticleEnergy, energy, theta, phi );
}

/**Default model implementation is valid only for electrons. */
G4bool
APrimeWWApproximation::IsApplicableToParticle( const G4ParticleDefinition & pDef ) {
    return G4Electron::Definition()->GetPDGEncoding() == pDef.GetPDGEncoding();
}

APrimeWWApproximation::AbstractGenerator *
APrimeWWApproximation::NewGenerator( G4double projELow, G4double projEUp
                                   , const G4ParticleDefinition * pDef
                                   , const G4Element * g4Element ) {
    // TODO: instantiate dumb von Nemann generator implementation
    throw std::runtime_error( "Not yet implemented." );
}

}
