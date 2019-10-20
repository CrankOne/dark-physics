# include "models/ww-approx.hh"

# include <G4ParticleDefinition.hh>

# include <cassert>

namespace dphmc {

bool
APrimeWWApproximation::GeneratorKey::operator==(
        const APrimeWWApproximation::GeneratorKey & other ) const {
    return materialZ == other.materialZ
        && energyRange == other.energyRange
        // && ...
        ;
}

::std::size_t
APrimeWWApproximation::HashGeneratorKey::operator()(
        const APrimeWWApproximation::GeneratorKey & k) const noexcept {
    return  std::hash<uint16_t>{}(k.materialZ)
         ^ (std::hash<uint16_t>{}(k.energyRange) << 1) >> 1
         // ^ ...
         ;
}

/** Constructs the model for given projectile energies.
 *
 * @param lowestZ Lowest element charge number. Lighter elements won't be considered
 * @param incidentEs Energy ranges bounds (may be unsorted, must be terminated with zero)
 * */
APrimeWWApproximation::APrimeWWApproximation( const G4double lowestZ
                                            , G4double * incidentEs ) {
    assert( incidentEs );
    for( G4double * prjE = incidentEs; *prjE; ++prjE ) {
        fRanges.insert( *prjE );
    }
}

G4double
APrimeWWApproximation::GetFullCrossSection( G4double incidentParticleEnergy
                                          , const G4ParticleDefinition * pDef
                                          , const G4Element * g4Element ) {
}

G4double
APrimeWWApproximation::GenerateOn( G4double incidentParticleEnergy
                                 , const G4ParticleDefinition * pDef
                                 , const G4Element * g4Element
                                 , G4double * energy
                                 , G4double * theta
                                 , G4double * phi ) {
}

/**@todo currently, returns `true` for all the charged particles, but model
 * is, actually, valid only for electrons. Consider to extend the application
 * to any charged lepton particle.
 * */
G4bool
APrimeWWApproximation::IsApplicableToParticle( const G4ParticleDefinition & pDef ) {
    if( !pDef.GetPDGCharge() ) {
        return false;
    }
    return true;
}

}

