# include "models/ww-approx.hh"

# include <G4Element.hh>
# include <G4NistManager.hh>
# include <Randomize.hh>
# include <G4Electron.hh>

# include <gtest/gtest.h>

namespace dphmc {
namespace test {

static G4double gMinZInMockModel = 40;
// Values are shuffled intentionally, 3 ranges are defined here:
//      4..12, 12..74, 74..123
static G4double gEnergyRanges[] = { 123, 12, 5, 74, 0 }
              , gEnergyRangesMin = 5
              , gEnergyRangesMax = 123
              ;

/// Special "mock" (pseudo) model to test indexing correctnesss
class MockModel : public APrimeWWApproximation {
public:
    MockModel( ) : APrimeWWApproximation( gMinZInMockModel, gEnergyRanges ) {}

    /// This mocking generator yields some unrealistic but easy-to-check
    /// values
    struct MockGenerator : AbstractGenerator {
        const G4double fullCS;

        const G4double thisERange[2];
        const G4ParticleDefinition * thisParticleDef;
        const G4Element * thisG4Element;

        MockGenerator( G4double fcs
                     , G4double eMin, G4double eMax
                     , const G4ParticleDefinition * pDef
                     , const G4Element * g4El
                     ) : fullCS(fcs)
                       , thisERange{eMin, eMax}
                       , thisParticleDef( pDef )
                       , thisG4Element( g4El )
                       {
            assert( eMin < eMax );
        }

        /// Returns fullCS value
        virtual G4double GetFullCrossSection() const override { return fullCS; }

        void AssureProjectileEnergyIsValid( G4double E ) const {
            ASSERT_GE( E, thisERange[0] );
            ASSERT_LT( E, thisERange[1] );
        }

        /// Returns uniformly distributed values.
        virtual G4double ShootKinematics( const G4double projectileEnergy
                                        , G4double * energy
                                        , G4double * theta
                                        , G4double * phi ) const override {
            AssureProjectileEnergyIsValid( projectileEnergy );
            *energy = G4RandFlat::shoot( thisERange[0], projectileEnergy );
            *theta = G4RandFlat::shoot(M_PI);
            *phi = G4RandFlat::shoot(0., 2.*M_PI);
            return G4RandFlat::shoot();
        }
    };

    // due to limitation of ASSERT_* macro this check has to be in
    // void-returning functions.
    void CheckNewGenCallIsValid( G4double projELow, G4double projEUp
                               , const G4ParticleDefinition * pDef
                               , const G4Element * g4Element ) const {
        ASSERT_GE( g4Element->GetZ(), gMinZInMockModel );
        ASSERT_GT( projEUp, projELow );
        ASSERT_GE( projELow, gEnergyRangesMin );
        ASSERT_LE( projEUp, gEnergyRangesMax );
    }

    /// Creates `MockGenerator` instance with `fullCS` member 
    virtual AbstractGenerator * NewGenerator( G4double projELow, G4double projEUp
                                            , const G4ParticleDefinition * pDef
                                            , const G4Element * g4Element ) override {
        CheckNewGenCallIsValid( projELow, projEUp, pDef, g4Element );
        return new MockGenerator( g4Element->GetA()/(projEUp - projELow)
                                , projELow, projEUp, pDef, g4Element );
    }

    /// Extends mock model to any charged particle
    virtual G4bool IsApplicableToParticle( const G4ParticleDefinition & pDef ) override {
        if( !pDef.GetPDGCharge() ) {
            return false;
        }
        return true;
    }
};

}  // namespace ::dphmc::test
}  // namespace dphmc

TEST( G4Infra, modelIndexing ) {
    dphmc::test::MockModel * mPtr = new dphmc::test::MockModel();

    ASSERT_FALSE( mPtr->GetEnergyTabulation().empty() );
    
    // Get some elements for testing from NIST materials DB
    G4NistManager* nistManager = G4NistManager::Instance();
    G4bool fromIsotopes = false;
    G4Element * elPb = nistManager->FindOrBuildElement("Pb", false)
            , * elAu = nistManager->FindOrBuildElement("Au", false)
            , * elXe = nistManager->FindOrBuildElement("Xe", false)
            , * elC =  nistManager->FindOrBuildElement("C", false)
            , * els[] = { elPb, elAu, elXe, elC, nullptr }
            ;

    // Generate some data using mock generator
    for( size_t nEv = 0; nEv < 1e5; ++nEv ) {
        // pick an element
        G4Element * choosenEl = els[G4RandFlat::shootInt( sizeof(els)/sizeof(G4Element*) - 1 )];
        assert(choosenEl);
        // Take some projectile
        G4ParticleDefinition * pDef = G4Electron::Definition();  // "e-"
        // choose an energy
        G4double projE = G4RandFlat::shoot( 0., 200. );

        // Get the "full cross section"
        G4double fullCS = mPtr->GetFullCrossSection( projE, pDef, choosenEl );

        bool isValidForConditions = true;
        // Require full cross section be NaN (or 0, if G4double type does not
        // support quiet NaN) for out of scope of definition and have some
        // value, defined by A/(E_{up} - E_{low}) otherwise
        if( projE < dphmc::test::gEnergyRangesMin
         || projE > dphmc::test::gEnergyRangesMax ) {
            ASSERT_TRUE( std::isnan(fullCS) );
            isValidForConditions = false;
        }
        if( choosenEl->GetZ() < dphmc::test::gMinZInMockModel ) {
            ASSERT_TRUE( std::isnan(fullCS) );
            isValidForConditions = false;
        }
        if( isValidForConditions ) {
            ASSERT_FALSE( std::isnan(fullCS) );
            ASSERT_TRUE( std::isfinite(fullCS) );
        }

        // Shoot event
        G4double E, theta, phi;
        G4double prob = mPtr->GenerateOn( projE, pDef, choosenEl, &E, &theta, &phi );
        if( ! isValidForConditions ) {
            ASSERT_TRUE( std::isnan(prob) );
            continue;
        }
        // naively get upper an lower energy limits from global list
        G4double eLow = *dphmc::test::gEnergyRanges
               ,  eUp = *dphmc::test::gEnergyRanges;
        for( size_t i = 0
           ; i < sizeof(dphmc::test::gEnergyRanges)/sizeof(double) - 1
           ; ++i ) {
            if( eLow > dphmc::test::gEnergyRanges[i] ) eLow = dphmc::test::gEnergyRanges[i];
            if(  eUp < dphmc::test::gEnergyRanges[i] )  eUp = dphmc::test::gEnergyRanges[i];
        }
        ASSERT_GT( E, eLow );
        ASSERT_LT( E, eUp );
        ASSERT_GT( theta, 0. );
        ASSERT_LT( theta, M_PI );
        ASSERT_GT( theta, 0. );
        ASSERT_LT( theta, 2*M_PI );
    }
    ASSERT_FALSE( mPtr->GetGeneratorsIndex().empty() );
}

