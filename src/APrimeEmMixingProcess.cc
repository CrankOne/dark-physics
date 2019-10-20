# include <Randomize.hh>

# include <limits>

# include "APrimeEmMixingProcess.hh"
# include "APrimeScatteringModel.hh"

// TODO: use CLHEP's RandGeneral instead of TFoam
// https://proj-clhep.web.cern.ch/proj-clhep/doc/CLHEP_2_2_0_4/html/RandGeneral_8h_source.html

namespace dphmc {

APrimeEmMixingProcess::APrimeEmMixingProcess( APrimeScatteringModel * modelPtr )
            : G4VDiscreteProcess( "WW-based A' EM mixing", fUserDefined )
            , theAPrimePtr( dphmc::APrime::Definition() )
            , fModelPtr(modelPtr) {
    SetProcessSubType( 83 );  // TODO: correct this
}

/** Delegates call to current model's
 * APrimeScatteringModel::IsApplicableToParticle().
 * */
G4bool
APrimeEmMixingProcess::IsApplicable(const G4ParticleDefinition & pDef) {
    return fModelPtr->IsApplicableToParticle(pDef);
}

/** Calculates mean free path \f$\lambda (E)\f$ based on material compound and
 * an incident particle type and energy.
 *
 * The mean free path \f$\lambda\f$ for \f$n\f$ elements with relative
 * fractions \f$n_i, i = 0 ... n\f$ is:
 * \f[
 *      \lambda (E) = 1 / \sum\limits_0^n [ n_i \cdot \sigma(E, Z_i) ],
 * \f]
 * where \f$Z_i\f$ stands for charge number of element \f$i\f$. The material
 * properties (charge and abundance) is supplied with Geant4 stepping
 * information, the full (integrated) cross-section \f$sigma(E, Z_i)\f$ is a
 * subject of a particular model.
 *
 * @todo When no elements were appliable for the media, particle or energy,
 * returns a `DBL_MAX` value. It is not clear, whether it is legal for Geant4
 * API (or, for instance, `inf` or zero have to be returned).
 * @todo `condition` is always set to `NotForced` currently.
 */
G4double
APrimeEmMixingProcess::GetMeanFreePath( const G4Track & aTrack,
                                 G4double /*previousStepSize*/,
                                 G4ForceCondition * condition ) {
    G4double lambdaInv = 0.;

    // not applicable for neutral particles; this shall be discriminated by the
    // Geant4 API before, so triggering of this assertion evidents the
    // malformed process definition
    assert( aTrack.GetDefinition()->GetPDGCharge() );

    // TODO: consider forced case
    *condition = NotForced;

    // Get the energy of incident particle
    const G4double incidentE = aTrack.GetKineticEnergy();
    // Get the lements list
    const G4ElementVector & elList = *( aTrack.GetMaterial()->GetElementVector() );
    // For each element of the current media, compute the mean free path
    size_t nEl = 0;
    for( auto it = elList.cbegin(); elList.cend() != it; ++it, ++nEl ) {
        G4Element & el = **it;
        G4double val = fModelPtr->GetFullCrossSection( incidentE, aTrack.GetDefinition(), &el );
        if( (std::numeric_limits<G4double>::has_quiet_NaN && std::isnan( val ))
         || val <= 0. )
            continue;  // not appliable for this element
        val *= aTrack.GetMaterial()->GetVecNbOfAtomsPerVolume()[nEl];
        lambdaInv += val;
    }
    if( lambdaInv ) return 1/lambdaInv;
    // None chemical elements were appliable for the process
    return DBL_MAX;  // TODO: check that it is legal
}


/** Routes the A' production among generators associated with materials in
 * step's media. Returned object is a proposal for Geant4 tracking management
 * algorithm to produce new particle and change momenta of the incident
 * projectile -- an instance of G4DynamicParticle class allocated on heap. The
 * projectile recoil is calculated by the following formulas:
 *
 * \f[
 * E_R = (1-x) E_0,
 * \f]
 * for the new energy of leaving projectile and
 * \f[
 * \theta_R = \arctan( \sqrt{ \frac{m_{A'}}{E_0} } \cdot ( 1 + \frac{m_{A'}}{2 E_0} ) ), \quad
 * \phi_R = - \phi_{A'},
 * \f]
 * for its angles in lab frame.
 *
 * @todo More precise recoil kinematics. Now just took it from \cite{Bjorken}
 * (C3-C5), but some elaboration is possible as we have mostly elastic reaction.
 * @todo Make use of `prob` variable. Consider to move the kinematical calculus
 * out of this class (according to distributed logic of Geant4 API) .
 * */
G4VParticleChange *
APrimeEmMixingProcess::PostStepDoIt( const G4Track & aTrack
			                       , const G4Step & /*aStep*/ ) {
    // Get the energy of incident particle
    const G4double incidentE = aTrack.GetKineticEnergy();
    // Get the lements list
    const G4ElementVector & elList = *( aTrack.GetMaterial()->GetElementVector() );
    // Having integrated cross sections for multiple media, make a
    // probabilistic competition among the elements within to randomly determine
    // particular nuclei.
    std::map<G4double, G4Element *> fractions;
    G4double norm = 0;
    size_t nEl = 0;
    for( auto it = elList.cbegin(); elList.cend() != it; ++it, ++nEl ) {
        G4Element & el = **it;
        double val = fModelPtr->GetFullCrossSection( incidentE, aTrack.GetDefinition(), &el );
        if( ! std::isfinite( val ) ) continue;  // not appliable for this element
        val *= aTrack.GetMaterial()->GetVecNbOfAtomsPerVolume()[nEl];
        fractions.emplace( val, *it );
        norm += val;
    }
    assert( ! fractions.empty() );  // at least one element in list
    G4double uRandVal = G4RandFlat::shoot( norm );
    auto elIt = fractions.lower_bound( uRandVal );  // returns first >=
    G4double apTheta, apE, apPhi;
    /*const G4double prob =*/ fModelPtr->GenerateOn( incidentE, aTrack.GetDefinition(), elIt->second
                                               , &apE, &apTheta, &apPhi );
    assert( apE < incidentE );
    assert( apTheta < M_PI && apTheta > 0 );
    // Compute other kinematic parameters of the system w.r.t. generated
    // particle
    const G4double recoilE = incidentE - apE
                 , recoilTheta = sqrt( theAPrimePtr->GetPDGMass() / incidentE )
                               * (1 + theAPrimePtr->GetPDGMass() / (2*incidentE)
                                 /* + ... todo: series? */ )
                 ;
    // Initialize A' direction vector:
    G4ThreeVector aprimeDirection(0., 0., 1); {
        aprimeDirection.setMag(1.);
        aprimeDirection.setTheta( apTheta );
        aprimeDirection.setPhi( apPhi );
    }
    // Initialize new projectile particle direction vector:
    G4ThreeVector projDirection(0., 0., 1); {
        projDirection.setMag(1.);
        projDirection.setTheta( recoilTheta );
        projDirection.setPhi( - apPhi );
    }
    G4DynamicParticle * movingAPrime =  new G4DynamicParticle( theAPrimePtr
                                                             , aprimeDirection
                                                             , apE );
    aParticleChange.Initialize( aTrack );

    // Set A':
    aParticleChange.SetNumberOfSecondaries( 1 );
    aParticleChange.AddSecondary( movingAPrime );
    // Set projectile changes:
    # if 0
    // Note: if "MomentumDirection" means, actually, momentum, one
    // need to calculate momentum precisely. Possibly, using this
    // method:
    G4ThreeVector p_e = aParticleChange.CalcMomentum( E_recoil,  );
    # else
    aParticleChange.ProposeEnergy( recoilE );
    aParticleChange.ProposeMomentumDirection( projDirection );
    # endif
    return &aParticleChange;
}

}

