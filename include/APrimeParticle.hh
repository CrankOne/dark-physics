# ifndef H_A_PRIME_PARTICLE_H
# define H_A_PRIME_PARTICLE_H

# include "dphmc-config.h"
# include <G4ParticleDefinition.hh>

namespace dphmc {

/**@class APrime
 * @brief Defines an A' particle («dark photon») in the Geant4 API terms.
 *
 * Definition of a custom particle in Geant4 is required to represent
 * all connected process. A' does not participate in any interactions
 * except gravitational and reveals itself by absence of energy in
 * electromagnetic showers when producted or by decay in quite
 * sophisticated schemas involving post-production of other exotic
 * particles.
 *
 * All the particles in Geant4 are implemented as singletons.
 *
 * @todo In the Definition() there may be some of unknown/wrong information
 * about this particle that must be corrected further (isospin numbers,
 * mass, leptonic number, etc).
 * @todo Decay physics is incomplete and not checked yet.
 */
class APrime : public G4ParticleDefinition {
public:
    /// Parameters of A' particle that might be changed prior to instantiation.
    static struct Parameters {
        double mass  ///< (hypothetical) mass of A' particle. Default is 16.7 MeV
             , mixingFactor  ///< Mixing constant, \f$\epsilon\f$ regulating the interaction probability.
             ;
        bool enableDecays;  ///< Whether to enable decay physics; default is false.
        int codePDG;  ///< PDG-like code, has to be a number in 81-100; default is 88.
    } parameters;
private:
    static APrime * theInstance;
    APrime();
    ~APrime();
public:
    /// Geant4's std way of particle definition declaration.
    static APrime * Definition();
    /// Geant4's std alias for Definition().
    static APrime * APrimeDefinition();
    // Note: another alias for Defention() method which we decide
    // not to implement due to the ctr name collision.
    //static APrime * APrime();
};

}

# endif  // H_A_PRIME_PARTICLE_H

