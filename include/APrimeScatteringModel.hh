# ifndef H_DPHMC_APRIME_MIXING_H
# define H_DPHMC_APRIME_MIXING_H

# include "dphmc-config.h"

# include <G4Types.hh>

class G4ParticleDefinition;
class G4Element;

namespace dphmc {

class APrimeGenerator;

///@brief  Interface for A' model used by APrimeEMMixingProcess.
///@details  Defines auxiliary class for picking up concrete configured dark photon
///physics model.
class APrimeScatteringModel {
public:
    /// Based on material parameters and incident particle energy and type,
    /// returns a full (integrated) cross section value.
    virtual G4double GetFullCrossSection( G4double incidentParticleEnergy
                                      , const G4ParticleDefinition * pDef
                                      , const G4Element * g4Element ) = 0;
    ///@brief Generates A' parameters.
    ///@details Based on material parameters and incident particle energy and
    /// type, generate the angles and energy fraction of produced particle.
    /// For the purpose of weighted tests, must return a statistical weight
    /// when possible. For the models that does not supply this information,
    /// 0 has to be returned.
    virtual G4double GenerateOn( G4double incidentParticleEnergy
                               , const G4ParticleDefinition * pDef
                               , const G4Element * g4Element
                               , G4double * energy
                               , G4double * theta
                               , G4double * phi ) = 0;
    /// Shall return `true` for incident particles that may produce A' via
    /// eloctromagnetic mixing while being scattered on nuclei within this
    /// model.
    virtual G4bool IsApplicableToParticle( const G4ParticleDefinition & ) = 0;
};

}  // namespace dphmc

# endif  // H_DPHMC_APRIME_MIXING_H
