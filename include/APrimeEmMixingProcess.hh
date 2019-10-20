# ifndef H_A_PRIME_PRODUCTION_PROCESS_H
# define H_A_PRIME_PRODUCTION_PROCESS_H

# include "dphmc-config.h"

# include <Geant4/G4VDiscreteProcess.hh>
# include "APrimeParticle.hh"

namespace dphmc {

class APrimeScatteringModel;

/**\class AMixingProcess
 * \brief Implements A' production in the Geant4 API terms.
 *
 * Following the distirbuted logic of Geant4 API, this class is only
 * responsible for proper routing probabilistic invocations of the model. It
 * is not responsible for the cross-section calculus, neither for the random
 * number generation.
 *
 * The Geant4 discrete process implements a following interface:
 *  1. Based on material compound of curent Geant4 step (which lies completely
 *  within a homogeneous volume filled with uniform media) the discrete process
 *  interface must supply tracking algorithm with mean free path. It shall be
 *  done in GetMeanFreePath() below and for our case it implis calulation of
 *  full cross section values on all the nuclei types of media.
 *  2. Taking into account this value, the tracking algorithm "throws a dice"
 *  for all the processes taken into account for this step. Winning process is
 *  then has to propose for the tracking algorithm the changes for the current
 *  track on current step with PostStepDoIt() method.
 *
 * This class implements invocation of methods of
 * APrimeScatteringModel interface in order to obtain corresponding values.
 * */
class APrimeEmMixingProcess : public G4VDiscreteProcess {
public:
    /// Ptr to common definition of A' particle
    APrime * theAPrimePtr;
private:
    /// no explicit ctr
    APrimeEmMixingProcess() = delete;
    /// no copy ctr
    APrimeEmMixingProcess( const APrimeEmMixingProcess & ) = delete;
public:
    /// Only allowed ctr -- instantiates process bound to particular model.
    APrimeEmMixingProcess( APrimeScatteringModel * modelPtr );

    ///@brief Implements final state parameters when process won.
    virtual G4VParticleChange* PostStepDoIt( const G4Track & ,
			                                 const G4Step & ) override;
    ///@brief Returns event probability, \f$mm^2\f$.
    virtual G4double GetMeanFreePath( const G4Track & aTrack,
                                      G4double previousStepSize,
                                      G4ForceCondition * condition ) override;
    ///@brief Returns `true` for incident particles under consideration.
    virtual G4bool IsApplicable(const G4ParticleDefinition &) override;
private:
    /// Pointer to the current model used.
    APrimeScatteringModel * fModelPtr;
};  // class AMixingProcess

}  // namespace dphmc

# endif  // H_A_PRIME_PRODUCTION_PROCESS_H

