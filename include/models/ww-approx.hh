# ifndef H_DPHMC_MODELS_APRIME_WW_APPROX_H
# define H_DPHMC_MODELS_APRIME_WW_APPROX_H

# include "APrimeScatteringModel.hh"

# include <set>
# include <unordered_map>

namespace dphmc {

/**@brief A Weizsacker-Williams based approximation model for A' particle
 * production on nuclei.
 *
 * This class uses C-routines defined in eponymous .c file to calculate the
 * Weiszacker-Williams approximation according to \cite{Bjorken}. It maintains
 * a set of event generators bound with cross-section calculation workspaces,
 * spanned (tabulated) for some range of incident particles.
 *
 * Within this implementation of the model, following approximation
 * takes place: for an energy \f$E_0\f$ there is a range
 * \f$(E_i, E_{i+1}]\f$ defined so that \f$E_i < E_0 < E_{i+1}\f$
 * where we imply \f$ d \sigma / d E_0 \simeq 0 \f$. Lowest and highest values
 * correspond to cut-off values, below (beyond) which the model will not be
 * applicable. The significant drawback of any tabulated model is that these
 * bounds have to be, generally, defined emperically, w.r.t. to projectiles
 * spectra.
 *
 * Current design encourages users to customize this model and
 * generators used by inheriting this class and subsidiary AbstractGenerator
 * interface. This way one can define various generators for wide range of
 * volatile conditions. Default implementation is pretty strightforward and
 * is restricted to electrons only.
 *
 * @note the particle mass and mixing constant is related to particle
 * definition itself rather than to WW model. See APrime class definition
 * for these parameters.
 * */
class APrimeWWApproximation : public APrimeScatteringModel {
public:
    /// Common parameters of subsidiary numerical procedures
    struct Parameters {
        /// Absolute error used for Gauss-Kronrod integration of chi
        double chiGKIntAbsErr;
        /// Relative error used for Gauss-Kronrod integration of chi
        double chiGKIntRelErr;
        /// Inreasing factor to loose the relative integration requirements for
        /// Gauss-Kronrod integration on chi error
        double chiGKIntRelErrIncFt;
        /// Maximum number of subintervals used for Gauss-Kronrod integration
        size_t chiGKIntLimit;
        /// Integration workspace nodes number used for Gauss-Kronrod integration
        size_t chiGKIntNNodes;
    };

    /// General Weiszacker-Williams event generator interface used by model.
    struct AbstractGenerator {
        /// Shall return full (integrated) cross section value for conditions
        /// this generator is dedicated to
        virtual G4double GetFullCrossSection() const = 0;
        /// Shall returns probability of the event if possible after writing
        /// the kinematics to variables given by ptrs.
        virtual G4double ShootKinematics( const G4double projectileEnergy
                                        , G4double * energy
                                        , G4double * theta
                                        , G4double * phi ) const = 0;
    };

    /// A key-like structure uniquely identifying the particular generator obj
    struct GeneratorKey {
        /// Defined by charge number of target nucleus.
        G4double materialZ;
        /// Defined by projectile energy; refers to APrimeWWApproximation::fRanges index.
        std::set<G4double>::difference_type energyRange;
        /// Projectile particle type.
        G4int pType;
        /// is-equal needed for std::unordered_map.
        bool operator==(const GeneratorKey &) const;
    };

    /// Hash function for custom hash key.
    struct HashGeneratorKey {
        ::std::size_t operator()(const GeneratorKey &) const noexcept;
    };

    /// Type alias for generators registry container type.
    typedef std::unordered_map< GeneratorKey
                              , AbstractGenerator*
                              , HashGeneratorKey> GeneratorsIndex;
private:
    /// Lowest charge number of nucleus under consideration. In Geant4's
    /// G4Element it is a float, so for compat it is float here as well.
    G4double fLowestZ;
    /// Sorted set of energy ranges. Immutable, as its index is used to choose CS.
    std::set<G4double> fRanges;
    /// Parameters instance
    Parameters fModelParameters;
    /// Associative container indexing all the generators created during the
    /// session by their energy range, material parameters, etc.
    GeneratorsIndex fGenerators;
protected:
    /// Tries to use existing generator instances to retrieve the full cross
    /// section value. If it does not exist, creates and emplaces new one with
    /// NewGenerator() method.
    virtual GeneratorsIndex::iterator RetrieveGeneratorFor( G4double incidentParticleEnergy
                                , const G4ParticleDefinition * pDef
                                , const G4Element * g4Element );
public:
    /// Model constructor. Takes most important parameterisation.
    APrimeWWApproximation( const G4double lowestZ
                         , G4double * incidentEs );

    /// Returns estimation yielded by analytic estimation
    virtual G4double GetFullCrossSection( G4double incidentParticleEnergy
                                        , const G4ParticleDefinition * pDef
                                        , const G4Element * g4Element ) override;
    /// Generates energy and angles of produced A'
    virtual G4double GenerateOn( G4double incidentParticleEnergy
                               , const G4ParticleDefinition * pDef
                               , const G4Element * g4Element
                               , G4double * energy
                               , G4double * theta
                               , G4double * phi ) override;
    /// Returns `true` if particle is an electron.
    virtual G4bool IsApplicableToParticle( const G4ParticleDefinition & ) override;

    /// Instantiates new generator instance for given conditions.
    virtual AbstractGenerator * NewGenerator( G4double projELow, G4double projEUp
                                            , const G4ParticleDefinition * pDef
                                            , const G4Element * g4Element );

    /// Returns const instance of current index container. This getter is for
    /// inspection purposes only.
    const GeneratorsIndex & GetGeneratorsIndex() const { return fGenerators; }

    /// Returns energy tabulation bounds. This getter is for inspection
    /// purposes only.
    const std::set<G4double> & GetEnergyTabulation() const { return fRanges; }
};

}  // namespace dphmc

# endif  // H_DPHMC_MODELS_APRIME_WW_APPROX_H

