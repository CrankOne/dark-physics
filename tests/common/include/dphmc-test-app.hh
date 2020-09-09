# ifndef H_DPHMC_TEST_FIXTURE_H
# define H_DPHMC_TEST_FIXTURE_H

# include "models/ww-approx.h"

# include <map>
# include <string>

# include "models/ww-approx.hh"

namespace dphmc {
namespace test {

/**\brief A parameters dictionary for testing applications
 *
 * This primitive parameters registry implements a (very basic) storage for
 * user-defined parameters to shorten application code. It indexes parameters
 * name versus the (value,description) pairs for setting and retrieval.
 * */
class ParametersSet : protected std::map<std::string, std::pair<double, std::string>> {
public:
    /// Defines a parameter for further modification; may raise
    /// `std::runtime_error` on duplicating name.
    void define( const std::string & name
               , const std::string & description
               , double value );

    /// Prints list of defined parameters in the
    /// form "<name=value> -- description", one per line
    void print_references( std::ostream & ) const;

    /// Returns the parameter value reference (available for modification),
    /// raises `std::runtime_error` if no parameter defined with such name
    double & operator[](const std::string &);

    /// Returns the parameter value (const version),
    /// raises `std::runtime_error` if no parameter defined with such name
    double operator[](const std::string &) const;
};

/// Adds ma, E0
void define_aprime_particle_parameters( ParametersSet & );
/// Adds A,Z
void define_media_parameters( ParametersSet & );
/// Adds epsilon, chiEpsAbs, chiEpsRel, chiEpsRelInc, chiLimit, chiNNodes
void define_cross_section_parameters( ParametersSet & ps );
/// Adds seed
void define_rng_parameters( ParametersSet & ps );

/// Reads the `APrimePhysParameters` values from parameters set (must have
/// Z, A, ma, E0, epsilon)
void configure_phys_parameters( const ParametersSet & ps
                              , aprime::PhysParameters & );

/// Reads the `IterativeQAGSParameters` values from parameters set (must have
/// chiEpsAbs, chiEpsRel, chiEpsRelInc, chiLimit, chiNNodes)
void configure_chiint_parameters( const ParametersSet & ps
                                , IterativeQAGSParameters & );

/// Initializes random number generator
void configure_rng( const ParametersSet & ps
                  , URandomState & rgs );

}
}

# endif  // H_DPHMC_TEST_FIXTURE_H

