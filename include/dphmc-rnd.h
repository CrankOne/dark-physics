# ifndef H_DPHMC_RANDOM_GEN_H
# define H_DPHMC_RANDOM_GEN_H

#include <stdlib.h>
#include <gsl/gsl_rng.h>

# ifdef __cplusplus
extern "C" {
# endif

/**\brief A reference to choosen random number generator state
 *
 * This structure keeps a primitive API for uniform random number sampling. As
 * library itself is not restricted to Geant4 usage we do provide this
 * shimmering layer to handle some other random number generators.
 *
 * The uniform random number sampling function is kept within this structure as
 * a callback to minimize a performance loss due to indirection and to give
 * the user code an opportunity to utilize their own random number generators.
 */
struct dphmc_URandomState {
    /** Use this for other */
    void * statePtr;

    /** Random number generating callback function: uniform for 64-bit float */
    int (*urandom_f)(void *, double * );
};

/*                             ______________
 * __________________________/ GSL Generator \_________________________________
 */

/** Initializes GSL random number generator */
int dphmc_rnd_gen_gsl_init( struct dphmc_URandomState * statePtr
                          , const gsl_rng_type * gslRngTypePtr
                          , int seed );

/** Frees GSL random number generator */
int dphmc_rnd_gen_gsl_free( struct dphmc_URandomState * statePtr );

/** Generates a uniformly distributed random number \f$[0:1]\f$ */
int dphmc_rnd_gen_gsl_urandom_f( void *, double * );

# ifdef __cplusplus
}
# endif

# endif  /* H_DPHMC_RANDOM_GEN_H */
