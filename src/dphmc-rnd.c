#include "dphmc-rnd.h"
#include "dphmc-inspection.h"

#include <assert.h>
#include <string.h>

#if 0  // XXX
/*                           ___________________
 * ________________________/ C stdlib Generator \_________________________________
 */

struct CStdlibRandGen {
    char * stateBuf;
    struct random_data * buf;
};

int
dphmc_rnd_gen_cstdlib_init( struct dphmc_URandomState * stateWPtr
                          , size_t stateSize
                          , int seed ) {
    assert(stateWPtr);
    assert(stateSize >= 8);

    struct CStdlibRandGen * r = malloc(sizeof(struct CStdlibRandGen));
    r->stateBuf = malloc(stateSize);
    r->buf = malloc(sizeof(struct random_data));

    r->buf->state = NULL;
    int rc = initstate_r( seed, r->stateBuf, stateSize, r->buf );
    if( 0 != rc ) {
        int errno_ = errno;
        DPhMC_error( DPHMC_THIRD_PARTY_ERROR
                   , "initstate_r() returned status %d: %s"
                   , errno_
                   , strerror(errno_) );
    }
    assert(rc);  // TODO: real check

    stateWPtr->statePtr = r;
    stateWPtr->urandomi = dphmc_rnd_gen_cstdlib_urandomi;
    stateWPtr->urandomf = dphmc_rnd_gen_cstdlib_urandomf;

    return 0;
}

int
dphmc_rnd_gen_cstdlib_free( struct dphmc_URandomState * stateWPtr ) {
    assert(stateWPtr);
    assert(stateWPtr->statePtr);
    struct CStdlibRandGen * r = (struct CStdlibRandGen *) stateWPtr->statePtr;

    assert(r->buf);
    assert(r->stateBuf);

    free(r->buf);
    free(r->stateBuf);
    free(r);

    stateWPtr->statePtr = NULL;
    stateWPtr->urandom = NULL;

    return 0;
}

int
dphmc_rnd_gen_cstdlib_urandom( void * state, int32_t * ) {
    assert(state);

    struct CStdlibRandGen * r = (struct CStdlibRandGen *) state;
    int32_t ui;
    int rc = random_r( r->buf, &ui );
    if( 0 != rc )
}
#endif

/*                             ______________
 * __________________________/ GSL Generator \_________________________________
 */

/** Runs `gsl_rng_alloc()` and `gsl_rng_set()` functions to initialize the
 * generator state instance and set the seed.
 * See also `dphmc_rnd_gen_gsl_free()` to delete the instance afterwards. */
int
dphmc_rnd_gen_gsl_init( struct dphmc_URandomState * stateWPtr
                      , const gsl_rng_type * gslRngTypePtr
                      , int seed ) {
    gsl_rng * r = gsl_rng_alloc( gslRngTypePtr );
    if( NULL == r ) {
        DPhMC_error( DPHMC_THIRD_PARTY_ERROR
                   , "Unable to allocate GSL random number generator." );
    }
    gsl_rng_set( r, seed );
    stateWPtr->statePtr = r;
    stateWPtr->urandom_f = dphmc_rnd_gen_gsl_urandom_f;
    return 0;
}

/** Runs `gsl_rng_free()` to delete the GSL generator state instance previously
 * allocated with `dphmc_rnd_gen_gsl_init()`. */
int
dphmc_rnd_gen_gsl_free( struct dphmc_URandomState * stateWPtr ) {
    assert( stateWPtr->statePtr );
    gsl_rng_free( (gsl_rng *) stateWPtr->statePtr );
    stateWPtr->statePtr = NULL;
    stateWPtr->urandom_f = NULL;
    return 0;
}

/** Note that contrary to native GSL functions like `gsl_rng_uniform()` this
 * function returns random value with closed bounds, designed for case when
 * generated value *may* be 0 or 1 as well as all the values beyond. */
int
dphmc_rnd_gen_gsl_urandom_f( void * statePtr, double * vPtr ) {
    assert(statePtr);
    gsl_rng * r = (gsl_rng *) statePtr;
    *vPtr = ((double) gsl_rng_get(r))/gsl_rng_max(r);
    return 0;
}

