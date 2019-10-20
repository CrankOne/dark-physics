# ifndef H_DPHMC_INSPECTION_H
# define H_DPHMC_INSPECTION_H

/**@file inspection.h
 * @brief Pure-C utilities for inspection of ongoing numerical calculus.
 *
 * Since some numerical procedures imply invokations of complex routines, the
 * error-reporting mechanism has to be foreseen. This file is steered by
 * standard `NDEBUG` macro, conditionally enabling basic stream-based messaging
 * and error-reporting mechanics for pure C routines.
 *
 * For convenience, the default message streaming destinations are set to 
 * standard library [error] output stream(s). User may override it with their
 * custom destination to avoid undesired output or handle messages in custom
 * manner.
 */

# include <stdio.h>

# include "dphmc-config.h"

# ifdef __cplusplus
extern "C" {
# endif

/** Stream descriptor for errors of numeric calculus. Default is stderr. */
extern FILE * dphmc_stderr;
/** Stream descriptor for errors of numeric calculus. Default is stdout. */
extern FILE * dphmc_stdout;
/** On calculus error, set to descriptive code */
extern int dphmc_error;

# ifdef __cplusplus
}  // extern "C"
# endif

/* no-debug macro is defined -- disable any output except from
 * critical ones which have to be put in stderr. */
# ifdef NDEBUG
#   define DPhMC_msg( ... )  /* disabled */

# else
#   define DPhMC_msg( fmt, ... )  \
    fprintf( dphmc_stdout ? dphmc_stdout : stdout, fmt "\n", ## __VA_ARGS__ );
# endif

# define DPhMC_msgw( fmt, ... ) \
    fprintf( dphmc_stderr ? dphmc_stderr : stderr \
           , "dphmc-warning: " fmt "\n", ## __VA_ARGS__ )
# define DPhMC_msge( fmt, ... ) \
    fprintf( dphmc_stderr ? dphmc_stderr : stderr \
           , "dphmc-error: " fmt "\n", ## __VA_ARGS__ )
# define DPhMC_error( code, fmt, ... ) { \
    fprintf( dphmc_stderr ? dphmc_stderr : stderr \
           , "dphmc-fatal: " fmt "\n", ## __VA_ARGS__ ); \
    dphmc_error = code; \
    return -1; }

/*
 * Error codes
 */

# define DPHMC_FOR_EACH_ERROR_CODE( m ) \
    m( DPHMC_BAD_VALUE,         -1, "Value is not a finite number, or outside a range" ) \
    m( DPHMC_THIRD_PARTY_ERROR, -2, "Third-party error appeared" ) \
    m( DPHMC_INTEGRATION_ERROR, -3, "Integration calculus did not converge" ) \
    /* ... */

# define DPHMC_DECLARE_CODE( name, code, description ) extern const int name;
DPHMC_FOR_EACH_ERROR_CODE( DPHMC_DECLARE_CODE )
# undef DPHMC_DECLARE_CODE

# endif  /* H_DPHMC_INSPECTION_H */
