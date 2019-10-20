# include "dphmc-inspection.h"

FILE * dphmc_stderr = NULL;
FILE * dphmc_stdout = NULL;
int dphmc_error = 0;

# define DPHMC_DECLARE_CODE( name, code, description ) const int name = code;
DPHMC_FOR_EACH_ERROR_CODE( DPHMC_DECLARE_CODE )
# undef DPHMC_DECLARE_CODE

