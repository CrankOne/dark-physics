#include "dphmc-rnd.h"

#include <gtest/gtest.h>

TEST( RandomGenAPI, gsl_rand_gen_works ) {
    int rc;  // return code (C diagnostic variable)
    double u;  // uniformely-distributed random variable
    // test sequences descriptions
    struct SeqDesc {
        int seed;
        size_t n;
    } sqs[] = { {    0,  500  }
              , {  123,  1000 }
              , { 8891,  100  }
              , { 9933,  100  }
              };
    // saved sequences
    std::map<int, std::vector<double>> saved;
    // fill the sequences
    for( int n = 0; n < sizeof(sqs)/sizeof(*sqs); ++n ) {
        struct dphmc_URandomState rgs;
        std::vector<double> seq;
        // create the generator and fill the sequence
        dphmc_rnd_gen_gsl_init( &rgs
                              , gsl_rng_ranlux
                              , sqs[n].seed
                              ); {
            for( size_t i = 0; i < sqs[n].n; ++i ) {
                rc = rgs.urandom_f(rgs.statePtr, &u);
                ASSERT_EQ(rc, 0);
                seq.push_back( u );
            }
        } dphmc_rnd_gen_gsl_free( &rgs );
        saved[sqs[n].seed] = seq;
    }
    // check the sequences for match
    for( int n = 0; n < sizeof(sqs)/sizeof(*sqs); ++n ) {
        struct dphmc_URandomState rgs;
        std::vector<double> seq;
        // create the generator and fill the sequence
        dphmc_rnd_gen_gsl_init( &rgs, gsl_rng_ranlux, sqs[n].seed ); {
            for( size_t i = 0; i < sqs[n].n; ++i ) {
                rc = rgs.urandom_f(rgs.statePtr, &u);
                ASSERT_EQ(rc, 0);
                seq.push_back( u );
            }
        } dphmc_rnd_gen_gsl_free( &rgs );
        // compare
        // TODO: instead of element-wise comparison, the checksum would be
        // preferable here
        const std::vector<double> & prev = saved[sqs[n].seed];
        ASSERT_EQ( prev.size(), seq.size() );
        for( int i = 0; i < prev.size(); ++i ) {
            ASSERT_NEAR( prev[i], seq[i], 1e-12 );
        }
    }
}

