# include "dphmc-test-fixture.hh"

# include <gtest/gtest.h>

#if 0
class MajorantTest : public ::testing::Test
                   , public MajorantTestCommon {
protected:
    virtual void SetUp() override {
        init_caches();
    }
    virtual void TearDown() override {
        free_caches();
    }
};

// Assure the x sampling of M_{x, 1} produces resonable results
TEST_F( MajorantTest, x_chiSq_test ) {
    size_t nSamples = 1e5;
    size_t histogram[100];

    size_t nBins = sizeof(histogram)/sizeof(*histogram);
    double x;
    memset(histogram, 0, sizeof(histogram));
    for( int i = 0; i < nSamples; ++i ) {
        dphmc_aprime_ww_mj1_sample_x( &_rgs, _cachesPtr, &x );
        ASSERT_LE( x, 1 );
        ASSERT_GE( x, 0 );
        if( x < 1. ) {
            ++histogram[int(x*nBins)];
        } else {
            ++histogram[int((nBins-1))];
        }
    }
    // ...
}

#if 0
TEST_F( MajorantTest, majorant2samples_majorant1 ) {
    // ...
}
#endif
#endif
