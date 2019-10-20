# include "dphmc-test-aprime-ww.hpp"

# include <gtest/gtest.h>

/// Ref. values for chi on W, beam E -- 1GeV, depending on mA, for check
static const double tungst200MeV_chi[] = {
    4.922800e-03, 7.792763e+00,         7.656892e-03, 6.973366e+00,
    1.331619e-02, 5.808514e+00,         2.195422e-02, 4.701498e+00,
    3.654872e-02, 3.429819e+00,         5.059613e-02, 2.529155e+00,
    7.038344e-02, 1.575898e+00,         9.304399e-02, 8.569128e-01,
    1.110798e-01, 5.041809e-01,         1.235988e-01, 3.423679e-01,
    //1.365309e-01, 2.259167e-01,   // violates cut-off
    //1.515500e-01, 1.377724e-01,   // violates cut-off
    //1.610303e-01, 1.001477e-01,   // violates cut-off
    0.
};

/// Ref. values for chi on W, beam E -- 1GeV, depending on mA, for check
static const double tungst1GeV_chi[] = {
    4.887086e-03, 9.907672e+00,         7.401088e-03, 9.907672e+00,
    1.057410e-02, 9.422838e+00,         1.578210e-02, 8.677289e+00,
    2.502866e-02, 7.654369e+00,         3.602034e-02, 6.703812e+00,
    5.376120e-02, 5.504511e+00,         7.317043e-02, 4.487482e+00,
    1.003147e-01, 3.356856e+00,         1.310122e-01, 2.396786e+00,
    1.637895e-01, 1.645148e+00,         1.979262e-01, 1.101250e+00,
    2.261919e-01, 7.668097e-01,         2.661336e-01, 4.463424e-01,
    2.975684e-01, 2.966444e-01,         3.279061e-01, 1.971534e-01,
    3.675282e-01, 1.185203e-01,         3.811551e-01, 1.005072e-01,
    0.
};

/// Ref. values for chi on W, beam E -- 1GeV, depending on mA, for check
static const double tungst6GeV_chi[] = {
    4.875239e-03, 1.023244e+01,         6.881349e-03, 1.087523e+01,
    1.004866e-02, 1.147585e+01,         1.547862e-02, 1.172528e+01,
    2.216839e-02, 1.135313e+01,         3.268775e-02, 1.056786e+01,
    5.337118e-02, 9.123762e+00,         8.301319e-02, 7.626986e+00,
    1.183153e-01, 6.307579e+00,         1.711037e-01, 4.855612e+00,
    2.149521e-01, 3.972688e+00,         2.760010e-01, 3.014672e+00,
    3.543885e-01, 2.106678e+00,         4.303334e-01, 1.466897e+00,
    5.038711e-01, 1.017759e+00,         5.986298e-01, 6.162344e-01,
    6.891190e-01, 3.771520e-01,         7.723878e-01, 2.462081e-01,
    8.167309e-01, 1.943474e-01,         8.999983e-01, 1.282431e-01,
    9.493606e-01, 1.005072e-01,         0.
};

static const double alu6GeV_chi[] = {
    4.995013e-03, 1.224054e+01,         6.491961e-03, 1.273272e+01,
    1.065138e-02, 1.324469e+01,         1.907138e-02, 1.264178e+01,
    3.398210e-02, 1.147585e+01,         6.173776e-02, 9.943243e+00,
    9.463827e-02, 8.646247e+00,         1.457776e-01, 7.124936e+00,
    2.368669e-01, 5.216418e+00,         3.630951e-01, 3.516950e+00,
    4.765143e-01, 2.475352e+00,         6.074111e-01, 1.627556e+00,
    7.686476e-01, 9.541776e-01,         9.447634e-01, 5.339359e-01,
    1.144444e+00, 2.841575e-01,         1.320640e+00, 1.720522e-01,
    1.462363e+00, 1.193728e-01,         1.531378e+00, 1.001477e-01,
    0.
};

// Tests effective photon flux from \cite{Bjorken} eq (A18) versus plot given
// on fig 10. (Schuster, Toro and others).
// See utils/apime-math for plotting utility.
TEST( Numerics, bjorkenChi ) {
    struct {
        std::string label;
        double values[3];
        const double * checkpoints;
    } ps[] = {
        { "W-74-183-.2",    {  0.2, 183.84, 74 },   tungst200MeV_chi },
        { "W-74-183-1",     {    1, 183.84, 74 },   tungst1GeV_chi   },
        { "W-74-183-6",     {    6, 183.84, 74 },   tungst6GeV_chi   },
        { "Al-13-26-6",     {    6,  26.98, 13 },   alu6GeV_chi      },
        { "", {-1} }
    };
    for( auto p = ps; p->values[0] > 0; ++p ) {
        for( const double * am = p->checkpoints
           ; *am && p->values[0] ; am += 2 ) {
            // We have to const_cast<> here since the signature of
            // dphmc_test_chi_bjorken() has to be kept valid for ROOT's TF1
            // that does not respect const validity, despite the value by ptr
            // is not changed
            const double tstVal = dphmc::test::chi_bjorken(
                                  const_cast<double *>(am), p->values )
                    , refVal = *(am + 1);
            // The original plot from \cite{Bjorken} is given in logarithmic
            // scale and reference points were then manually digitized with
            // poor precision. It seems to still be in reasonable limits but
            // not very precise.
            ASSERT_NEAR( tstVal, refVal, 1e-1 );
        }
    }
}

