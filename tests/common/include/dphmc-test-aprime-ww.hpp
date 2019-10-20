# ifndef H_DPHMC_TEST_APRIME_WW_H
# define H_DPHMC_TEST_APRIME_WW_H

# include "numerics/aprime-cs.h"

# include <vector>

namespace dphmc {
namespace test {

/// \cite{Bjorken} has plotted an important intermediate result of their \chi
/// function at fig 10 for a series of energies. It is important check-point to
/// validate against them
/// Signature corresponds to ROOT's TF1 with parameters:
///     * p_0 -- beam energy, GeV
///     * p_1 -- material A
///     * p_2 -- material Z
/// while `x` refers to A' mass in GeV.
double chi_bjorken( double *, double * );

struct DCrossSectionPoint {double x, theta, dCS;};

/// Calculates points for A' cross-section surface of w.r.t. expression given
/// in \cite{Bjorken}, formula (A12).
std::vector<DCrossSectionPoint>
aprime_dcs_ww( const dphmc_APrimeWSParameters & apcsPs
             , size_t xNPts, size_t tNPts
             , double * xRange, double * thetaRange );

}  // namespace ::dphmc::test
}  // namespace dphmc

# endif  // H_DPHMC_TEST_APRIME_WW_H
