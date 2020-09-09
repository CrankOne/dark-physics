# ifndef H_DPHMC_TEST_APRIME_WW_H
# define H_DPHMC_TEST_APRIME_WW_H

# include "models/ww-approx.hh"

# include <vector>

namespace dphmc {
namespace test {

/// Useful class to generate exponentially-spaced grids
class ExpValues {
private:
    size_t _nPoints;
    double _low, _up, _current, _step;
public:
    ExpValues( size_t n_, double low_, double up_ );
    double value() const;
    bool good() const;
    double operator++();
};

/// \cite{Bjorken} has plotted an important intermediate result of their \chi
/// function at fig 10 for a series of energies. It is important check-point to
/// validate against them
/// Signature corresponds to ROOT's TF1 with parameters:
///     * p_0 -- beam energy, GeV
///     * p_1 -- material A
///     * p_2 -- material Z
/// while `x` refers to A' mass in GeV.
double chi_bjorken( double *, double * );

struct DCrossSectionPoint {
    double x, theta, dCS
         , saDDTheta
         , numDDTheta, numDDThetaErr
         ;
};

/// Calculates points for A' cross-section surface of w.r.t. expression given
/// in \cite{Bjorken}, formula (A12).
std::vector<DCrossSectionPoint>
aprime_dcs_ww( const aprime::PhysParameters & apcsPs
             , const IterativeQAGSParameters & qagsp
             , size_t xNPts, size_t tNPts
             , double * xRange, double * thetaRange );

/// Calculates diff. CS maximum points
std::vector< std::tuple<double, double, double> >
aprime_dcs_max( const aprime::PhysParameters & apcsPs
              , const IterativeQAGSParameters & qagsp
              , size_t xNPts );

}  // namespace ::dphmc::test
}  // namespace dphmc

# endif  // H_DPHMC_TEST_APRIME_WW_H
