#ifndef H_DPHMC_TESTING_HISTOGRAM_H
#define H_DPHMC_TESTING_HISTOGRAM_H

#include <cstdlib>

namespace dphmc {
namespace test {

struct HistogramAxis {
    size_t nBins;
    double range[2];
    HistogramAxis( size_t nBins_, double vMin, double vMax );
    virtual size_t value_to_nbin( double v ) const;
};

template<int Dt> class Histogram;

template<>
class Histogram<1> {
protected:
    HistogramAxis _axis;
    size_t _underflow, _overflow, _sum;
    size_t * _bins;
public:
    Histogram( size_t nBins, double vMin, double vMax );
    ~Histogram();
    void fill(double v);
    const HistogramAxis & axis() const { return _axis; }
    size_t underflow() const { return _underflow; }
    size_t overflow() const { return _overflow; }
    const size_t * bins() { return _bins; }
    size_t sum() const { return _sum; }
};

}
}

#endif  // H_DPHMC_TESTING_HISTOGRAM_H
