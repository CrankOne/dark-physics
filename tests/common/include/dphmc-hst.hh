#ifndef H_DPHMC_TESTING_HISTOGRAM_H
#define H_DPHMC_TESTING_HISTOGRAM_H

#include <cstdlib>

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

/// Histogram axis representation -- uniform
struct HistogramAxis {
    size_t nBins;
    double range[2];
    HistogramAxis( size_t nBins_, double vMin, double vMax );
    /// Returns bin number for value
    virtual ssize_t value_to_nbin( double v ) const;
    /// Returns lower bound of n-th bin
    virtual double lower_bound( size_t n ) const;
};

/// Logarithmic scale log axis
struct HistogramLogAxis : public HistogramAxis {
    HistogramLogAxis(size_t n, double min, double max );
    virtual ssize_t value_to_nbin( double v ) const;
    virtual double lower_bound( size_t n ) const;
};

/// A histogram class template
template<int Dt> class Histogram;

/// 1D histogram specialization
template<> class Histogram<1> {
protected:
    HistogramAxis * _axis;
    size_t _underflow, _overflow, _sum;
    size_t * _bins;
public:
    Histogram( size_t nBins, double vMin, double vMax
             , bool isLog=false );
    ~Histogram();
    void fill(double v);
    const HistogramAxis & axis() const { return *_axis; }
    size_t underflow() const { return _underflow; }
    size_t overflow() const { return _overflow; }
    const size_t * bins() { return _bins; }
    size_t sum() const { return _sum; }
};


/// 2D histogram specialization
template<> class Histogram<2> {
protected:
    HistogramAxis * _axes[2];
    size_t _underflow, _overflow, _sum;
    size_t * _bins;
public:
    Histogram( size_t nBins1, double v1Min, double v1Max, bool is1Log
             , size_t nBins2, double v2Min, double v2Max, bool is2Log );
    ~Histogram();
    void fill(double v1, double v2);
    const HistogramAxis & axis1() const { return *_axes[0]; }
    const HistogramAxis & axis2() const { return *_axes[1]; }
    const size_t underflow() const { return _underflow; }
    const size_t overflow() const { return _overflow; }
    const size_t * bins() { return _bins; }
    size_t sum() const { return _sum; }

    size_t & at( size_t n1, size_t n2 );
    size_t at( size_t n1, size_t n2 ) const;

    size_t & operator()( size_t n1, size_t n2 ) { return at(n1, n2); }
    size_t operator()( size_t n1, size_t n2 ) const { return at(n1, n2); }
};

}
}

#endif  // H_DPHMC_TESTING_HISTOGRAM_H
