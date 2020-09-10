#include "dphmc-hst.hh"

#include <cassert>
#include <cstring>

namespace dphmc {
namespace test {

HistogramAxis::HistogramAxis( size_t nBins_, double vMin, double vMax ) 
            : nBins(nBins_)
            , range{vMin, vMax} {
    assert(nBins_);
    assert(vMin < vMax);
}

size_t
HistogramAxis::value_to_nbin( double v ) const {
    if( v == range[1] ) { return nBins-1; }
    return (v - range[0])*nBins/(range[1] - range[0]);
}


Histogram<1>::Histogram( size_t nBins, double vMin, double vMax )
        : _axis(nBins, vMin, vMax)
        , _underflow(0), _overflow(0), _sum(0) {
    _bins = new size_t [nBins];
    memset( _bins, 0, sizeof(size_t)*nBins );
}

Histogram<1>::~Histogram() {
    delete [] _bins;
}

void
Histogram<1>::fill(double v) {
    size_t n = _axis.value_to_nbin(v);
    if(n < 0) {
        ++_underflow;
        return;
    }
    if(n >= _axis.nBins) {
        ++_overflow;
        return;
    }
    ++_sum;
    ++_bins[n];
}

}
}

