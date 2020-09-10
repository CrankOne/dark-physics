#include "dphmc-hst.hh"

#include <cassert>
#include <cstring>
#include <cmath>

namespace dphmc {
namespace test {

ExpValues::ExpValues( size_t n_, double low_, double up_ )
        : _nPoints(n_), _low(low_), _up(up_), _current(_low)
        , _step( pow(_up/_low, 1./(_nPoints - 1)) ) {}

double
ExpValues::value() const {
    return _current;
}

bool
ExpValues::good() const {
    return _current <= _up + _step*1e-9;
}

double
ExpValues::operator++() {
    _current *= _step;
    return _current;
}


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


Histogram<2>::Histogram( size_t nBins1, double v1Min, double v1Max
                       , size_t nBins2, double v2Min, double v2Max )
        : _axes{ HistogramAxis(nBins1, v1Min, v1Max)
               , HistogramAxis(nBins2, v2Min, v2Max) }
        , _underflow(0), _overflow(0), _sum(0) {
    _bins = new size_t [nBins1*nBins2];
    memset( _bins, 0, sizeof(size_t)*nBins1*nBins2 );
}

Histogram<2>::~Histogram() {
    delete [] _bins;
}

void
Histogram<2>::fill(double v1, double v2) {
    size_t n1 = _axes[0].value_to_nbin(v1)
         , n2 = _axes[1].value_to_nbin(v2)
         ;
    if(n1 < 0 || n2 < 0) {
        ++_underflow;
        return;
    }
    if( n1 >= _axes[0].nBins
     || n2 >= _axes[1].nBins) {
        ++_overflow;
        return;
    }
    ++_sum;
    ++_bins[n1*_axes[1].nBins + n2];
}

size_t &
Histogram<2>::at( size_t n1, size_t n2 ) {
    assert( n1 < _axes[0].nBins );
    assert( n2 < _axes[1].nBins );
    return _bins[n1*_axes[1].nBins + n2];
}

size_t
Histogram<2>::at( size_t n1, size_t n2 ) const {
    assert( n1 < _axes[0].nBins );
    assert( n2 < _axes[1].nBins );
    return _bins[n1*_axes[1].nBins + n2];
}

}
}

