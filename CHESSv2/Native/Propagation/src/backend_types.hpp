#pragma once

#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>

#include <vector>

using Real = boost::multiprecision::mpfr_float;
using Complex = boost::multiprecision::mpc_complex;

// Row-compressed matrix sampled at one Lobatto node.  Real operators avoid
// allocating complex coefficients; complexValue is populated only when a
// genuinely complex entry is present.
struct SparseOperator {
    bool realOnly = true;
    std::vector<int> rowPtr;
    std::vector<int> column;
    std::vector<Real> realValue;
    std::vector<Complex> complexValue;
};
