#pragma once

#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/mpfr.hpp>

#include <vector>

using Real = boost::multiprecision::mpfr_float;
using Complex = boost::multiprecision::mpc_complex;

struct SparseOperator {
    bool real_only = true;
    std::vector<int> row_ptr;
    std::vector<int> column;
    std::vector<Real> real_value;
    std::vector<Complex> complex_value;
};
