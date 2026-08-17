#pragma once

#include "backend_types.hpp"

#include <flint/gr.h>
#include <flint/gr_mat.h>
#include <flint/nfloat.h>

#include <string>
#include <vector>

class FlintCoupledSolver {
public:
    FlintCoupledSolver(
        int requestedPrecision,
        int nodeCount,
        int fullDimension,
        const std::vector<int>& active,
        const std::vector<Real>& scalarMatrix,
        const std::vector<SparseOperator>& b0Operators,
        int flintBits
    );
    ~FlintCoupledSolver();

    FlintCoupledSolver(const FlintCoupledSolver&) = delete;
    FlintCoupledSolver& operator=(const FlintCoupledSolver&) = delete;

    std::vector<Complex> Solve(const std::vector<Complex>& rhs);

private:
    int decimalPrecision;
    int dimension;
    gr_ctx_t realContext;
    gr_ctx_t complexContext;
    gr_mat_t lu;
    std::vector<slong> permutation;
    bool realContextInitialized = false;
    bool complexContextInitialized = false;
    bool luInitialized = false;

    void Clear() noexcept;
    std::string RealString(const Real& value) const;
    void SetComponent(nfloat_ptr destination, const Real& value);
    Real GetComponent(nfloat_srcptr source);
    void SetComplex(gr_ptr destination, const Complex& value);
    Complex GetComplex(gr_srcptr source);
};
