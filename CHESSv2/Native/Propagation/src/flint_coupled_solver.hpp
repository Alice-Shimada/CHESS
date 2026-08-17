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
        int decimal_precision,
        int node_count,
        int full_dimension,
        const std::vector<int>& active,
        const std::vector<Real>& scalar_matrix,
        const std::vector<SparseOperator>& b0_operators,
        int flint_bits
    );
    ~FlintCoupledSolver();

    FlintCoupledSolver(const FlintCoupledSolver&) = delete;
    FlintCoupledSolver& operator=(const FlintCoupledSolver&) = delete;

    std::vector<Complex> solve(const std::vector<Complex>& rhs);

private:
    int decimal_precision_;
    int dimension_;
    gr_ctx_t real_context_;
    gr_ctx_t complex_context_;
    gr_mat_t lu_;
    std::vector<slong> permutation_;
    bool real_context_initialized_ = false;
    bool complex_context_initialized_ = false;
    bool lu_initialized_ = false;

    void clear() noexcept;
    std::string real_string(const Real& value) const;
    void set_component(nfloat_ptr destination, const Real& value);
    Real get_component(nfloat_srcptr source);
    void set_complex(gr_ptr destination, const Complex& value);
    Complex get_complex(gr_srcptr source);
};
