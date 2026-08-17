#include "flint_coupled_solver.hpp"

#include <flint/flint.h>

#include <ios>
#include <stdexcept>

FlintCoupledSolver::FlintCoupledSolver(
    int decimal_precision,
    int node_count,
    int full_dimension,
    const std::vector<int>& active,
    const std::vector<Real>& scalar_matrix,
    const std::vector<SparseOperator>& b0_operators,
    int flint_bits
) : decimal_precision_(decimal_precision),
    dimension_(node_count * static_cast<int>(active.size())) {
    try {
        if (flint_bits <= 0 || dimension_ <= 0 ||
            scalar_matrix.size() !=
                static_cast<std::size_t>(node_count) * node_count ||
            b0_operators.size() != static_cast<std::size_t>(node_count - 1)) {
            throw std::runtime_error("invalid coupled-solver dimensions");
        }
        if (nfloat_ctx_init(real_context_, flint_bits, 0) != GR_SUCCESS) {
            throw std::runtime_error("failed to initialise real nfloat context");
        }
        real_context_initialized_ = true;
        if (nfloat_complex_ctx_init(complex_context_, flint_bits, 0) !=
            GR_SUCCESS) {
            throw std::runtime_error(
                "failed to initialise complex nfloat context"
            );
        }
        complex_context_initialized_ = true;
        gr_mat_init(lu_, dimension_, dimension_, complex_context_);
        lu_initialized_ = true;
        if (gr_mat_zero(lu_, complex_context_) != GR_SUCCESS) {
            throw std::runtime_error("failed to zero coupled matrix");
        }

        const int active_count = static_cast<int>(active.size());
        for (int row_node = 0; row_node < node_count; ++row_node) {
            for (int column_node = 0; column_node < node_count; ++column_node) {
                const Real& derivative = scalar_matrix[
                    static_cast<std::size_t>(row_node) * node_count + column_node
                ];
                if (derivative == 0) continue;
                for (int index = 0; index < active_count; ++index) {
                    set_complex(
                        gr_mat_entry_ptr(
                            lu_, row_node * active_count + index,
                            column_node * active_count + index, complex_context_
                        ),
                        Complex(derivative, 0)
                    );
                }
            }
        }

        std::vector<int> active_position(full_dimension, -1);
        for (int index = 0; index < active_count; ++index) {
            active_position[active[index]] = index;
        }
        gr_ptr temporary = gr_heap_init(complex_context_);
        try {
            for (int node = 1; node < node_count; ++node) {
                const SparseOperator& op = b0_operators[node - 1];
                for (int row = 0; row < full_dimension; ++row) {
                    const int active_row = active_position[row];
                    if (active_row < 0) continue;
                    for (int entry = op.row_ptr[row];
                         entry < op.row_ptr[row + 1]; ++entry) {
                        const int active_column =
                            active_position[op.column[entry]];
                        if (active_column < 0) continue;
                        const Complex value = op.real_only
                            ? Complex(op.real_value[entry], 0)
                            : op.complex_value[entry];
                        set_complex(temporary, value);
                        gr_ptr destination = gr_mat_entry_ptr(
                            lu_, node * active_count + active_row,
                            node * active_count + active_column,
                            complex_context_
                        );
                        if (gr_sub(
                                destination, destination, temporary,
                                complex_context_
                            ) != GR_SUCCESS) {
                            throw std::runtime_error(
                                "failed to assemble coupled B0 block"
                            );
                        }
                    }
                }
            }
        } catch (...) {
            gr_heap_clear(temporary, complex_context_);
            throw;
        }
        gr_heap_clear(temporary, complex_context_);

        permutation_.resize(dimension_);
        slong rank = 0;
        const int status = gr_mat_lu(
            &rank, permutation_.data(), lu_, lu_, 1, complex_context_
        );
        if (status != GR_SUCCESS || rank != dimension_) {
            throw std::runtime_error("coupled active-support LU failed");
        }
    } catch (...) {
        clear();
        throw;
    }
}

FlintCoupledSolver::~FlintCoupledSolver() { clear(); }

void FlintCoupledSolver::clear() noexcept {
    if (lu_initialized_) {
        gr_mat_clear(lu_, complex_context_);
        lu_initialized_ = false;
    }
    if (complex_context_initialized_) {
        gr_ctx_clear(complex_context_);
        complex_context_initialized_ = false;
    }
    if (real_context_initialized_) {
        gr_ctx_clear(real_context_);
        real_context_initialized_ = false;
    }
}

std::string FlintCoupledSolver::real_string(const Real& value) const {
    return value.str(
        static_cast<std::streamsize>(decimal_precision_ + 20),
        std::ios_base::scientific
    );
}

void FlintCoupledSolver::set_component(
    nfloat_ptr destination, const Real& value) {
    const std::string text = real_string(value);
    if (gr_set_str(destination, text.c_str(), real_context_) != GR_SUCCESS) {
        throw std::runtime_error("failed to convert MPFR value to nfloat");
    }
}

Real FlintCoupledSolver::get_component(nfloat_srcptr source) {
    char* text = nullptr;
    if (gr_get_str_n(
            &text, source, decimal_precision_ + 20, real_context_
        ) != GR_SUCCESS || text == nullptr) {
        if (text != nullptr) flint_free(text);
        throw std::runtime_error("failed to convert nfloat value to MPFR");
    }
    const std::string copy(text);
    flint_free(text);
    return Real(copy);
}

void FlintCoupledSolver::set_complex(gr_ptr destination, const Complex& value) {
    nfloat_complex_ptr entry = static_cast<nfloat_complex_ptr>(destination);
    set_component(NFLOAT_COMPLEX_RE(entry, complex_context_), value.real());
    set_component(NFLOAT_COMPLEX_IM(entry, complex_context_), value.imag());
}

Complex FlintCoupledSolver::get_complex(gr_srcptr source) {
    nfloat_complex_srcptr entry = static_cast<nfloat_complex_srcptr>(source);
    return Complex(
        get_component(NFLOAT_COMPLEX_RE(entry, complex_context_)),
        get_component(NFLOAT_COMPLEX_IM(entry, complex_context_))
    );
}

std::vector<Complex> FlintCoupledSolver::solve(
    const std::vector<Complex>& right_hand_side) {
    if (right_hand_side.size() != static_cast<std::size_t>(dimension_)) {
        throw std::runtime_error("coupled RHS has incompatible dimensions");
    }
    gr_mat_t rhs;
    gr_mat_t solution;
    gr_mat_init(rhs, dimension_, 1, complex_context_);
    gr_mat_init(solution, dimension_, 1, complex_context_);
    try {
        if (gr_mat_zero(rhs, complex_context_) != GR_SUCCESS) {
            throw std::runtime_error("failed to initialise coupled RHS");
        }
        for (int row = 0; row < dimension_; ++row) {
            set_complex(
                gr_mat_entry_ptr(rhs, row, 0, complex_context_),
                right_hand_side[row]
            );
        }
        if (gr_mat_nonsingular_solve_lu_precomp(
                solution, permutation_.data(), lu_, rhs, complex_context_
            ) != GR_SUCCESS) {
            throw std::runtime_error("coupled active-support solve failed");
        }
        std::vector<Complex> result(dimension_);
        for (int row = 0; row < dimension_; ++row) {
            result[row] = get_complex(
                gr_mat_entry_srcptr(solution, row, 0, complex_context_)
            );
        }
        gr_mat_clear(solution, complex_context_);
        gr_mat_clear(rhs, complex_context_);
        return result;
    } catch (...) {
        gr_mat_clear(solution, complex_context_);
        gr_mat_clear(rhs, complex_context_);
        throw;
    }
}
