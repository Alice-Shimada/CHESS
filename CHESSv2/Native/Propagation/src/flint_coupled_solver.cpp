#include "flint_coupled_solver.hpp"

#include <flint/flint.h>

#include <ios>
#include <stdexcept>

FlintCoupledSolver::FlintCoupledSolver(
    int requestedPrecision,
    int nodeCount,
    int fullDimension,
    const std::vector<int>& active,
    const std::vector<Real>& scalarMatrix,
    const std::vector<SparseOperator>& b0Operators,
    int flintBits
) : decimalPrecision(requestedPrecision),
    dimension(nodeCount * static_cast<int>(active.size())) {
    try {
        if (flintBits <= 0 || dimension <= 0 ||
            scalarMatrix.size() !=
                static_cast<std::size_t>(nodeCount) * nodeCount ||
            b0Operators.size() != static_cast<std::size_t>(nodeCount - 1)) {
            throw std::runtime_error("invalid coupled-solver dimensions");
        }
        if (nfloat_ctx_init(realContext, flintBits, 0) != GR_SUCCESS) {
            throw std::runtime_error("failed to initialise real nfloat context");
        }
        realContextInitialized = true;
        if (nfloat_complex_ctx_init(complexContext, flintBits, 0) !=
            GR_SUCCESS) {
            throw std::runtime_error(
                "failed to initialise complex nfloat context"
            );
        }
        complexContextInitialized = true;
        gr_mat_init(lu, dimension, dimension, complexContext);
        luInitialized = true;
        if (gr_mat_zero(lu, complexContext) != GR_SUCCESS) {
            throw std::runtime_error("failed to zero coupled matrix");
        }

        const int activeCount = static_cast<int>(active.size());
        for (int rowNode = 0; rowNode < nodeCount; ++rowNode) {
            for (int columnNode = 0; columnNode < nodeCount; ++columnNode) {
                const Real& derivative = scalarMatrix[
                    static_cast<std::size_t>(rowNode) * nodeCount + columnNode
                ];
                if (derivative == 0) continue;
                for (int index = 0; index < activeCount; ++index) {
                    SetComplex(
                        gr_mat_entry_ptr(
                            lu, rowNode * activeCount + index,
                            columnNode * activeCount + index, complexContext
                        ),
                        Complex(derivative, 0)
                    );
                }
            }
        }

        std::vector<int> activePosition(fullDimension, -1);
        for (int index = 0; index < activeCount; ++index) {
            activePosition[active[index]] = index;
        }
        gr_ptr temporary = gr_heap_init(complexContext);
        try {
            for (int node = 1; node < nodeCount; ++node) {
                const SparseOperator& op = b0Operators[node - 1];
                for (int row = 0; row < fullDimension; ++row) {
                    const int activeRow = activePosition[row];
                    if (activeRow < 0) continue;
                    for (int entry = op.rowPtr[row];
                         entry < op.rowPtr[row + 1]; ++entry) {
                        const int activeColumn =
                            activePosition[op.column[entry]];
                        if (activeColumn < 0) continue;
                        const Complex value = op.realOnly
                            ? Complex(op.realValue[entry], 0)
                            : op.complexValue[entry];
                        SetComplex(temporary, value);
                        gr_ptr destination = gr_mat_entry_ptr(
                            lu, node * activeCount + activeRow,
                            node * activeCount + activeColumn,
                            complexContext
                        );
                        if (gr_sub(
                                destination, destination, temporary,
                                complexContext
                            ) != GR_SUCCESS) {
                            throw std::runtime_error(
                                "failed to assemble coupled B0 block"
                            );
                        }
                    }
                }
            }
        } catch (...) {
            gr_heap_clear(temporary, complexContext);
            throw;
        }
        gr_heap_clear(temporary, complexContext);

        permutation.resize(dimension);
        slong rank = 0;
        const int status = gr_mat_lu(
            &rank, permutation.data(), lu, lu, 1, complexContext
        );
        if (status != GR_SUCCESS || rank != dimension) {
            throw std::runtime_error("coupled active-support LU failed");
        }
    } catch (...) {
        Clear();
        throw;
    }
}

FlintCoupledSolver::~FlintCoupledSolver() { Clear(); }

void FlintCoupledSolver::Clear() noexcept {
    if (luInitialized) {
        gr_mat_clear(lu, complexContext);
        luInitialized = false;
    }
    if (complexContextInitialized) {
        gr_ctx_clear(complexContext);
        complexContextInitialized = false;
    }
    if (realContextInitialized) {
        gr_ctx_clear(realContext);
        realContextInitialized = false;
    }
}

std::string FlintCoupledSolver::RealString(const Real& value) const {
    return value.str(
        static_cast<std::streamsize>(decimalPrecision + 20),
        std::ios_base::scientific
    );
}

void FlintCoupledSolver::SetComponent(
    nfloat_ptr destination, const Real& value) {
    const std::string text = RealString(value);
    if (gr_set_str(destination, text.c_str(), realContext) != GR_SUCCESS) {
        throw std::runtime_error("failed to convert MPFR value to nfloat");
    }
}

Real FlintCoupledSolver::GetComponent(nfloat_srcptr source) {
    char* text = nullptr;
    if (gr_get_str_n(
            &text, source, decimalPrecision + 20, realContext
        ) != GR_SUCCESS || text == nullptr) {
        if (text != nullptr) flint_free(text);
        throw std::runtime_error("failed to convert nfloat value to MPFR");
    }
    const std::string copy(text);
    flint_free(text);
    return Real(copy);
}

void FlintCoupledSolver::SetComplex(gr_ptr destination, const Complex& value) {
    nfloat_complex_ptr entry = static_cast<nfloat_complex_ptr>(destination);
    SetComponent(NFLOAT_COMPLEX_RE(entry, complexContext), value.real());
    SetComponent(NFLOAT_COMPLEX_IM(entry, complexContext), value.imag());
}

Complex FlintCoupledSolver::GetComplex(gr_srcptr source) {
    nfloat_complex_srcptr entry = static_cast<nfloat_complex_srcptr>(source);
    return Complex(
        GetComponent(NFLOAT_COMPLEX_RE(entry, complexContext)),
        GetComponent(NFLOAT_COMPLEX_IM(entry, complexContext))
    );
}

std::vector<Complex> FlintCoupledSolver::Solve(
    const std::vector<Complex>& rightHandSide) {
    if (rightHandSide.size() != static_cast<std::size_t>(dimension)) {
        throw std::runtime_error("coupled RHS has incompatible dimensions");
    }
    gr_mat_t rhs;
    gr_mat_t solution;
    gr_mat_init(rhs, dimension, 1, complexContext);
    gr_mat_init(solution, dimension, 1, complexContext);
    try {
        if (gr_mat_zero(rhs, complexContext) != GR_SUCCESS) {
            throw std::runtime_error("failed to initialise coupled RHS");
        }
        for (int row = 0; row < dimension; ++row) {
            SetComplex(
                gr_mat_entry_ptr(rhs, row, 0, complexContext),
                rightHandSide[row]
            );
        }
        if (gr_mat_nonsingular_solve_lu_precomp(
                solution, permutation.data(), lu, rhs, complexContext
            ) != GR_SUCCESS) {
            throw std::runtime_error("coupled active-support solve failed");
        }
        std::vector<Complex> result(dimension);
        for (int row = 0; row < dimension; ++row) {
            result[row] = GetComplex(
                gr_mat_entry_srcptr(solution, row, 0, complexContext)
            );
        }
        gr_mat_clear(solution, complexContext);
        gr_mat_clear(rhs, complexContext);
        return result;
    } catch (...) {
        gr_mat_clear(solution, complexContext);
        gr_mat_clear(rhs, complexContext);
        throw;
    }
}
