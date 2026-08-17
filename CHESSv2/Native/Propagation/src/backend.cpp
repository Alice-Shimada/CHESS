#include "wstp.h"
#include "backend_types.hpp"
#include "flint_coupled_solver.hpp"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


namespace {

class Backend {
public:
    explicit Backend(const std::string& setup_file) { load(setup_file); }

    // Propagate one or more independent boundary columns through the same
    // cached collocation system.  Each column has `layers` sequential
    // fake-parameter/epsilon coefficients.  Keeping physical columns in one
    // request avoids repeating setup, WSTP parsing, and matrix storage for the
    // B0-only route used by CHESSv2.
    std::string run(
        const std::string& boundary_text, int layers, int columns,
        bool cache_states) {
        if (polynomial_mode_) {
            throw std::runtime_error(
                "sequential run requested for a polynomial backend"
            );
        }
        if (layers <= 0) throw std::runtime_error("layers must be positive");
        if (columns <= 0) throw std::runtime_error("columns must be positive");

        std::istringstream input(boundary_text);
        int input_dimension = 0;
        int input_layers = 0;
        int input_columns = 0;
        input >> input_dimension >> input_layers >> input_columns;
        if (!input || input_dimension != dimension_ || input_layers != layers ||
            input_columns != columns) {
            throw std::runtime_error("boundary dimensions do not match cached backend");
        }

        std::vector<Complex> boundary(
            static_cast<std::size_t>(dimension_) * layers * columns);
        for (Complex& value : boundary) value = read_complex(input);
        if (!input) throw std::runtime_error("failed to parse boundary coefficients");

        const bool boundary_is_real = std::all_of(
            boundary.begin(), boundary.end(),
            [](const Complex& value) { return value.imag() == 0; });
        const bool system_is_real = std::all_of(
            operators_.begin(), operators_.end(),
            [](const SparseOperator& op) { return op.real_only; });

        const auto start = std::chrono::steady_clock::now();
        cache_states_ = cache_states;
        if (boundary_is_real && system_is_real) {
            propagate_real(boundary, layers, columns);
        } else {
            propagate(boundary, layers, columns);
        }
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();

        std::ostringstream output;
        output << "CHESSNATIVE1 " << std::setprecision(17) << seconds << ' '
               << layers << ' ' << columns << ' ' << dimension_;
        output << std::scientific << std::setprecision(precision_);
        for (int column = 0; column < columns; ++column) {
            for (int layer = 0; layer < layers; ++layer) {
                for (int component = 0; component < dimension_; ++component) {
                    write_endpoint_value(output, column, layer, component);
                }
            }
        }
        return output.str();
    }

    // Direct mixed-polynomial transport.  B0 has already been incorporated
    // into one active-support collocation LU; Bp with p>0 only builds the RHS
    // from lower epsilon coefficients.  No auxiliary-delta expansion occurs.
    std::string run_polynomial(
        const std::string& boundary_text, int layers, int columns,
        bool cache_states) {
        if (!polynomial_mode_) {
            throw std::runtime_error(
                "polynomial run requested for a sequential backend"
            );
        }
        if (layers <= 0 || columns <= 0) {
            throw std::runtime_error("polynomial boundary dimensions are empty");
        }

        std::istringstream input(boundary_text);
        int input_dimension = 0;
        int input_layers = 0;
        int input_columns = 0;
        input >> input_dimension >> input_layers >> input_columns;
        if (!input || input_dimension != dimension_ || input_layers != layers ||
            input_columns != columns) {
            throw std::runtime_error(
                "polynomial boundary dimensions do not match cached backend"
            );
        }
        std::vector<Complex> boundary(
            static_cast<std::size_t>(dimension_) * layers * columns
        );
        for (Complex& value : boundary) value = read_complex(input);
        if (!input) {
            throw std::runtime_error("failed to parse polynomial boundary");
        }

        const auto start = std::chrono::steady_clock::now();
        cache_states_ = cache_states;
        propagate_polynomial(boundary, layers, columns);
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start
        ).count();

        std::ostringstream output;
        output << "CHESSPOLY1 " << std::setprecision(17) << seconds << ' '
               << layers << ' ' << columns << ' ' << dimension_ << ' '
               << active_.size();
        output << std::scientific << std::setprecision(precision_);
        for (int column = 0; column < columns; ++column) {
            for (int layer = 0; layer < layers; ++layer) {
                for (int component = 0; component < dimension_; ++component) {
                    write_endpoint_value(output, column, layer, component);
                }
            }
        }
        return output.str();
    }

    // Serialize data cached by the most recent run without repeating the
    // propagation.  mode=0 returns every coefficient at every node (canonical
    // output); mode=1 sums coefficients 0..selected_order at each node (the
    // physical delta=1 output needed by B0-only propagation).
    std::string fetch(int mode, int selected_order) const {
        if (last_layers_ <= 0 || last_columns_ <= 0) {
            throw std::runtime_error("no cached propagation result");
        }
        if (!cache_states_) {
            throw std::runtime_error("the last propagation did not cache node states");
        }
        if (selected_order < 0 || selected_order >= last_layers_) {
            throw std::runtime_error("selected order is outside cached layers");
        }
        if (mode != 0 && mode != 1) {
            throw std::runtime_error("unknown native fetch mode");
        }

        std::ostringstream output;
        output << "CHESSFETCH1 " << mode << ' ' << selected_order << ' '
               << last_layers_ << ' ' << last_columns_ << ' ' << dimension_ << ' '
               << node_count_;
        output << std::scientific << std::setprecision(precision_);

        if (mode == 0) {
            for (int column = 0; column < last_columns_; ++column) {
                for (int layer = 0; layer <= selected_order; ++layer) {
                    for (int node = 0; node < node_count_; ++node) {
                        for (int component = 0; component < dimension_; ++component) {
                            write_cached_value(output, column, layer, component, node);
                        }
                    }
                }
            }
        } else {
            for (int column = 0; column < last_columns_; ++column) {
                for (int node = 0; node < node_count_; ++node) {
                    for (int component = 0; component < dimension_; ++component) {
                        write_cached_sum(
                            output, column, selected_order, component, node);
                    }
                }
            }
        }
        return output.str();
    }

private:
    int precision_ = 0;
    int dimension_ = 0;
    int node_count_ = 0;
    std::vector<Real> lu_;
    std::vector<int> pivots_;
    std::vector<SparseOperator> operators_;
    bool polynomial_mode_ = false;
    int polynomial_degree_ = 0;
    std::vector<std::vector<SparseOperator>> positive_operators_;
    std::vector<int> active_;
    std::vector<int> inactive_;
    std::unique_ptr<FlintCoupledSolver> coupled_solver_;
    int last_layers_ = 0;
    int last_columns_ = 0;
    bool last_result_is_real_ = false;
    bool cache_states_ = false;
    std::vector<Real> last_real_states_;
    std::vector<Complex> last_complex_states_;
    std::vector<Real> last_real_endpoints_;
    std::vector<Complex> last_complex_endpoints_;

    static Complex read_complex(std::istream& input) {
        std::string real_text;
        std::string imag_text;
        input >> real_text >> imag_text;
        if (!input) throw std::runtime_error("invalid complex number in backend input");
        return Complex(Real(real_text), Real(imag_text));
    }

    std::size_t state_index(
        int column, int layer, int component, int node) const {
        return (((static_cast<std::size_t>(column) * last_layers_ + layer) *
                 dimension_ + component) * node_count_ + node);
    }

    std::size_t endpoint_index(
        int column, int layer, int component) const {
        return ((static_cast<std::size_t>(column) * last_layers_ + layer) *
                dimension_ + component);
    }

    void write_endpoint_value(
        std::ostringstream& output, int column, int layer,
        int component) const {
        const std::size_t index = endpoint_index(column, layer, component);
        if (last_result_is_real_) {
            output << ' ' << last_real_endpoints_.at(index) << " 0";
        } else {
            const Complex& value = last_complex_endpoints_.at(index);
            output << ' ' << value.real() << ' ' << value.imag();
        }
    }

    void write_cached_value(
        std::ostringstream& output, int column, int layer,
        int component, int node) const {
        const std::size_t index = state_index(column, layer, component, node);
        if (last_result_is_real_) {
            output << ' ' << last_real_states_.at(index) << " 0";
        } else {
            const Complex& value = last_complex_states_.at(index);
            output << ' ' << value.real() << ' ' << value.imag();
        }
    }

    void write_cached_sum(
        std::ostringstream& output, int column, int selected_order,
        int component, int node) const {
        if (last_result_is_real_) {
            Real sum(0);
            for (int layer = 0; layer <= selected_order; ++layer) {
                sum += last_real_states_.at(
                    state_index(column, layer, component, node));
            }
            output << ' ' << sum << " 0";
        } else {
            Complex sum(0);
            for (int layer = 0; layer <= selected_order; ++layer) {
                sum += last_complex_states_.at(
                    state_index(column, layer, component, node));
            }
            output << ' ' << sum.real() << ' ' << sum.imag();
        }
    }

    void load(const std::string& setup_file) {
        std::ifstream input(setup_file);
        if (!input) throw std::runtime_error("cannot open CHESS C++ setup file");

        std::string magic;
        input >> magic >> precision_ >> dimension_ >> node_count_;
        if (!input || (magic != "CHESSCPP1" && magic != "CHESSCPP2") ||
            precision_ <= 0 ||
            dimension_ <= 0 || node_count_ < 2) {
            throw std::runtime_error("invalid CHESS C++ setup header");
        }
        polynomial_mode_ = magic == "CHESSCPP2";

        Real::default_precision(precision_);
        Complex::default_precision(precision_);
        Real::thread_default_precision(precision_);
        Complex::thread_default_precision(precision_);

        lu_.resize(static_cast<std::size_t>(node_count_) * node_count_);
        for (Real& entry : lu_) {
            const Complex value = read_complex(input);
            if (value.imag() != 0) {
                throw std::runtime_error("scalar collocation matrix must be real");
            }
            entry = value.real();
        }

        if (!polynomial_mode_) {
            operators_.resize(node_count_ - 1);
            for (SparseOperator& op : operators_) read_operator(input, op);
        } else {
            int active_count = 0;
            int flint_bits = 0;
            input >> polynomial_degree_ >> active_count >> flint_bits;
            if (!input || polynomial_degree_ <= 0 || active_count <= 0 ||
                active_count > dimension_) {
                throw std::runtime_error("invalid polynomial setup header");
            }
            active_.resize(active_count);
            std::vector<bool> active_mask(dimension_, false);
            for (int& component : active_) {
                input >> component;
                if (!input || component < 0 || component >= dimension_ ||
                    active_mask[component]) {
                    throw std::runtime_error("invalid active component list");
                }
                active_mask[component] = true;
            }
            for (int component = 0; component < dimension_; ++component) {
                if (!active_mask[component]) inactive_.push_back(component);
            }
            std::vector<SparseOperator> b0_operators(node_count_ - 1);
            for (SparseOperator& op : b0_operators) read_operator(input, op);
            coupled_solver_ = std::make_unique<FlintCoupledSolver>(
                precision_, node_count_, dimension_, active_, lu_,
                b0_operators, flint_bits
            );

            positive_operators_.resize(polynomial_degree_);
            for (auto& operators : positive_operators_) {
                operators.resize(node_count_ - 1);
                for (SparseOperator& op : operators) read_operator(input, op);
            }
        }
        if (!input) throw std::runtime_error("truncated CHESS C++ setup file");

        factorize();
    }

    void read_operator(std::istream& input, SparseOperator& op) const {
        int nonzero_count = 0;
        input >> nonzero_count;
        if (!input || nonzero_count < 0) {
            throw std::runtime_error("invalid sparse-operator header");
        }

        struct Entry {
            int row;
            int column;
            Complex value;
        };
        std::vector<Entry> entries;
        entries.reserve(nonzero_count);
        op.row_ptr.assign(dimension_ + 1, 0);

        for (int index = 0; index < nonzero_count; ++index) {
            int row = 0;
            int column = 0;
            input >> row >> column;
            Complex value = read_complex(input);
            if (!input || row < 0 || row >= dimension_ ||
                column < 0 || column >= dimension_) {
                throw std::runtime_error("invalid sparse-operator entry");
            }
            entries.push_back({row, column, std::move(value)});
        }

        std::sort(entries.begin(), entries.end(), [](const Entry& left, const Entry& right) {
            return std::pair<int, int>(left.row, left.column) <
                   std::pair<int, int>(right.row, right.column);
        });
        for (const Entry& entry : entries) {
            ++op.row_ptr[entry.row + 1];
            if (entry.value.imag() != 0) op.real_only = false;
        }
        for (int row = 0; row < dimension_; ++row) {
            op.row_ptr[row + 1] += op.row_ptr[row];
        }

        op.column.reserve(entries.size());
        if (op.real_only) {
            op.real_value.reserve(entries.size());
            for (const Entry& entry : entries) {
                op.column.push_back(entry.column);
                op.real_value.push_back(entry.value.real());
            }
        } else {
            op.complex_value.reserve(entries.size());
            for (const Entry& entry : entries) {
                op.column.push_back(entry.column);
                op.complex_value.push_back(entry.value);
            }
        }
    }

    void factorize() {
        pivots_.resize(node_count_);
        for (int column = 0; column < node_count_; ++column) {
            int pivot = column;
            Real largest = abs(lu_[static_cast<std::size_t>(column) * node_count_ + column]);
            for (int row = column + 1; row < node_count_; ++row) {
                const Real candidate = abs(lu_[static_cast<std::size_t>(row) * node_count_ + column]);
                if (candidate > largest) {
                    largest = candidate;
                    pivot = row;
                }
            }
            if (largest == 0) throw std::runtime_error("singular scalar collocation matrix");
            pivots_[column] = pivot;
            if (pivot != column) {
                for (int j = 0; j < node_count_; ++j) {
                    std::swap(lu_[static_cast<std::size_t>(column) * node_count_ + j],
                              lu_[static_cast<std::size_t>(pivot) * node_count_ + j]);
                }
            }

            const Real diagonal = lu_[static_cast<std::size_t>(column) * node_count_ + column];
            for (int row = column + 1; row < node_count_; ++row) {
                Real& multiplier = lu_[static_cast<std::size_t>(row) * node_count_ + column];
                multiplier /= diagonal;
                for (int j = column + 1; j < node_count_; ++j) {
                    lu_[static_cast<std::size_t>(row) * node_count_ + j] -=
                        multiplier * lu_[static_cast<std::size_t>(column) * node_count_ + j];
                }
            }
        }
    }

    void solve_many(std::vector<Complex>& right_hand_side) const {
        const int q = node_count_;
        if (right_hand_side.size() % q != 0) {
            throw std::runtime_error("scalar RHS has incompatible dimensions");
        }
        const int n = static_cast<int>(right_hand_side.size() / q);

#pragma omp parallel
        {
            Real::thread_default_precision(precision_);
            Complex::thread_default_precision(precision_);
#pragma omp for schedule(static)
            for (int component = 0; component < n; ++component) {
                Complex* vector = right_hand_side.data() + static_cast<std::size_t>(component) * q;
                for (int column = 0; column < q; ++column) {
                    if (pivots_[column] != column) std::swap(vector[column], vector[pivots_[column]]);
                }
                for (int row = 1; row < q; ++row) {
                    for (int column = 0; column < row; ++column) {
                        vector[row] -= lu_[static_cast<std::size_t>(row) * q + column] * vector[column];
                    }
                }
                for (int row = q - 1; row >= 0; --row) {
                    for (int column = row + 1; column < q; ++column) {
                        vector[row] -= lu_[static_cast<std::size_t>(row) * q + column] * vector[column];
                    }
                    vector[row] /= lu_[static_cast<std::size_t>(row) * q + row];
                }
            }
        }
    }

    void solve_many_real(std::vector<Real>& right_hand_side) const {
        const int q = node_count_;
        if (right_hand_side.size() % q != 0) {
            throw std::runtime_error("real scalar RHS has incompatible dimensions");
        }
        const int n = static_cast<int>(right_hand_side.size() / q);

#pragma omp parallel
        {
            Real::thread_default_precision(precision_);
#pragma omp for schedule(static)
            for (int component = 0; component < n; ++component) {
                Real* vector = right_hand_side.data() + static_cast<std::size_t>(component) * q;
                for (int column = 0; column < q; ++column) {
                    if (pivots_[column] != column) std::swap(vector[column], vector[pivots_[column]]);
                }
                for (int row = 1; row < q; ++row) {
                    for (int column = 0; column < row; ++column) {
                        vector[row] -= lu_[static_cast<std::size_t>(row) * q + column] * vector[column];
                    }
                }
                for (int row = q - 1; row >= 0; --row) {
                    for (int column = row + 1; column < q; ++column) {
                        vector[row] -= lu_[static_cast<std::size_t>(row) * q + column] * vector[column];
                    }
                    vector[row] /= lu_[static_cast<std::size_t>(row) * q + row];
                }
            }
        }
    }

    void apply_operator(const SparseOperator& op, int node,
                        const std::vector<Complex>& previous,
                        std::vector<Complex>& right_hand_side) const {
        const int q = node_count_;
        if (op.real_only) {
            for (int row = 0; row < dimension_; ++row) {
                Complex sum(0);
                for (int index = op.row_ptr[row]; index < op.row_ptr[row + 1]; ++index) {
                    sum += op.real_value[index] *
                           previous[static_cast<std::size_t>(op.column[index]) * q + node];
                }
                right_hand_side[static_cast<std::size_t>(row) * q + node] = std::move(sum);
            }
        } else {
            for (int row = 0; row < dimension_; ++row) {
                Complex sum(0);
                for (int index = op.row_ptr[row]; index < op.row_ptr[row + 1]; ++index) {
                    sum += op.complex_value[index] *
                           previous[static_cast<std::size_t>(op.column[index]) * q + node];
                }
                right_hand_side[static_cast<std::size_t>(row) * q + node] = std::move(sum);
            }
        }
    }

    void add_operator(const SparseOperator& op, int node,
                      const std::vector<Complex>& source,
                      std::vector<Complex>& right_hand_side) const {
        const int q = node_count_;
        if (op.real_only) {
            for (int row = 0; row < dimension_; ++row) {
                Complex sum(0);
                for (int index = op.row_ptr[row];
                     index < op.row_ptr[row + 1]; ++index) {
                    sum += op.real_value[index] *
                           source[static_cast<std::size_t>(
                               op.column[index]) * q + node];
                }
                right_hand_side[static_cast<std::size_t>(row) * q + node] +=
                    sum;
            }
        } else {
            for (int row = 0; row < dimension_; ++row) {
                Complex sum(0);
                for (int index = op.row_ptr[row];
                     index < op.row_ptr[row + 1]; ++index) {
                    sum += op.complex_value[index] *
                           source[static_cast<std::size_t>(
                               op.column[index]) * q + node];
                }
                right_hand_side[static_cast<std::size_t>(row) * q + node] +=
                    sum;
            }
        }
    }

    void propagate_polynomial(
        const std::vector<Complex>& boundary, int layers, int columns) {
        const int q = node_count_;
        const int n = dimension_;
        const int active_count = static_cast<int>(active_.size());
        std::vector<std::vector<Complex>> solved(
            layers, std::vector<Complex>(static_cast<std::size_t>(n) * q)
        );
        last_layers_ = layers;
        last_columns_ = columns;
        last_result_is_real_ = false;
        last_real_states_.clear();
        last_real_endpoints_.clear();
        last_complex_endpoints_.assign(
            static_cast<std::size_t>(columns) * layers * n, Complex(0)
        );
        if (cache_states_) {
            last_complex_states_.assign(
                static_cast<std::size_t>(columns) * layers * n * q,
                Complex(0)
            );
        } else {
            last_complex_states_.clear();
        }

        for (int physical_column = 0; physical_column < columns;
             ++physical_column) {
            for (auto& layer : solved) {
                std::fill(layer.begin(), layer.end(), Complex(0));
            }
            for (int layer = 0; layer < layers; ++layer) {
                std::vector<Complex> rhs(
                    static_cast<std::size_t>(n) * q, Complex(0)
                );
                for (int component = 0; component < n; ++component) {
                    const std::size_t boundary_index =
                        (static_cast<std::size_t>(physical_column) * layers +
                         layer) * n + component;
                    rhs[static_cast<std::size_t>(component) * q] =
                        boundary[boundary_index];
                }

#pragma omp parallel
                {
                    Real::thread_default_precision(precision_);
                    Complex::thread_default_precision(precision_);
#pragma omp for schedule(static)
                    for (int node = 1; node < q; ++node) {
                        const int maximum_power = std::min(
                            polynomial_degree_, layer
                        );
                        for (int power = 1; power <= maximum_power; ++power) {
                            add_operator(
                                positive_operators_[power - 1][node - 1],
                                node, solved[layer - power], rhs
                            );
                        }
                    }
                }

                std::vector<Complex> current(
                    static_cast<std::size_t>(n) * q, Complex(0)
                );
                std::vector<Complex> active_rhs(
                    static_cast<std::size_t>(q) * active_count
                );
                for (int node = 0; node < q; ++node) {
                    for (int index = 0; index < active_count; ++index) {
                        active_rhs[
                            static_cast<std::size_t>(node) * active_count + index
                        ] = rhs[
                            static_cast<std::size_t>(active_[index]) * q + node
                        ];
                    }
                }
                const std::vector<Complex> active_solution =
                    coupled_solver_->solve(active_rhs);
                for (int node = 0; node < q; ++node) {
                    for (int index = 0; index < active_count; ++index) {
                        current[
                            static_cast<std::size_t>(active_[index]) * q + node
                        ] = active_solution[
                            static_cast<std::size_t>(node) * active_count + index
                        ];
                    }
                }

                if (!inactive_.empty()) {
                    std::vector<Complex> inactive_rhs(
                        static_cast<std::size_t>(inactive_.size()) * q
                    );
                    for (std::size_t index = 0; index < inactive_.size();
                         ++index) {
                        for (int node = 0; node < q; ++node) {
                            inactive_rhs[index * q + node] = rhs[
                                static_cast<std::size_t>(inactive_[index]) * q +
                                node
                            ];
                        }
                    }
                    solve_many(inactive_rhs);
                    for (std::size_t index = 0; index < inactive_.size();
                         ++index) {
                        for (int node = 0; node < q; ++node) {
                            current[
                                static_cast<std::size_t>(inactive_[index]) * q +
                                node
                            ] = inactive_rhs[index * q + node];
                        }
                    }
                }

                solved[layer] = current;
                for (int component = 0; component < n; ++component) {
                    last_complex_endpoints_[endpoint_index(
                        physical_column, layer, component
                    )] = current[static_cast<std::size_t>(component) * q + q - 1];
                    if (cache_states_) {
                        for (int node = 0; node < q; ++node) {
                            last_complex_states_[state_index(
                                physical_column, layer, component, node
                            )] = current[
                                static_cast<std::size_t>(component) * q + node
                            ];
                        }
                    }
                }
            }
        }
    }

    void propagate(
        const std::vector<Complex>& boundary, int layers, int columns) {
        const int q = node_count_;
        const int n = dimension_;
        std::vector<Complex> previous(static_cast<std::size_t>(n) * q);
        std::vector<Complex> current(static_cast<std::size_t>(n) * q);
        last_layers_ = layers;
        last_columns_ = columns;
        last_result_is_real_ = false;
        last_real_states_.clear();
        last_real_endpoints_.clear();
        last_complex_endpoints_.assign(
            static_cast<std::size_t>(columns) * layers * n, Complex(0));
        if (cache_states_) {
            last_complex_states_.assign(
                static_cast<std::size_t>(columns) * layers * n * q,
                Complex(0));
        } else {
            last_complex_states_.clear();
        }

        for (int physical_column = 0; physical_column < columns; ++physical_column) {
            std::fill(previous.begin(), previous.end(), Complex(0));
            for (int layer = 0; layer < layers; ++layer) {
                std::fill(current.begin(), current.end(), Complex(0));
                for (int component = 0; component < n; ++component) {
                    const std::size_t boundary_index =
                        (static_cast<std::size_t>(physical_column) * layers + layer) * n +
                        component;
                    current[static_cast<std::size_t>(component) * q] =
                        boundary[boundary_index];
                }
                if (layer > 0) {
#pragma omp parallel
                    {
                        Real::thread_default_precision(precision_);
                        Complex::thread_default_precision(precision_);
#pragma omp for schedule(static)
                        for (int node = 1; node < q; ++node) {
                            apply_operator(
                                operators_[node - 1], node, previous, current);
                        }
                    }
                }
                solve_many(current);
                for (int component = 0; component < n; ++component) {
                    last_complex_endpoints_[endpoint_index(
                        physical_column, layer, component)] =
                        current[static_cast<std::size_t>(component) * q + q - 1];
                    if (cache_states_) {
                        for (int node = 0; node < q; ++node) {
                            last_complex_states_[state_index(
                                physical_column, layer, component, node)] =
                                current[
                                    static_cast<std::size_t>(component) * q + node];
                        }
                    }
                }
                previous.swap(current);
            }
        }
    }

    void propagate_real(
        const std::vector<Complex>& boundary, int layers, int columns) {
        const int q = node_count_;
        const int n = dimension_;
        std::vector<Real> previous(static_cast<std::size_t>(n) * q);
        std::vector<Real> current(static_cast<std::size_t>(n) * q);
        last_layers_ = layers;
        last_columns_ = columns;
        last_result_is_real_ = true;
        last_complex_states_.clear();
        last_complex_endpoints_.clear();
        last_real_endpoints_.assign(
            static_cast<std::size_t>(columns) * layers * n, Real(0));
        if (cache_states_) {
            last_real_states_.assign(
                static_cast<std::size_t>(columns) * layers * n * q, Real(0));
        } else {
            last_real_states_.clear();
        }

        for (int physical_column = 0; physical_column < columns; ++physical_column) {
            std::fill(previous.begin(), previous.end(), Real(0));
            for (int layer = 0; layer < layers; ++layer) {
                std::fill(current.begin(), current.end(), Real(0));
                for (int component = 0; component < n; ++component) {
                    const std::size_t boundary_index =
                        (static_cast<std::size_t>(physical_column) * layers + layer) * n +
                        component;
                    current[static_cast<std::size_t>(component) * q] =
                        boundary[boundary_index].real();
                }
                if (layer > 0) {
#pragma omp parallel
                    {
                        Real::thread_default_precision(precision_);
#pragma omp for schedule(static)
                        for (int node = 1; node < q; ++node) {
                            const SparseOperator& op = operators_[node - 1];
                            for (int row = 0; row < n; ++row) {
                                Real sum(0);
                                for (int index = op.row_ptr[row];
                                     index < op.row_ptr[row + 1]; ++index) {
                                    sum += op.real_value[index] *
                                           previous[
                                               static_cast<std::size_t>(
                                                   op.column[index]) * q + node];
                                }
                                current[static_cast<std::size_t>(row) * q + node] =
                                    std::move(sum);
                            }
                        }
                    }
                }
                solve_many_real(current);
                for (int component = 0; component < n; ++component) {
                    last_real_endpoints_[endpoint_index(
                        physical_column, layer, component)] =
                        current[static_cast<std::size_t>(component) * q + q - 1];
                    if (cache_states_) {
                        for (int node = 0; node < q; ++node) {
                            last_real_states_[state_index(
                                physical_column, layer, component, node)] =
                                current[
                                    static_cast<std::size_t>(component) * q + node];
                        }
                    }
                }
                previous.swap(current);
            }
        }
    }
};

std::unique_ptr<Backend> backend;

template <class Function>
void guarded(Function&& function) {
    try {
        function();
    } catch (const std::exception& error) {
        std::cerr << "CHESS C++ backend: " << error.what() << std::endl;
        WSPutSymbol(stdlink, "$Failed");
    } catch (...) {
        std::cerr << "CHESS C++ backend: unknown error" << std::endl;
        WSPutSymbol(stdlink, "$Failed");
    }
}

}  // namespace

void chess_native_load(const char* setup_file) {
    guarded([&]() {
        backend = std::make_unique<Backend>(setup_file);
        WSPutSymbol(stdlink, "Null");
    });
}

void chess_native_run(
    const char* boundary, int layers, int columns, int cache_states) {
    guarded([&]() {
        if (!backend) throw std::runtime_error("backend is not loaded");
        const std::string result = backend->run(
            boundary, layers, columns, cache_states != 0);
        WSPutString(stdlink, result.c_str());
    });
}

void chess_native_polynomial_run(
    const char* boundary, int layers, int columns, int cache_states) {
    guarded([&]() {
        if (!backend) throw std::runtime_error("backend is not loaded");
        const std::string result = backend->run_polynomial(
            boundary, layers, columns, cache_states != 0
        );
        WSPutString(stdlink, result.c_str());
    });
}

void chess_native_fetch(int mode, int selected_order) {
    guarded([&]() {
        if (!backend) throw std::runtime_error("backend is not loaded");
        const std::string result = backend->fetch(mode, selected_order);
        WSPutString(stdlink, result.c_str());
    });
}

void chess_native_clear() {
    guarded([&]() {
        backend.reset();
        WSPutSymbol(stdlink, "Null");
    });
}

int main(int argc, char** argv) {
    return WSMain(argc, argv);
}
