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
    explicit Backend(const std::string& setupFile) { Load(setupFile); }

    // Propagate one or more independent boundary columns through the same
    // cached collocation system.  Each column has `layers` sequential
    // fake-parameter/epsilon coefficients.  Keeping physical columns in one
    // request avoids repeating setup, WSTP parsing, and matrix storage for the
    // B0-only route used by CHESSv2.
    std::string Run(
        const std::string& boundaryText, int layers, int columns,
        bool shouldCacheStates) {
        if (polynomialMode) {
            throw std::runtime_error(
                "sequential run requested for a polynomial backend"
            );
        }
        if (layers <= 0) throw std::runtime_error("layers must be positive");
        if (columns <= 0) throw std::runtime_error("columns must be positive");

        std::istringstream input(boundaryText);
        int inputDimension = 0;
        int inputLayers = 0;
        int inputColumns = 0;
        input >> inputDimension >> inputLayers >> inputColumns;
        if (!input || inputDimension != dimension || inputLayers != layers ||
            inputColumns != columns) {
            throw std::runtime_error("boundary dimensions do not match cached backend");
        }

        std::vector<Complex> boundary(
            static_cast<std::size_t>(dimension) * layers * columns);
        for (Complex& value : boundary) value = ReadComplex(input);
        if (!input) throw std::runtime_error("failed to parse boundary coefficients");

        const bool boundaryIsReal = std::all_of(
            boundary.begin(), boundary.end(),
            [](const Complex& value) { return value.imag() == 0; });
        const bool systemIsReal = std::all_of(
            operators.begin(), operators.end(),
            [](const SparseOperator& op) { return op.realOnly; });

        const auto start = std::chrono::steady_clock::now();
        cacheStates = shouldCacheStates;
        if (boundaryIsReal && systemIsReal) {
            PropagateReal(boundary, layers, columns);
        } else {
            Propagate(boundary, layers, columns);
        }
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();

        std::ostringstream output;
        output << "CHESSNATIVE1 " << std::setprecision(17) << seconds << ' '
               << layers << ' ' << columns << ' ' << dimension;
        output << std::scientific << std::setprecision(precision);
        for (int column = 0; column < columns; ++column) {
            for (int layer = 0; layer < layers; ++layer) {
                for (int component = 0; component < dimension; ++component) {
                    WriteEndpointValue(output, column, layer, component);
                }
            }
        }
        return output.str();
    }

    // Direct mixed-polynomial transport.  B0 has already been incorporated
    // into one active-support collocation LU; Bp with p>0 only builds the RHS
    // from lower epsilon coefficients.  No auxiliary-delta expansion occurs.
    std::string RunPolynomial(
        const std::string& boundaryText, int layers, int columns,
        bool shouldCacheStates) {
        if (!polynomialMode) {
            throw std::runtime_error(
                "polynomial run requested for a sequential backend"
            );
        }
        if (layers <= 0 || columns <= 0) {
            throw std::runtime_error("polynomial boundary dimensions are empty");
        }

        std::istringstream input(boundaryText);
        int inputDimension = 0;
        int inputLayers = 0;
        int inputColumns = 0;
        input >> inputDimension >> inputLayers >> inputColumns;
        if (!input || inputDimension != dimension || inputLayers != layers ||
            inputColumns != columns) {
            throw std::runtime_error(
                "polynomial boundary dimensions do not match cached backend"
            );
        }
        std::vector<Complex> boundary(
            static_cast<std::size_t>(dimension) * layers * columns
        );
        for (Complex& value : boundary) value = ReadComplex(input);
        if (!input) {
            throw std::runtime_error("failed to parse polynomial boundary");
        }

        const auto start = std::chrono::steady_clock::now();
        cacheStates = shouldCacheStates;
        PropagatePolynomial(boundary, layers, columns);
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start
        ).count();

        std::ostringstream output;
        output << "CHESSPOLY1 " << std::setprecision(17) << seconds << ' '
               << layers << ' ' << columns << ' ' << dimension << ' '
               << active.size();
        output << std::scientific << std::setprecision(precision);
        for (int column = 0; column < columns; ++column) {
            for (int layer = 0; layer < layers; ++layer) {
                for (int component = 0; component < dimension; ++component) {
                    WriteEndpointValue(output, column, layer, component);
                }
            }
        }
        return output.str();
    }

    // Serialize data cached by the most recent Run without repeating the
    // propagation.  mode=0 returns every coefficient at every node (canonical
    // output); mode=1 sums coefficients 0..selectedOrder at each node (the
    // physical delta=1 output needed by B0-only propagation).
    std::string Fetch(int mode, int selectedOrder) const {
        if (lastLayers <= 0 || lastColumns <= 0) {
            throw std::runtime_error("no cached propagation result");
        }
        if (!cacheStates) {
            throw std::runtime_error("the last propagation did not cache node states");
        }
        if (selectedOrder < 0 || selectedOrder >= lastLayers) {
            throw std::runtime_error("selected order is outside cached layers");
        }
        if (mode != 0 && mode != 1) {
            throw std::runtime_error("unknown native fetch mode");
        }

        std::ostringstream output;
        output << "CHESSFETCH1 " << mode << ' ' << selectedOrder << ' '
               << lastLayers << ' ' << lastColumns << ' ' << dimension << ' '
               << nodeCount;
        output << std::scientific << std::setprecision(precision);

        if (mode == 0) {
            for (int column = 0; column < lastColumns; ++column) {
                for (int layer = 0; layer <= selectedOrder; ++layer) {
                    for (int node = 0; node < nodeCount; ++node) {
                        for (int component = 0; component < dimension; ++component) {
                            WriteCachedValue(output, column, layer, component, node);
                        }
                    }
                }
            }
        } else {
            for (int column = 0; column < lastColumns; ++column) {
                for (int node = 0; node < nodeCount; ++node) {
                    for (int component = 0; component < dimension; ++component) {
                        WriteCachedSum(
                            output, column, selectedOrder, component, node);
                    }
                }
            }
        }
        return output.str();
    }

private:
    int precision = 0;
    int dimension = 0;
    int nodeCount = 0;
    std::vector<Real> lu;
    std::vector<int> pivots;
    std::vector<SparseOperator> operators;
    bool polynomialMode = false;
    int polynomialDegree = 0;
    std::vector<std::vector<SparseOperator>> positiveOperators;
    std::vector<int> active;
    std::vector<int> inactive;
    std::unique_ptr<FlintCoupledSolver> coupledSolver;
    int lastLayers = 0;
    int lastColumns = 0;
    bool lastResultIsReal = false;
    bool cacheStates = false;
    std::vector<Real> lastRealStates;
    std::vector<Complex> lastComplexStates;
    std::vector<Real> lastRealEndpoints;
    std::vector<Complex> lastComplexEndpoints;

    static Complex ReadComplex(std::istream& input) {
        std::string realText;
        std::string imagText;
        input >> realText >> imagText;
        if (!input) throw std::runtime_error("invalid complex number in backend input");
        return Complex(Real(realText), Real(imagText));
    }

    std::size_t StateIndex(
        int column, int layer, int component, int node) const {
        return (((static_cast<std::size_t>(column) * lastLayers + layer) *
                 dimension + component) * nodeCount + node);
    }

    std::size_t EndpointIndex(
        int column, int layer, int component) const {
        return ((static_cast<std::size_t>(column) * lastLayers + layer) *
                dimension + component);
    }

    void WriteEndpointValue(
        std::ostringstream& output, int column, int layer,
        int component) const {
        const std::size_t index = EndpointIndex(column, layer, component);
        if (lastResultIsReal) {
            output << ' ' << lastRealEndpoints.at(index) << " 0";
        } else {
            const Complex& value = lastComplexEndpoints.at(index);
            output << ' ' << value.real() << ' ' << value.imag();
        }
    }

    void WriteCachedValue(
        std::ostringstream& output, int column, int layer,
        int component, int node) const {
        const std::size_t index = StateIndex(column, layer, component, node);
        if (lastResultIsReal) {
            output << ' ' << lastRealStates.at(index) << " 0";
        } else {
            const Complex& value = lastComplexStates.at(index);
            output << ' ' << value.real() << ' ' << value.imag();
        }
    }

    void WriteCachedSum(
        std::ostringstream& output, int column, int selectedOrder,
        int component, int node) const {
        if (lastResultIsReal) {
            Real sum(0);
            for (int layer = 0; layer <= selectedOrder; ++layer) {
                sum += lastRealStates.at(
                    StateIndex(column, layer, component, node));
            }
            output << ' ' << sum << " 0";
        } else {
            Complex sum(0);
            for (int layer = 0; layer <= selectedOrder; ++layer) {
                sum += lastComplexStates.at(
                    StateIndex(column, layer, component, node));
            }
            output << ' ' << sum.real() << ' ' << sum.imag();
        }
    }

    void Load(const std::string& setupFile) {
        std::ifstream input(setupFile);
        if (!input) throw std::runtime_error("cannot open CHESS C++ setup file");

        std::string magic;
        input >> magic >> precision >> dimension >> nodeCount;
        if (!input || (magic != "CHESSCPP1" && magic != "CHESSCPP2") ||
            precision <= 0 ||
            dimension <= 0 || nodeCount < 2) {
            throw std::runtime_error("invalid CHESS C++ setup header");
        }
        polynomialMode = magic == "CHESSCPP2";

        Real::default_precision(precision);
        Complex::default_precision(precision);
        Real::thread_default_precision(precision);
        Complex::thread_default_precision(precision);

        lu.resize(static_cast<std::size_t>(nodeCount) * nodeCount);
        for (Real& entry : lu) {
            const Complex value = ReadComplex(input);
            if (value.imag() != 0) {
                throw std::runtime_error("scalar collocation matrix must be real");
            }
            entry = value.real();
        }

        if (!polynomialMode) {
            operators.resize(nodeCount - 1);
            for (SparseOperator& op : operators) ReadOperator(input, op);
        } else {
            int activeCount = 0;
            int flintBits = 0;
            input >> polynomialDegree >> activeCount >> flintBits;
            if (!input || polynomialDegree <= 0 || activeCount <= 0 ||
                activeCount > dimension) {
                throw std::runtime_error("invalid polynomial setup header");
            }
            active.resize(activeCount);
            std::vector<bool> activeMask(dimension, false);
            for (int& component : active) {
                input >> component;
                if (!input || component < 0 || component >= dimension ||
                    activeMask[component]) {
                    throw std::runtime_error("invalid active component list");
                }
                activeMask[component] = true;
            }
            for (int component = 0; component < dimension; ++component) {
                if (!activeMask[component]) inactive.push_back(component);
            }
            std::vector<SparseOperator> b0Operators(nodeCount - 1);
            for (SparseOperator& op : b0Operators) ReadOperator(input, op);
            coupledSolver = std::make_unique<FlintCoupledSolver>(
                precision, nodeCount, dimension, active, lu,
                b0Operators, flintBits
            );

            positiveOperators.resize(polynomialDegree);
            for (auto& operatorSet : positiveOperators) {
                operatorSet.resize(nodeCount - 1);
                for (SparseOperator& op : operatorSet) ReadOperator(input, op);
            }
        }
        if (!input) throw std::runtime_error("truncated CHESS C++ setup file");

        Factorize();
    }

    void ReadOperator(std::istream& input, SparseOperator& op) const {
        int nonzeroCount = 0;
        input >> nonzeroCount;
        if (!input || nonzeroCount < 0) {
            throw std::runtime_error("invalid sparse-operator header");
        }

        struct Entry {
            int row;
            int column;
            Complex value;
        };
        std::vector<Entry> entries;
        entries.reserve(nonzeroCount);
        op.rowPtr.assign(dimension + 1, 0);

        for (int index = 0; index < nonzeroCount; ++index) {
            int row = 0;
            int column = 0;
            input >> row >> column;
            Complex value = ReadComplex(input);
            if (!input || row < 0 || row >= dimension ||
                column < 0 || column >= dimension) {
                throw std::runtime_error("invalid sparse-operator entry");
            }
            entries.push_back({row, column, std::move(value)});
        }

        std::sort(entries.begin(), entries.end(), [](const Entry& left, const Entry& right) {
            return std::pair<int, int>(left.row, left.column) <
                   std::pair<int, int>(right.row, right.column);
        });
        for (const Entry& entry : entries) {
            ++op.rowPtr[entry.row + 1];
            if (entry.value.imag() != 0) op.realOnly = false;
        }
        for (int row = 0; row < dimension; ++row) {
            op.rowPtr[row + 1] += op.rowPtr[row];
        }

        op.column.reserve(entries.size());
        if (op.realOnly) {
            op.realValue.reserve(entries.size());
            for (const Entry& entry : entries) {
                op.column.push_back(entry.column);
                op.realValue.push_back(entry.value.real());
            }
        } else {
            op.complexValue.reserve(entries.size());
            for (const Entry& entry : entries) {
                op.column.push_back(entry.column);
                op.complexValue.push_back(entry.value);
            }
        }
    }

    void Factorize() {
        pivots.resize(nodeCount);
        for (int column = 0; column < nodeCount; ++column) {
            int pivot = column;
            Real largest = abs(lu[static_cast<std::size_t>(column) * nodeCount + column]);
            for (int row = column + 1; row < nodeCount; ++row) {
                const Real candidate = abs(lu[static_cast<std::size_t>(row) * nodeCount + column]);
                if (candidate > largest) {
                    largest = candidate;
                    pivot = row;
                }
            }
            if (largest == 0) throw std::runtime_error("singular scalar collocation matrix");
            pivots[column] = pivot;
            if (pivot != column) {
                for (int j = 0; j < nodeCount; ++j) {
                    std::swap(lu[static_cast<std::size_t>(column) * nodeCount + j],
                              lu[static_cast<std::size_t>(pivot) * nodeCount + j]);
                }
            }

            const Real diagonal = lu[static_cast<std::size_t>(column) * nodeCount + column];
            for (int row = column + 1; row < nodeCount; ++row) {
                Real& multiplier = lu[static_cast<std::size_t>(row) * nodeCount + column];
                multiplier /= diagonal;
                for (int j = column + 1; j < nodeCount; ++j) {
                    lu[static_cast<std::size_t>(row) * nodeCount + j] -=
                        multiplier * lu[static_cast<std::size_t>(column) * nodeCount + j];
                }
            }
        }
    }

    void SolveMany(std::vector<Complex>& rightHandSide) const {
        const int q = nodeCount;
        if (rightHandSide.size() % q != 0) {
            throw std::runtime_error("scalar RHS has incompatible dimensions");
        }
        const int n = static_cast<int>(rightHandSide.size() / q);

#pragma omp parallel
        {
            Real::thread_default_precision(precision);
            Complex::thread_default_precision(precision);
#pragma omp for schedule(static)
            for (int component = 0; component < n; ++component) {
                Complex* vector = rightHandSide.data() + static_cast<std::size_t>(component) * q;
                for (int column = 0; column < q; ++column) {
                    if (pivots[column] != column) std::swap(vector[column], vector[pivots[column]]);
                }
                for (int row = 1; row < q; ++row) {
                    for (int column = 0; column < row; ++column) {
                        vector[row] -= lu[static_cast<std::size_t>(row) * q + column] * vector[column];
                    }
                }
                for (int row = q - 1; row >= 0; --row) {
                    for (int column = row + 1; column < q; ++column) {
                        vector[row] -= lu[static_cast<std::size_t>(row) * q + column] * vector[column];
                    }
                    vector[row] /= lu[static_cast<std::size_t>(row) * q + row];
                }
            }
        }
    }

    void SolveManyReal(std::vector<Real>& rightHandSide) const {
        const int q = nodeCount;
        if (rightHandSide.size() % q != 0) {
            throw std::runtime_error("real scalar RHS has incompatible dimensions");
        }
        const int n = static_cast<int>(rightHandSide.size() / q);

#pragma omp parallel
        {
            Real::thread_default_precision(precision);
#pragma omp for schedule(static)
            for (int component = 0; component < n; ++component) {
                Real* vector = rightHandSide.data() + static_cast<std::size_t>(component) * q;
                for (int column = 0; column < q; ++column) {
                    if (pivots[column] != column) std::swap(vector[column], vector[pivots[column]]);
                }
                for (int row = 1; row < q; ++row) {
                    for (int column = 0; column < row; ++column) {
                        vector[row] -= lu[static_cast<std::size_t>(row) * q + column] * vector[column];
                    }
                }
                for (int row = q - 1; row >= 0; --row) {
                    for (int column = row + 1; column < q; ++column) {
                        vector[row] -= lu[static_cast<std::size_t>(row) * q + column] * vector[column];
                    }
                    vector[row] /= lu[static_cast<std::size_t>(row) * q + row];
                }
            }
        }
    }

    void ApplyOperator(const SparseOperator& op, int node,
                        const std::vector<Complex>& previous,
                        std::vector<Complex>& rightHandSide) const {
        const int q = nodeCount;
        if (op.realOnly) {
            for (int row = 0; row < dimension; ++row) {
                Complex sum(0);
                for (int index = op.rowPtr[row]; index < op.rowPtr[row + 1]; ++index) {
                    sum += op.realValue[index] *
                           previous[static_cast<std::size_t>(op.column[index]) * q + node];
                }
                rightHandSide[static_cast<std::size_t>(row) * q + node] = std::move(sum);
            }
        } else {
            for (int row = 0; row < dimension; ++row) {
                Complex sum(0);
                for (int index = op.rowPtr[row]; index < op.rowPtr[row + 1]; ++index) {
                    sum += op.complexValue[index] *
                           previous[static_cast<std::size_t>(op.column[index]) * q + node];
                }
                rightHandSide[static_cast<std::size_t>(row) * q + node] = std::move(sum);
            }
        }
    }

    void AddOperator(const SparseOperator& op, int node,
                      const std::vector<Complex>& source,
                      std::vector<Complex>& rightHandSide) const {
        const int q = nodeCount;
        if (op.realOnly) {
            for (int row = 0; row < dimension; ++row) {
                Complex sum(0);
                for (int index = op.rowPtr[row];
                     index < op.rowPtr[row + 1]; ++index) {
                    sum += op.realValue[index] *
                           source[static_cast<std::size_t>(
                               op.column[index]) * q + node];
                }
                rightHandSide[static_cast<std::size_t>(row) * q + node] +=
                    sum;
            }
        } else {
            for (int row = 0; row < dimension; ++row) {
                Complex sum(0);
                for (int index = op.rowPtr[row];
                     index < op.rowPtr[row + 1]; ++index) {
                    sum += op.complexValue[index] *
                           source[static_cast<std::size_t>(
                               op.column[index]) * q + node];
                }
                rightHandSide[static_cast<std::size_t>(row) * q + node] +=
                    sum;
            }
        }
    }

    void PropagatePolynomial(
        const std::vector<Complex>& boundary, int layers, int columns) {
        const int q = nodeCount;
        const int n = dimension;
        const int activeCount = static_cast<int>(active.size());
        std::vector<std::vector<Complex>> solved(
            layers, std::vector<Complex>(static_cast<std::size_t>(n) * q)
        );
        lastLayers = layers;
        lastColumns = columns;
        lastResultIsReal = false;
        lastRealStates.clear();
        lastRealEndpoints.clear();
        lastComplexEndpoints.assign(
            static_cast<std::size_t>(columns) * layers * n, Complex(0)
        );
        if (cacheStates) {
            lastComplexStates.assign(
                static_cast<std::size_t>(columns) * layers * n * q,
                Complex(0)
            );
        } else {
            lastComplexStates.clear();
        }

        for (int physicalColumn = 0; physicalColumn < columns;
             ++physicalColumn) {
            for (auto& layer : solved) {
                std::fill(layer.begin(), layer.end(), Complex(0));
            }
            for (int layer = 0; layer < layers; ++layer) {
                std::vector<Complex> rhs(
                    static_cast<std::size_t>(n) * q, Complex(0)
                );
                for (int component = 0; component < n; ++component) {
                    const std::size_t boundaryIndex =
                        (static_cast<std::size_t>(physicalColumn) * layers +
                         layer) * n + component;
                    rhs[static_cast<std::size_t>(component) * q] =
                        boundary[boundaryIndex];
                }

#pragma omp parallel
                {
                    Real::thread_default_precision(precision);
                    Complex::thread_default_precision(precision);
#pragma omp for schedule(static)
                    for (int node = 1; node < q; ++node) {
                        const int maximumPower = std::min(
                            polynomialDegree, layer
                        );
                        for (int power = 1; power <= maximumPower; ++power) {
                            AddOperator(
                                positiveOperators[power - 1][node - 1],
                                node, solved[layer - power], rhs
                            );
                        }
                    }
                }

                std::vector<Complex> current(
                    static_cast<std::size_t>(n) * q, Complex(0)
                );
                std::vector<Complex> activeRhs(
                    static_cast<std::size_t>(q) * activeCount
                );
                for (int node = 0; node < q; ++node) {
                    for (int index = 0; index < activeCount; ++index) {
                        activeRhs[
                            static_cast<std::size_t>(node) * activeCount + index
                        ] = rhs[
                            static_cast<std::size_t>(active[index]) * q + node
                        ];
                    }
                }
                const std::vector<Complex> activeSolution =
                    coupledSolver->Solve(activeRhs);
                for (int node = 0; node < q; ++node) {
                    for (int index = 0; index < activeCount; ++index) {
                        current[
                            static_cast<std::size_t>(active[index]) * q + node
                        ] = activeSolution[
                            static_cast<std::size_t>(node) * activeCount + index
                        ];
                    }
                }

                if (!inactive.empty()) {
                    std::vector<Complex> inactiveRhs(
                        static_cast<std::size_t>(inactive.size()) * q
                    );
                    for (std::size_t index = 0; index < inactive.size();
                         ++index) {
                        for (int node = 0; node < q; ++node) {
                            inactiveRhs[index * q + node] = rhs[
                                static_cast<std::size_t>(inactive[index]) * q +
                                node
                            ];
                        }
                    }
                    SolveMany(inactiveRhs);
                    for (std::size_t index = 0; index < inactive.size();
                         ++index) {
                        for (int node = 0; node < q; ++node) {
                            current[
                                static_cast<std::size_t>(inactive[index]) * q +
                                node
                            ] = inactiveRhs[index * q + node];
                        }
                    }
                }

                solved[layer] = current;
                for (int component = 0; component < n; ++component) {
                    lastComplexEndpoints[EndpointIndex(
                        physicalColumn, layer, component
                    )] = current[static_cast<std::size_t>(component) * q + q - 1];
                    if (cacheStates) {
                        for (int node = 0; node < q; ++node) {
                            lastComplexStates[StateIndex(
                                physicalColumn, layer, component, node
                            )] = current[
                                static_cast<std::size_t>(component) * q + node
                            ];
                        }
                    }
                }
            }
        }
    }

    void Propagate(
        const std::vector<Complex>& boundary, int layers, int columns) {
        const int q = nodeCount;
        const int n = dimension;
        std::vector<Complex> previous(static_cast<std::size_t>(n) * q);
        std::vector<Complex> current(static_cast<std::size_t>(n) * q);
        lastLayers = layers;
        lastColumns = columns;
        lastResultIsReal = false;
        lastRealStates.clear();
        lastRealEndpoints.clear();
        lastComplexEndpoints.assign(
            static_cast<std::size_t>(columns) * layers * n, Complex(0));
        if (cacheStates) {
            lastComplexStates.assign(
                static_cast<std::size_t>(columns) * layers * n * q,
                Complex(0));
        } else {
            lastComplexStates.clear();
        }

        for (int physicalColumn = 0; physicalColumn < columns; ++physicalColumn) {
            std::fill(previous.begin(), previous.end(), Complex(0));
            for (int layer = 0; layer < layers; ++layer) {
                std::fill(current.begin(), current.end(), Complex(0));
                for (int component = 0; component < n; ++component) {
                    const std::size_t boundaryIndex =
                        (static_cast<std::size_t>(physicalColumn) * layers + layer) * n +
                        component;
                    current[static_cast<std::size_t>(component) * q] =
                        boundary[boundaryIndex];
                }
                if (layer > 0) {
#pragma omp parallel
                    {
                        Real::thread_default_precision(precision);
                        Complex::thread_default_precision(precision);
#pragma omp for schedule(static)
                        for (int node = 1; node < q; ++node) {
                            ApplyOperator(
                                operators[node - 1], node, previous, current);
                        }
                    }
                }
                SolveMany(current);
                for (int component = 0; component < n; ++component) {
                    lastComplexEndpoints[EndpointIndex(
                        physicalColumn, layer, component)] =
                        current[static_cast<std::size_t>(component) * q + q - 1];
                    if (cacheStates) {
                        for (int node = 0; node < q; ++node) {
                            lastComplexStates[StateIndex(
                                physicalColumn, layer, component, node)] =
                                current[
                                    static_cast<std::size_t>(component) * q + node];
                        }
                    }
                }
                previous.swap(current);
            }
        }
    }

    void PropagateReal(
        const std::vector<Complex>& boundary, int layers, int columns) {
        const int q = nodeCount;
        const int n = dimension;
        std::vector<Real> previous(static_cast<std::size_t>(n) * q);
        std::vector<Real> current(static_cast<std::size_t>(n) * q);
        lastLayers = layers;
        lastColumns = columns;
        lastResultIsReal = true;
        lastComplexStates.clear();
        lastComplexEndpoints.clear();
        lastRealEndpoints.assign(
            static_cast<std::size_t>(columns) * layers * n, Real(0));
        if (cacheStates) {
            lastRealStates.assign(
                static_cast<std::size_t>(columns) * layers * n * q, Real(0));
        } else {
            lastRealStates.clear();
        }

        for (int physicalColumn = 0; physicalColumn < columns; ++physicalColumn) {
            std::fill(previous.begin(), previous.end(), Real(0));
            for (int layer = 0; layer < layers; ++layer) {
                std::fill(current.begin(), current.end(), Real(0));
                for (int component = 0; component < n; ++component) {
                    const std::size_t boundaryIndex =
                        (static_cast<std::size_t>(physicalColumn) * layers + layer) * n +
                        component;
                    current[static_cast<std::size_t>(component) * q] =
                        boundary[boundaryIndex].real();
                }
                if (layer > 0) {
#pragma omp parallel
                    {
                        Real::thread_default_precision(precision);
#pragma omp for schedule(static)
                        for (int node = 1; node < q; ++node) {
                            const SparseOperator& op = operators[node - 1];
                            for (int row = 0; row < n; ++row) {
                                Real sum(0);
                                for (int index = op.rowPtr[row];
                                     index < op.rowPtr[row + 1]; ++index) {
                                    sum += op.realValue[index] *
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
                SolveManyReal(current);
                for (int component = 0; component < n; ++component) {
                    lastRealEndpoints[EndpointIndex(
                        physicalColumn, layer, component)] =
                        current[static_cast<std::size_t>(component) * q + q - 1];
                    if (cacheStates) {
                        for (int node = 0; node < q; ++node) {
                            lastRealStates[StateIndex(
                                physicalColumn, layer, component, node)] =
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
void Guarded(Function&& function) {
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

void CHESSNativeLoad(const char* setupFile) {
    Guarded([&]() {
        backend = std::make_unique<Backend>(setupFile);
        WSPutSymbol(stdlink, "Null");
    });
}

void CHESSNativeRun(
    const char* boundary, int layers, int columns, int cacheStates) {
    Guarded([&]() {
        if (!backend) throw std::runtime_error("backend is not loaded");
        const std::string result = backend->Run(
            boundary, layers, columns, cacheStates != 0);
        WSPutString(stdlink, result.c_str());
    });
}

void CHESSNativePolynomialRun(
    const char* boundary, int layers, int columns, int cacheStates) {
    Guarded([&]() {
        if (!backend) throw std::runtime_error("backend is not loaded");
        const std::string result = backend->RunPolynomial(
            boundary, layers, columns, cacheStates != 0
        );
        WSPutString(stdlink, result.c_str());
    });
}

void CHESSNativeFetch(int mode, int selectedOrder) {
    Guarded([&]() {
        if (!backend) throw std::runtime_error("backend is not loaded");
        const std::string result = backend->Fetch(mode, selectedOrder);
        WSPutString(stdlink, result.c_str());
    });
}

void CHESSNativeClear() {
    Guarded([&]() {
        backend.reset();
        WSPutSymbol(stdlink, "Null");
    });
}

int main(int argc, char** argv) {
    return WSMain(argc, argv);
}
