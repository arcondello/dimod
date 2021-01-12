// Copyright 2020 D-Wave Systems Inc.
//
//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at
//
//        http://www.apache.org/licenses/LICENSE-2.0
//
//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

#pragma once

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dimod/utils.h"

namespace dimod {

template <class Bias>
class AdjVectorBQM2 {
 public:
    using bias_type = Bias;
    using size_type = std::size_t;
    using variable_type = std::size_t;

    using neighborhood_type =
            typename std::vector<std::pair<variable_type, bias_type>>;
    using neighborhood_iterator = typename neighborhood_type::iterator;
    using const_neighborhood_iterator =
            typename neighborhood_type::const_iterator;

    bias_type offset;

    AdjVectorBQM2(): offset(0) {}

    /**
     * Construct a BQM from a dense array.
     *
     * @param dense An array containing the biases. Assumed to contain
     *     `num_variables`^2 elements. The upper and lower triangle are summed.
     * @param num_variables The number of variables.
     */
    template <class T>
    AdjVectorBQM2(const T dense[], size_type num_variables,
                  bool ignore_diagonal = false): offset(0) {
        // we know how big our linear is going to be
        adj_.resize(num_variables);

        bias_type qbias;

        if (!ignore_diagonal) {
            for (size_type v = 0; v < num_variables; ++v) {
                adj_[v].second = dense[v * (num_variables + 1)];
            }
        }

        for (size_type u = 0; u < num_variables; ++u) {
            for (size_type v = u + 1; v < num_variables; ++v) {
                qbias = dense[u * num_variables + v] +
                        dense[v * num_variables + u];

                if (qbias != 0) {
                    adj_[u].first.emplace_back(v, qbias);
                    adj_[v].first.emplace_back(u, qbias);
                }
            }
        }
    }

    /**
     * Construct a BQM from COO-formated iterators.
     *
     * A sparse BQM encoded in [COOrdinate] format is specified by three
     * arrays of (row, column, value).
     *
     * [COOrdinate]: https://w.wiki/n$L
     *
     * @param row_iterator Iterator pointing to the beginning of the row data.
     *     Must be a random access iterator.
     * @param col_iterator Iterator pointing to the beginning of the column
     *     data. Must be a random access iterator.
     * @param bias_iterator Iterator pointing to the beginning of the bias data.
     *     Must be a random access iterator.
     * @param length The number of (row, column, bias) entries.
     * @param ignore_diagonal If true, entries on the diagonal of the sparse
     *     matrix are ignored.
     */
    template <class ItRow, class ItCol, class ItBias>
    AdjVectorBQM2(ItRow row_iterator, ItCol col_iterator, ItBias bias_iterator,
                 size_type length, bool ignore_diagonal = false): offset(0) {
        // determine the number of variables so we can allocate adj
        if (length > 0) {
            size_type max_label = std::max(
                    *std::max_element(row_iterator, row_iterator + length),
                    *std::max_element(col_iterator, col_iterator + length));
            adj_.resize(max_label + 1);
        }

        // Count the degrees and use that to reserve the neighborhood vectors
        std::vector<size_type> degrees(adj_.size());
        ItRow rit(row_iterator);
        ItRow cit(col_iterator);
        for (size_type i = 0; i < length; ++i, ++rit, ++cit) {
            if (*rit != *cit) {
                degrees[*rit] += 1;
                degrees[*cit] += 1;
            }
        }
        for (size_type i = 0; i < degrees.size(); ++i) {
            adj_[i].first.reserve(degrees[i]);
        }

        // add the values to the adjacency, not worrying about order or
        // duplicates
        for (size_type i = 0; i < length; i++) {
            if (*row_iterator == *col_iterator) {
                // linear bias
                if (!ignore_diagonal) {
                    linear(*row_iterator) += *bias_iterator;
                }
            } else {
                // quadratic bias
                adj_[*row_iterator].first.emplace_back(*col_iterator,
                                                      *bias_iterator);
                adj_[*col_iterator].first.emplace_back(*row_iterator,
                                                      *bias_iterator);
            }

            ++row_iterator;
            ++col_iterator;
            ++bias_iterator;
        }

        // now sort each neighborhood and remove duplicates
        for (variable_type v = 0; v < adj_.size(); ++v) {
            auto begin = adj_[v].first.begin();
            auto end = adj_[v].first.end();
            if (!std::is_sorted(begin, end)) {
                std::sort(begin, end);
            }

            // now remove any duplicate variables, adding the biases
            auto it = adj_[v].first.begin();
            while (it + 1 < adj_[v].first.end()) {
                if (it->first == (it + 1)->first) {
                    it->second += (it + 1)->second;
                    adj_[v].first.erase(it + 1);
                } else {
                    ++it;
                }
            }
        }
    }

    /**
     * Add quadratic bias to the interaction `u`, `v`.
     *
     * @param u A variable in the binary quadratic model.
     * @param v A variable in the binary quadratic model.
     * @param b The bias to be added to the interaction between `u` and `v`.
     *
     * The behavior of this function is undefined when either `u` or `v` is
     * less than 0 or greater than or equal to the number of variables in the
     * BQM.
     */
    void add_interaction(variable_type u, variable_type v, bias_type b = 0) {
        auto uv = find_quadratic(u, v);
        if (uv.second) {
            uv.first->second += b;
        } else {
            adj_[u].first.emplace(uv.first, v, b);
        }
        auto vu = find_quadratic(v, u);
        if (vu.second) {
            vu.first->second += b;
        } else {
            adj_[v].first.emplace(vu.first, u, b);
        }
    }

    /// Add one (disconnected) variable to the BQM and return its index.
    variable_type add_variable() {
        adj_.resize(adj_.size() + 1);
        return adj_.size() - 1;
    }

    /// Return the degree of variable `v`
    size_type degree(variable_type v) const { return adj_[v].first.size(); }

    /**
     * Return the energy of the given sample.
     *
     * @param Iter A random access iterator pointing to the beginning of the
     *     sample.
     * @return The energy of the given sample.
     *
     * The behavior of this function is undefined when the sample is not
     * <num_variables>"()" long.
     */
    template <class Iter>
    bias_type energy(Iter sample_start) {
        bias_type energy = offset;
        for (variable_type u = 0; u < num_variables(); ++u) {
            auto u_val = *(sample_start + u);

            energy += u_val * linear(u);

            // only look at the lower triangle
            auto span = neighborhood(u);
            while (span.first != span.second && span.first->first < u) {
                variable_type v = span.first->first;
                auto v_val = *(sample_start + v);

                energy += u_val * v_val * span.first->second;

                ++span.first;
            }
        }
    }

    /// Return a reference to the linear bias associated with `v`.
    bias_type& linear(variable_type v) { return adj_[v].second; }

    /// Return a reference to the linear bias associated with `v`.
    const bias_type& linear(variable_type v) const { return adj_[v].second; }

    /// Return a pair of iterators - the start and end of the neighborhood
    std::pair<const_neighborhood_iterator, const_neighborhood_iterator>
    neighborhood(variable_type u) const {
        return std::make_pair(adj_[u].first.cbegin(), adj_[u].first.cend());
    }

    /// Return the number of interactions in the binary quadratic model.
    size_type num_interactions() const {
        size_type count = 0;
        bool remainder = 0;
        for (auto it = adj_.begin(); it != adj_.end(); ++it) {
            if (it->first.size() % 2 == 0) {
                count += it->first.size() / 2;
            } else {
                count += (it->first.size() - 1) / 2;

                if (remainder) {
                    remainder = false;
                    count += 1;
                } else {
                    remainder = true;
                }
            }
        }
        return count;
    }

    /// Return the number of variables in the binary quadratic model.
    size_type num_variables() const { return adj_.size(); }

    /**
     * Return the quadratic bias associated with `u`, `v`.
     *
     * Note that this function does not return a reference, this is because
     * each quadratic bias is stored twice.
     *
     * @param u A variable in the binary quadratic model.
     * @param v A variable in the binary quadratic model.
     * @return The quadratic bias if it exists, otherwise 0.
     */
    bias_type quadratic(variable_type u, variable_type v) const {
        auto low = find_quadratic(u, v);
        if (low.second) {
            return low.first->second;
        } else {
            return 0;
        }
    }

    /**
     * Return the quadratic bias associated with `u`, `v`.
     *
     * Note that this function does not return a reference, this is because
     * each quadratic bias is stored twice.
     *
     * @param u A variable in the binary quadratic model.
     * @param v A variable in the binary quadratic model.
     * @return The quadratic bias if it exists, otherwise 0.
     * @exception out_of_range If either `u` or `v` are not variables or if
     *     they do not have an interaction then the function throws an
     *     exception.
     */
    bias_type quadratic_at(variable_type u, variable_type v) const {
        if (!(u >= 0 && u < adj_.size() && v >= 0 && v < adj_.size()))
            throw std::out_of_range("out of range variable");

        auto low = find_quadratic(u, v);
        if (low.second) {
            return low.first->second;
        } else {
            throw std::out_of_range("given variables have no interaction");
        }
    }

    /// Remove the interaction between `u` and `v`. Return true if it existed.
    bool remove_interaction(variable_type u, variable_type v) {
        auto uv = find_quadratic(u, v);
        if (uv.second) {
            adj_[u].first.erase(uv.first);
            // the mirror must exist
            adj_[v].first.erase(quadratic_lb(v, u));
        }
        return uv.second;
    }

    /// Resize the binary quadratic model to contain n variables.
    void resize(size_type n) {
        if (n < adj_.size()) {
            // Clean out any of the to-be-deleted variables from the
            // neighborhoods.
            // This approach is better in the dense case. In the sparse case
            // we could determine which neighborhoods need to be trimmed rather
            // than just doing them all.
            for (size_type v = 0; v < n; ++v) {
                auto low = quadratic_lb(v, n);
                adj_[v].first.erase(low, adj_[v].first.end());
            }
        }

        adj_.resize(n);
    }

    /**
     * Set the quadratic bias of the interaction `u`, `v`.
     *
     * @param u A variable in the binary quadratic model.
     * @param v A variable in the binary quadratic model.
     * @param b The bias.
     *
     * The behavior of this function is undefined when either `u` or `v` is
     * less than 0 or greater than or equal to the number of variables in the
     * BQM.
     */
    void set_quadratic(variable_type u, variable_type v, bias_type b) {
        auto uv = find_quadratic(u, v);
        if (uv.second) {
            uv.first->second = b;
        } else {
            adj_[u].first.emplace(uv.first, v, b);
        }
        auto vu = find_quadratic(v, u);
        if (vu.second) {
            vu.first->second = b;
        } else {
            adj_[v].first.emplace(vu.first, u, b);
        }
    }

 private:
    std::vector<std::pair<neighborhood_type, bias_type>> adj_;

    std::pair<neighborhood_iterator, bool> find_quadratic(variable_type u,
                                                          variable_type v) {
        auto low = quadratic_lb(u, v);
        return std::make_pair(low,
                              low != adj_[u].first.end() && low->first == v);
    }

    std::pair<const_neighborhood_iterator, bool> find_quadratic(
            variable_type u, variable_type v) const {
        auto low = quadratic_lb(u, v);
        return std::make_pair(low,
                              low != adj_[u].first.cend() && low->first == v);
    }

    inline neighborhood_iterator quadratic_lb(variable_type u,
                                              variable_type v) {
        return std::lower_bound(adj_[u].first.begin(), adj_[u].first.end(), v,
                                utils::comp_v<variable_type, bias_type>);
    }

    inline const_neighborhood_iterator quadratic_lb(variable_type u,
                                                    variable_type v) const {
        return std::lower_bound(adj_[u].first.cbegin(), adj_[u].first.cend(), v,
                                utils::comp_v<variable_type, bias_type>);
    }
};

}  // namespace dimod
