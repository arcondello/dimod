// Copyright 2021 D-Wave Systems Inc.
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
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dimod/abc.h"
#include "dimod/neighborhood.h"
#include "dimod/vartypes.h"

namespace dimod {

// template <class Bias, class Index = int>
// class QuadraticModelBase {
//  public:
//     /// The first template parameter (Bias)
//     using bias_type = Bias;

//     /// The second template parameter (Index).
//     using index_type = Index;

//     /// Unsigned integral type that can represent non-negative values.
//     using size_type = std::size_t;

//     /// A random access iterator to `pair<const index_type&, const bias_type&>`
//     using const_neighborhood_iterator =
//             typename Neighborhood<bias_type, index_type>::const_iterator;

//     /// A forward iterator pointing to the quadratic biases.
//     using const_quadratic_iterator = ConstQuadraticIterator<Bias, Index>;

//     template <class B, class I>
//     friend class BinaryQuadraticModel;

//     QuadraticModelBase() : offset_(0) {}

//     /// Add quadratic bias for the given variables.
//     void add_quadratic(index_type u, index_type v, bias_type bias) {
//         assert(0 <= u && static_cast<size_t>(u) < this->num_variables());
//         assert(0 <= v && static_cast<size_t>(v) < this->num_variables());

//         if (u == v) {
//             switch (this->vartype(u)) {
//                 case Vartype::BINARY: {
//                     // 1*1 == 1 and 0*0 == 0 so this is linear
//                     this->linear(u) += bias;
//                     break;
//                 }
//                 case Vartype::SPIN: {
//                     // -1*-1 == +1*+1 == 1 so this is a constant offset
//                     this->offset_ += bias;
//                     break;
//                 }
//                 default: {
//                     // self-loop
//                     this->adj_[u][v] += bias;
//                     break;
//                 }
//             }
//         } else {
//             this->adj_[u][v] += bias;
//             this->adj_[v][u] += bias;
//         }
//     }

//     /*
//      * Add quadratic biases from a dense matrix.
//      *
//      * `dense` must be an array of length `num_variables^2`.
//      *
//      * Values on the diagonal are treated differently depending on the variable
//      * type.
//      *
//      * # Exceptions
//      * The behavior of this method is undefined when the model has fewer than
//      * `num_variables` variables.
//      */
//     template <class T>
//     void add_quadratic(const T dense[], index_type num_variables) {
//         assert(0 <= num_variables && static_cast<size_t>(num_variables) <= this->num_variables());

//         if (this->is_linear()) {
//             for (index_type u = 0; u < num_variables; ++u) {
//                 // diagonal
//                 this->add_quadratic_back(u, u, dense[u * (num_variables + 1)]);

//                 // off-diagonal
//                 for (index_type v = u + 1; v < num_variables; ++v) {
//                     bias_type qbias = dense[u * num_variables + v] + dense[v * num_variables + u];

//                     if (qbias) {
//                         this->add_quadratic_back(u, v, qbias);
//                     }
//                 }
//             }
//         } else {
//             // we cannot rely on the ordering
//             for (index_type u = 0; u < num_variables; ++u) {
//                 // diagonal
//                 this->add_quadratic(u, u, dense[u * (num_variables + 1)]);

//                 // off-diagonal
//                 for (index_type v = u + 1; v < num_variables; ++v) {
//                     bias_type qbias = dense[u * num_variables + v] + dense[v * num_variables + u];

//                     if (qbias) {
//                         this->add_quadratic(u, v, qbias);
//                     }
//                 }
//             }
//         }
//     }

//     /**
//      * Add quadratic bias for the given variables at the end of eachother's neighborhoods.
//      *
//      * # Parameters
//      * - `u` - a variable.
//      * - `v` - a variable.
//      * - `bias` - the quadratic bias associated with `u` and `v`.
//      *
//      * # Exceptions
//      * When `u` is less than the largest neighbor in `v`'s neighborhood,
//      * `v` is less than the largest neighbor in `u`'s neighborhood, or either
//      * `u` or `v` is greater than ``num_variables()`` then the behavior of
//      * this method is undefined.
//      */
//     void add_quadratic_back(index_type u, index_type v, bias_type bias) {
//         assert(0 <= u && static_cast<size_t>(u) <= this->num_variables());
//         assert(0 <= v && static_cast<size_t>(v) <= this->num_variables());

//         // check the condition for adding at the back
//         assert(this->adj_[v].empty() || this->adj_[v].back().first <= u);
//         assert(this->adj_[u].empty() || this->adj_[u].back().first <= v);

//         if (u == v) {
//             switch (this->vartype(u)) {
//                 case Vartype::BINARY: {
//                     // 1*1 == 1 and 0*0 == 0 so this is linear
//                     this->linear(u) += bias;
//                     break;
//                 }
//                 case Vartype::SPIN: {
//                     // -1*-1 == +1*+1 == 1 so this is a constant offset
//                     this->offset() += bias;
//                     break;
//                 }
//                 default: {
//                     // self-loop
//                     this->adj_[u].emplace_back(v, bias);
//                     break;
//                 }
//             }
//         } else {
//             this->adj_[u].emplace_back(v, bias);
//             this->adj_[v].emplace_back(u, bias);
//         }
//     }

//     /// Remove the offset and all variables and interactions from the model.
//     void clear() {
//         this->adj_.clear();
//         this->linear_biases_.clear();
//         this->offset_ = 0;
//     }

//     /// Return True if the model has no quadratic biases.
//     bool is_linear() const {
//         for (auto it = adj_.begin(); it != adj_.end(); ++it) {
//             if ((*it).size()) {
//                 return false;
//             }
//         }
//         return true;
//     }

//     const_quadratic_iterator cbegin_quadratic() const { return const_quadratic_iterator(this, 0); }

//     const_quadratic_iterator cend_quadratic() const {
//         return const_quadratic_iterator(this, this->num_variables());
//     }

//     /**
//      * Return the energy of the given sample.
//      *
//      * The `sample_start` must be random access iterator pointing to the
//      * beginning of the sample.
//      *
//      * The behavior of this function is undefined when the sample is not
//      * `num_variables()` long.
//      */
//     template <class Iter>  // todo: allow different return types
//     bias_type energy(Iter sample_start) const {
//         bias_type en = offset();

//         for (index_type u = 0; u < static_cast<index_type>(num_variables()); ++u) {
//             auto u_val = *(sample_start + u);

//             en += u_val * linear(u);

//             auto span = neighborhood(u);
//             for (auto nit = span.first; nit != span.second && (*nit).first <= u; ++nit) {
//                 auto v_val = *(sample_start + (*nit).first);
//                 auto bias = (*nit).second;
//                 en += bias * u_val * v_val;
//             }
//         }

//         return en;
//     }

//     /**
//      * Return the energy of the given sample.
//      *
//      * The `sample` must be a vector containing the sample.
//      *
//      * The behavior of this function is undefined when the sample is not
//      * `num_variables()` long.
//      */
//     template <class T>
//     bias_type energy(const std::vector<T>& sample) const {
//         // todo: check length?
//         return energy(sample.cbegin());
//     }

//     virtual const bias_type& lower_bound(index_type) const = 0;

//     /// Return a reference to the linear bias associated with `v`.
//     bias_type& linear(index_type v) { return linear_biases_[v]; }

//     /// Return a reference to the linear bias associated with `v`.
//     const bias_type& linear(index_type v) const { return linear_biases_[v]; }

//     /// Return a pair of iterators - the start and end of the neighborhood
//     std::pair<const_neighborhood_iterator, const_neighborhood_iterator> neighborhood(
//             index_type u) const {
//         return std::make_pair(adj_[u].cbegin(), adj_[u].cend());
//     }

//     *
//      * The neighborhood of variable `v`.
//      *
//      * @param A variable `v`.
//      * @param The neighborhood will start with the first out variable that
//      * does not compare less than `start`.
//      *
//      * @returns A pair of iterators pointing to the start and end of the
//      *     neighborhood.
     
//     std::pair<const_neighborhood_iterator, const_neighborhood_iterator> neighborhood(
//             index_type u, index_type start) const {
//         return std::make_pair(adj_[u].lower_bound(start), adj_[u].cend());
//     }

//     /**
//      * Return the quadratic bias associated with `u`, `v`.
//      *
//      * If `u` and `v` do not have a quadratic bias, returns 0.
//      *
//      * Note that this function does not return a reference, this is because
//      * each quadratic bias is stored twice.
//      *
//      */
//     bias_type quadratic(index_type u, index_type v) const { return adj_[u].get(v); }

//     /**
//      * Return the quadratic bias associated with `u`, `v`.
//      *
//      * Note that this function does not return a reference, this is because
//      * each quadratic bias is stored twice.
//      *
//      * Raises an `out_of_range` error if either `u` or `v` are not variables or
//      * if they do not have an interaction then the function throws an exception.
//      */
//     bias_type quadratic_at(index_type u, index_type v) const { return adj_[u].at(v); }

//     /**
//      * Total bytes consumed by the biases and indices.
//      *
//      * If `capacity` is true, use the capacity of the underlying vectors rather
//      * than the size.
//      */
//     size_type nbytes(bool capacity = false) const noexcept {
//         size_type count = sizeof(bias_type);  // offset
//         if (capacity) {
//             count += this->linear_biases_.capacity() * sizeof(bias_type);
//         } else {
//             count += this->linear_biases_.size() * sizeof(bias_type);
//         }
//         for (size_type v = 0; v < this->num_variables(); ++v) {
//             count += this->adj_[v].nbytes(capacity);
//         }
//         return count;
//     }

//     /// Return the number of variables in the quadratic model.
//     size_type num_variables() const { return linear_biases_.size(); }

//     /**
//      *  Return the number of interactions in the quadratic model.
//      *
//      * `O(num_variables*log(num_variables))` complexity.
//      */
//     size_type num_interactions() const {
//         size_type count = 0;
//         for (index_type v = 0; v < static_cast<index_type>(this->num_variables()); ++v) {
//             count += this->adj_[v].size();

//             // account for self-loops
//             auto lb = this->adj_[v].lower_bound(v);
//             if (lb != this->adj_[v].cend() && lb->first == v) {
//                 count += 1;
//             }
//         }
//         return count / 2;
//     }

//     /// The number of other variables `v` interacts with.
//     size_type num_interactions(index_type v) const { return adj_[v].size(); }

//     /// Return a reference to the offset
//     bias_type& offset() { return offset_; }

//     /// Return a reference to the offset
//     const bias_type& offset() const { return offset_; }

//     /// Remove the interaction if it exists
//     bool remove_interaction(index_type u, index_type v) {
//         if (adj_[u].erase(v)) {
//             if (u != v) {
//                 adj_[v].erase(u);
//             }
//             return true;
//         } else {
//             return false;
//         }
//     }

//     /// Resize model to contain n variables.
//     void resize(index_type n) {
//         if (n < (index_type)this->num_variables()) {
//             // Clean out any of the to-be-deleted variables from the
//             // neighborhoods.
//             // This approach is better in the dense case. In the sparse case
//             // we could determine which neighborhoods need to be trimmed rather
//             // than just doing them all.
//             for (index_type v = 0; v < n; ++v) {
//                 this->adj_[v].erase(this->adj_[v].lower_bound(n), this->adj_[v].end());
//             }
//         }

//         this->linear_biases_.resize(n);
//         this->adj_.resize(n);
//     }

//     /// Scale offset, linear biases, and interactions by a factor
//     void scale(double scale_factor) {
//         // adjust offset
//         offset() *= scale_factor;

//         // adjust linear biases and quadratic interactions
//         for (size_type u = 0; u < num_variables(); u++) {
//             linear_biases_[u] *= scale_factor;

//             auto begin = adj_[u].begin();
//             auto end = adj_[u].end();
//             for (auto nit = begin; nit != end; ++nit) {
//                 (*nit).second *= scale_factor;
//             }
//         }
//     }

//     /// Set the quadratic bias for the given variables.
//     void set_quadratic(index_type u, index_type v, bias_type bias) {
//         assert(0 <= u && static_cast<size_t>(u) < this->num_variables());
//         assert(0 <= v && static_cast<size_t>(v) < this->num_variables());

//         if (u == v) {
//             switch (this->vartype(u)) {
//                 // unlike add_quadratic, setting is not really defined.
//                 case Vartype::BINARY: {
//                     throw std::domain_error(
//                             "Cannot set the quadratic bias of a binary variable with itself");
//                 }
//                 case Vartype::SPIN: {
//                     throw std::domain_error(
//                             "Cannot set the quadratic bias of a spin variable with itself");
//                 }
//                 default: {
//                     this->adj_[u][v] = bias;
//                     break;
//                 }
//             }
//         } else {
//             this->adj_[u][v] = bias;
//             this->adj_[v][u] = bias;
//         }
//     }

//     /// Exchange the contents of the quadratic model with the contents of `other`.
//     void swap(QuadraticModelBase<bias_type, index_type>& other) {
//         std::swap(this->linear_biases_, other.linear_biases_);
//         std::swap(this->adj_, other.adj_);
//         std::swap(this->offset_, other.offset_);
//     }

//     /// Swap the linear and quadratic biases between two variables.
//     void swap_variables(index_type u, index_type v) {
//         assert(0 <= u && static_cast<size_t>(u) < this->num_variables());
//         assert(0 <= v && static_cast<size_t>(v) < this->num_variables());

//         if (u == v) return;  // nothing to do!

//         // developer note, this is pretty expensive since we remove an
//         // element and then add an element to each of the neighborhoods that
//         // touch u or v. We could do it more efficiently by just swapping
//         // their values if they are both present, or only moving the elements
//         // between u and v if only one is there. But let's use the simple
//         // implementation for now.

//         // remove any references to u or v in any other neighborhoods (don't
//         // worry, we'll put them back)
//         for (auto it = this->adj_[u].begin(); it < this->adj_[u].end(); ++it) {
//             if (it->first != v) {
//                 this->adj_[it->first].erase(u);
//             }
//         }
//         for (auto it = this->adj_[v].begin(); it < this->adj_[v].end(); ++it) {
//             if (it->first != u) {
//                 this->adj_[it->first].erase(v);
//             }
//         }

//         // swap the neighborhoods of u and v
//         std::swap(this->adj_[u], this->adj_[v]);
//         std::swap(this->linear_biases_[u], this->linear_biases_[v]);

//         // now put u and v back into their neighbor's neighborhoods (according
//         // to their new indices)
//         for (auto it = this->adj_[u].begin(); it < this->adj_[u].end(); ++it) {
//             if (it->first != v) {
//                 this->adj_[it->first][u] = it->second;
//             }
//         }
//         for (auto it = this->adj_[v].begin(); it < this->adj_[v].end(); ++it) {
//             if (it->first != u) {
//                 this->adj_[it->first][v] = it->second;
//             }
//         }

//         // finally fix u and v themselves
//         if (this->adj_[u].erase(u)) {
//             auto bias = this->adj_[v][v];
//             this->adj_[u][v] = bias;
//             this->adj_[v][u] = bias;
//             this->adj_[v].erase(v);
//         }
//     }

//     virtual const bias_type& upper_bound(index_type) const = 0;

//     virtual const Vartype& vartype(index_type) const = 0;

//  protected:
//     std::vector<bias_type> linear_biases_;
//     std::vector<Neighborhood<bias_type, index_type>> adj_;

//     bias_type offset_;

//     friend class ConstQuadraticIterator<Bias, Index>;
// };


/// A Binary Quadratic Model is a quadratic polynomial over binary variables.
template <class Bias, class Index = int>
class BinaryQuadraticModel : public abc::QuadraticModelBase<Bias, Index> {
 public:
    /// The type of the base class.
    using base_type = abc::QuadraticModelBase<Bias, Index>;

    /// The first template parameter (Bias).
    using bias_type = typename base_type::bias_type;

    /// The second template parameter (Index).
    using index_type = typename base_type::index_type;

    /// Unsigned integral that can represent non-negative values.
    using size_type = typename base_type::size_type;

    /// Empty constructor. The vartype defaults to `Vartype::BINARY`.
    BinaryQuadraticModel() : BinaryQuadraticModel(Vartype::BINARY) {}

    /// Create a BQM of the given `vartype`.
    explicit BinaryQuadraticModel(Vartype vartype): base_type(), vartype_(vartype) {}

    BinaryQuadraticModel(index_type n, Vartype vartype): base_type(n), vartype_(vartype) {}

    template <class T>
    BinaryQuadraticModel(const T dense[], index_type num_variables, Vartype vartype)
            : BinaryQuadraticModel(num_variables, vartype) {
        this->add_quadratic_from_dense(dense, num_variables);
    }

    bias_type lower_bound() const { return vartype_info<bias_type>::min(this->vartype_); }

    bias_type lower_bound(index_type v) const { return vartype_info<bias_type>::min(this->vartype_); }

    bias_type upper_bound() const { return vartype_info<bias_type>::max(this->vartype_); }

    bias_type upper_bound(index_type v) const { return vartype_info<bias_type>::max(this->vartype_); }

    Vartype vartype() const { return this->vartype_; }

    Vartype vartype(index_type v) const { return this->vartype_; }

 private:
    // The vartype of the BQM
    Vartype vartype_;




 // public:
 //    /// The type of the base class.
 //    using base_type = QuadraticModelBase<Bias, Index>;

 //    /// The first template parameter (Bias).
 //    using bias_type = typename base_type::bias_type;

 //    /// The second template parameter (Index).
 //    using index_type = Index;

 //    /// Unsigned integral that can represent non-negative values.
 //    using size_type = typename base_type::size_type;

 //    /// Empty constructor. The vartype defaults to `Vartype::BINARY`.
 //    BinaryQuadraticModel() : BinaryQuadraticModel(Vartype::BINARY) {}

 //    /// Create a BQM of the given `vartype`.
 //    explicit BinaryQuadraticModel(Vartype vartype)
 //            : base_type(),
 //              vartype_(vartype),
 //              lower_bound_(vartype_info<bias_type>::default_min(vartype)),
 //              upper_bound_(vartype_info<bias_type>::default_max(vartype)) {}

 //    /// Create a BQM with `n` variables of the given `vartype`.
 //    BinaryQuadraticModel(index_type n, Vartype vartype) : BinaryQuadraticModel(vartype) {
 //        this->resize(n);
 //    }

 //    /**
 //     * Create a BQM from a dense matrix.
 //     *
 //     * `dense` must be an array of length `num_variables^2`.
 //     *
 //     * Values on the diagonal are treated differently depending on the variable
 //     * type.
 //     * If the BQM is SPIN-valued, then the values on the diagonal are
 //     * added to the offset.
 //     * If the BQM is BINARY-valued, then the values on the diagonal are added
 //     * as linear biases.
 //     *
 //     */
 //    template <class T>
 //    BinaryQuadraticModel(const T dense[], index_type num_variables, Vartype vartype)
 //            : BinaryQuadraticModel(num_variables, vartype) {
 //        add_quadratic(dense, num_variables);
 //    }

 //    /**
 //     * Add the variables, interactions and biases from another BQM.
 //     *
 //     * The size of the updated BQM will be adjusted appropriately.
 //     *
 //     * If the other BQM does not have the same vartype, the biases are adjusted
 //     * accordingly.
 //     */
 //    template <class B, class I>
 //    void add_bqm(const BinaryQuadraticModel<B, I>& bqm) {
 //        if (bqm.vartype() != this->vartype()) {
 //            // we could do this without the copy, but for now let's just do
 //            // it simply
 //            auto bqm_copy = BinaryQuadraticModel<B, I>(bqm);
 //            bqm_copy.change_vartype(vartype());
 //            this->add_bqm(bqm_copy);
 //            return;
 //        }

 //        // make sure we're big enough
 //        if (bqm.num_variables() > this->num_variables()) {
 //            this->resize(bqm.num_variables());
 //        }

 //        // offset
 //        this->offset() += bqm.offset();

 //        // linear
 //        for (size_type v = 0; v < bqm.num_variables(); ++v) {
 //            base_type::linear(v) += bqm.linear(v);
 //        }

 //        // quadratic
 //        for (size_type v = 0; v < bqm.num_variables(); ++v) {
 //            if (bqm.adj_[v].size() == 0) continue;

 //            this->adj_[v].reserve(this->adj_[v].size() + bqm.adj_[v].size());

 //            auto span = bqm.neighborhood(v);
 //            for (auto it = span.first; it != span.second; ++it) {
 //                this->adj_[v].emplace_back(it->first, it->second);
 //            }

 //            this->adj_[v].sort_and_sum();
 //        }
 //    }

 //    /**
 //     * Add the variables, interactions and biases from another BQM.
 //     *
 //     * The size of the updated BQM will be adjusted appropriately.
 //     *
 //     * If the other BQM does not have the same vartype, the biases are adjusted
 //     * accordingly.
 //     *
 //     * `mapping` must be a vector at least as long as the given BQM. It
 //     * should map the variables of `bqm` to new labels.
 //     */
 //    template <class B, class I, class T>
 //    void add_bqm(const BinaryQuadraticModel<B, I>& bqm, const std::vector<T>& mapping) {
 //        if (bqm.vartype() != this->vartype()) {
 //            // we could do this without the copy, but for now let's just do
 //            // it simply
 //            auto bqm_copy = BinaryQuadraticModel<B, I>(bqm);
 //            bqm_copy.change_vartype(vartype());
 //            this->add_bqm(bqm_copy, mapping);
 //            return;
 //        }

 //        // offset
 //        this->offset() += bqm.offset();

 //        if (bqm.num_variables() == 0)
 //            // nothing else to do, other BQM is empty
 //            return;

 //        // make sure we're big enough
 //        T max_v = *std::max_element(mapping.begin(), mapping.begin() + bqm.num_variables());
 //        if (max_v < 0) {
 //            throw std::out_of_range("contents of mapping must be non-negative");
 //        } else if ((size_type)max_v >= this->num_variables()) {
 //            this->resize(max_v + 1);
 //        }

 //        // linear
 //        for (size_type v = 0; v < bqm.num_variables(); ++v) {
 //            this->linear(mapping[v]) += bqm.linear(v);
 //        }

 //        // quadratic
 //        for (size_type v = 0; v < bqm.num_variables(); ++v) {
 //            if (bqm.adj_[v].size() == 0) continue;

 //            index_type this_v = mapping[v];

 //            this->adj_[this_v].reserve(this->adj_[this_v].size() + bqm.adj_[v].size());

 //            auto span = bqm.neighborhood(v);
 //            for (auto it = span.first; it != span.second; ++it) {
 //                this->adj_[this_v].emplace_back(mapping[it->first], it->second);
 //            }

 //            this->adj_[this_v].sort_and_sum();
 //        }
 //    }

 //    using base_type::add_quadratic;

 //    *
 //     * Construct a BQM from COO-formated iterators.
 //     *
 //     * A sparse BQM encoded in [COOrdinate] format is specified by three
 //     * arrays of (row, column, value).
 //     *
 //     * [COOrdinate]: https://w.wiki/n$L
 //     *
 //     * `row_iterator` must be a random access iterator  pointing to the
 //     * beginning of the row data. `col_iterator` must be a random access
 //     * iterator pointing to the beginning of the column data. `bias_iterator`
 //     * must be a random access iterator pointing to the beginning of the bias
 //     * data. `length` must be the number of (row, column, bias) entries.
     
 //    template <class ItRow, class ItCol, class ItBias>
 //    void add_quadratic(ItRow row_iterator, ItCol col_iterator, ItBias bias_iterator,
 //                       index_type length) {
 //        // determine the number of variables so we can resize ourself if needed
 //        if (length > 0) {
 //            index_type max_label = std::max(*std::max_element(row_iterator, row_iterator + length),
 //                                            *std::max_element(col_iterator, col_iterator + length));

 //            if ((size_t)max_label >= base_type::num_variables()) {
 //                this->resize(max_label + 1);
 //            }
 //        } else if (length < 0) {
 //            throw std::out_of_range("length must be positive");
 //        }

 //        // count the number of elements to be inserted into each
 //        std::vector<index_type> counts(base_type::num_variables(), 0);
 //        ItRow rit(row_iterator);
 //        ItCol cit(col_iterator);
 //        for (index_type i = 0; i < length; ++i, ++rit, ++cit) {
 //            if (*rit != *cit) {
 //                counts[*rit] += 1;
 //                counts[*cit] += 1;
 //            }
 //        }

 //        // reserve the neighborhoods
 //        for (size_type i = 0; i < counts.size(); ++i) {
 //            base_type::adj_[i].reserve(counts[i]);
 //        }

 //        // add the values to the neighborhoods, not worrying about order
 //        rit = row_iterator;
 //        cit = col_iterator;
 //        ItBias bit(bias_iterator);
 //        for (index_type i = 0; i < length; ++i, ++rit, ++cit, ++bit) {
 //            if (*rit == *cit) {
 //                // let add_quadratic handle this case based on vartype
 //                add_quadratic(*rit, *cit, *bit);
 //            } else {
 //                base_type::adj_[*rit].emplace_back(*cit, *bit);
 //                base_type::adj_[*cit].emplace_back(*rit, *bit);
 //            }
 //        }

 //        // finally sort and sum the neighborhoods we touched
 //        for (size_type i = 0; i < counts.size(); ++i) {
 //            if (counts[i] > 0) {
 //                base_type::adj_[i].sort_and_sum();
 //            }
 //        }
 //    }

 //    /// Add one (disconnected) variable to the BQM and return its index.
 //    index_type add_variable() {
 //        index_type vi = this->num_variables();
 //        this->resize(vi + 1);
 //        return vi;
 //    }

 //    /// Change the vartype of the binary quadratic model
 //    void change_vartype(Vartype vartype) {
 //        if (vartype == vartype_) return;  // nothing to do

 //        bias_type lin_mp, lin_offset_mp, quad_mp, quad_offset_mp, lin_quad_mp;
 //        if (vartype == Vartype::BINARY) {
 //            lin_mp = 2;
 //            lin_offset_mp = -1;
 //            quad_mp = 4;
 //            lin_quad_mp = -2;
 //            quad_offset_mp = .5;
 //        } else if (vartype == Vartype::SPIN) {
 //            lin_mp = .5;
 //            lin_offset_mp = .5;
 //            quad_mp = .25;
 //            lin_quad_mp = .25;
 //            quad_offset_mp = .125;
 //        } else {
 //            throw std::logic_error("unexpected vartype");
 //        }

 //        for (size_type ui = 0; ui < base_type::num_variables(); ++ui) {
 //            bias_type lbias = base_type::linear_biases_[ui];

 //            base_type::linear_biases_[ui] *= lin_mp;
 //            base_type::offset_ += lin_offset_mp * lbias;

 //            auto begin = base_type::adj_[ui].begin();
 //            auto end = base_type::adj_[ui].end();
 //            for (auto nit = begin; nit != end; ++nit) {
 //                bias_type qbias = (*nit).second;

 //                (*nit).second *= quad_mp;
 //                base_type::linear_biases_[ui] += lin_quad_mp * qbias;
 //                base_type::offset_ += quad_offset_mp * qbias;
 //            }
 //        }

 //        this->vartype_ = vartype;
 //        this->lower_bound_ = vartype_info<bias_type>::default_min(vartype);
 //        this->upper_bound_ = vartype_info<bias_type>::default_max(vartype);
 //    }

 //    const bias_type& lower_bound(index_type v) const { return this->lower_bound_; }

 //    /**
 //     *  Return the number of interactions in the quadratic model.
 //     *
 //     * `O(num_variables)` complexity.
 //     */
 //    size_type num_interactions() const {
 //        // we can do better than QuadraticModelBase by not needing to account
 //        // for self-loops
 //        size_type count = 0;
 //        for (auto it = this->adj_.begin(); it != this->adj_.end(); ++it) {
 //            count += it->size();
 //        }
 //        return count / 2;
 //    }

 //    /// The number of other variables `v` interacts with.
 //    size_type num_interactions(index_type v) const { return base_type::num_interactions(v); }

 //    /// Exchange the contents of the binary quadratic model with the contents of `other`.
 //    void swap(BinaryQuadraticModel<bias_type, index_type>& other) {
 //        base_type::swap(other);
 //        std::swap(this->vartype_, other.vartype_);
 //        std::swap(this->lower_bound_, other.lower_bound_);
 //        std::swap(this->upper_bound_, other.upper_bound_);
 //    }

 //    const bias_type& upper_bound(index_type v) const { return this->upper_bound_; }

 //    /// Return the vartype of the binary quadratic model.
 //    const Vartype& vartype() const { return vartype_; }

 //    /// Return the vartype of `v`.
 //    const Vartype& vartype(index_type v) const { return vartype_; }

 // private:
 //    Vartype vartype_;

 //    // hold the lower and upper bounds for the vartype
 //    bias_type lower_bound_;
 //    bias_type upper_bound_;
};

// template <class T>
// struct VarInfo {
//     Vartype vartype;
//     T lb;
//     T ub;

//     VarInfo(Vartype vartype, T lb, T ub) : vartype(vartype), lb(lb), ub(ub) {}
// };

template <class Bias, class Index = int>
class QuadraticModel : public abc::QuadraticModelBase<Bias, Index> {
 public:
    /// The type of the base class.
    using base_type = abc::QuadraticModelBase<Bias, Index>;

    /// The first template parameter (Bias).
    using bias_type = Bias;

    /// The second template parameter (Index).
    using index_type = Index;

    /// Unsigned integral that can represent non-negative values.
    using size_type = typename base_type::size_type;

    QuadraticModel(): base_type(), varinfo_() {}

    QuadraticModel(const QuadraticModel& qm) : base_type(qm), varinfo_(qm.varinfo_) {}

    QuadraticModel(const QuadraticModel&& qm) : base_type(std::move(qm)), varinfo_(std::move(qm.varinfo_)) {}

    explicit QuadraticModel(const BinaryQuadraticModel<bias_type, index_type>& bqm)
            : base_type(bqm), varinfo_(bqm.num_variables(), varinfo_type(bqm.vartype())) {}

    template<class B, class I>
    explicit QuadraticModel(const BinaryQuadraticModel<B, I>& bqm) {
        // in the future we could speed this up
        this->resize(bqm.num_variables(), bqm.vartype(), bqm.lower_bound(), bqm.upper_bound());

        for (size_type v = 0; v < bqm.num_variables(); ++v) {
            this->set_linear(v, bqm.linear(v));
        }

        for (auto qit = bqm.cbegin_quadratic(); qit != bqm.cend_quadratic(); ++qit) {
            this->add_quadratic_back(qit->u, qit->v, qit->bias);
        }

        this->set_offset(bqm.offset());
    }

    QuadraticModel& operator=(const QuadraticModel& other) {
        base_type::operator=(other);
        this->varinfo_ = other.varinfo_;
        return *this;
    }

    QuadraticModel& operator=(QuadraticModel&& other) noexcept {
        using std::swap;
        base_type::operator=(std::move(other));
        this->varinfo_ = std::move(other.varinfo_);
        return *this;
    }

    index_type add_variable(Vartype vartype) {
        return this->add_variable(vartype, vartype_info<bias_type>::default_min(vartype),
                                  vartype_info<bias_type>::default_max(vartype));
    }

    index_type add_variable(Vartype vartype, bias_type lb, bias_type ub) {
        assert(lb <= ub);

        assert(lb >= vartype_info<bias_type>::min(vartype));
        assert(ub <= vartype_info<bias_type>::max(vartype));

        assert(vartype != Vartype::BINARY || lb == 0);
        assert(vartype != Vartype::BINARY || ub == 1);

        assert(vartype != Vartype::SPIN || lb == -1);
        assert(vartype != Vartype::SPIN || ub == +1);

        index_type v = this->num_variables();

        this->varinfo_.emplace_back(vartype, lb, ub);
        base_type::add_variable();

        return v;
    }

//     void clear() {
//         this->varinfo_.clear();
//         base_type::clear();
//     }

//     /// Change the vartype of `v`, updating the biases appropriately.
//     void change_vartype(Vartype vartype, index_type v) {
//         if (vartype == this->vartype(v)) return;  // nothing to do

//         if (this->vartype(v) == Vartype::BINARY && vartype == Vartype::SPIN) {
//             // binary to spin
//             for (auto it = this->adj_[v].begin(); it != this->adj_[v].end(); ++it) {
//                 this->linear(it->first) += it->second / 2;
//                 this->adj_[it->first][v] /= 2;  // log(n)
//                 it->second /= 2;
//             }
//             this->offset() += this->linear(v) / 2;
//             this->linear(v) /= 2;

//             this->vartype(v) = Vartype::SPIN;
//             this->lower_bound(v) = -1;
//             this->upper_bound(v) = +1;
//         } else if (this->vartype(v) == Vartype::SPIN && vartype == Vartype::BINARY) {
//             // spin to binary
//             for (auto it = this->adj_[v].begin(); it != this->adj_[v].end(); ++it) {
//                 this->linear(it->first) -= it->second;
//                 this->adj_[it->first][v] *= 2;  // log(n)
//                 it->second *= 2;
//             }
//             this->offset() -= this->linear(v);
//             this->linear(v) *= 2;

//             this->vartype(v) = Vartype::BINARY;
//             this->lower_bound(v) = 0;
//             this->upper_bound(v) = 1;
//         } else if (this->vartype(v) == Vartype::BINARY && vartype == Vartype::INTEGER) {
//             // binary to integer
//             this->varinfo_[v].vartype = Vartype::INTEGER;
//         } else if (this->vartype(v) == Vartype::SPIN && vartype == Vartype::INTEGER) {
//             // spin to integer (via spin to binary)
//             this->change_vartype(Vartype::BINARY, v);
//             this->change_vartype(Vartype::INTEGER, v);
//         } else {
//             // todo: support integer to real and vice versa, need to figure
//             // out how to handle bounds in that case though
//             throw std::logic_error("invalid vartype change");
//         }
//     }

    bias_type lower_bound(index_type v) const { return this->varinfo_[v].lb; }

//     const bias_type& lower_bound(index_type v) const { return varinfo_[v].lb; }

//     constexpr bias_type max_integer() { return vartype_limits<bias_type, Vartype::INTEGER>::max(); }

//     /**
//      * Total bytes consumed by the biases, vartype info, bounds, and indices.
//      *
//      * If `capacity` is true, use the capacity of the underlying vectors rather
//      * than the size.
//      */
//     size_type nbytes(bool capacity = false) const noexcept {
//         size_type count = base_type::nbytes(capacity);
//         if (capacity) {
//             count += this->varinfo_.capacity() * sizeof(VarInfo<bias_type>);
//         } else {
//             count += this->varinfo_.size() * sizeof(VarInfo<bias_type>);
//         }
//         return count;
//     }

    // Resize the model to contain `n` variables.
    void resize(index_type n) {
        // we could do this as an assert, but let's be careful since
        // we're often calling this from python
        if (n > static_cast<index_type>(this->num_variables())) {
            throw std::logic_error(
                    "n must be smaller than the number of variables when no "
                    "`vartype` is specified");
        }
        // doesn't matter what vartype we specify since we're shrinking
        base_type::resize(n);
        this->varinfo_.erase(this->varinfo_.begin() + n, this->varinfo_.end());
    }

    /**
     * Resize the model to contain `n` variables.
     *
     * The `vartype` is used to any new variables added.
     *
     * The `vartype` must be `Vartype::BINARY` or `Vartype::SPIN`.
     */
    void resize(index_type n, Vartype vartype) {
        if (vartype == Vartype::BINARY) {
            this->resize(n, vartype, 0, 1);
        } else if (vartype == Vartype::SPIN) {
            this->resize(n, vartype, -1, +1);
        } else {
            throw std::logic_error("must provide bounds for integer vartypes when resizing");
        }
    }

    /**
     * Resize the model to contain `n` variables.
     *
     * The `vartype` is used to any new variables added.
     */
    void resize(index_type n, Vartype vartype, bias_type lb, bias_type ub) {
        assert(n > 0);

        assert(lb <= ub);

        assert(lb >= vartype_info<bias_type>::min(vartype));
        assert(ub <= vartype_info<bias_type>::max(vartype));

        assert(vartype != Vartype::BINARY || lb == 0);
        assert(vartype != Vartype::BINARY || ub == 1);

        assert(vartype != Vartype::SPIN || lb == -1);
        assert(vartype != Vartype::SPIN || ub == +1);

        this->varinfo_.resize(n, varinfo_type(vartype, lb, ub));
        base_type::resize(n);
    }

    void set_lower_bound(index_type v, bias_type lb) { this->varinfo_[v].lb = lb; }

    void set_upper_bound(index_type v, bias_type ub) { this->varinfo_[v].ub = ub; }

    bias_type upper_bound(index_type v) const { return this->varinfo_[v].ub; }

//     const bias_type& upper_bound(index_type v) const { return varinfo_[v].ub; }

//     Vartype& vartype(index_type v) { return varinfo_[v].vartype; }

    Vartype vartype(index_type v) const { return this->varinfo_[v].vartype; }

    // friend void swap(QuadraticModel& first, QuadraticModel& second) {
    //     using std::swap;
    //     throw std::logic_error("hi");
    //     base_type::swap(first, second);
    //     swap(first.varinfo_, second.varinfo_);
    // }

    // void swap(QuadraticModel<bias_type, index_type>& other) {
    //     base_type::swap(other);
    //     std::swap(this->varinfo_, other.varinfo_);
    // }

//     /// Exchange the contents of the quadratic model with the contents of `other`.
//     void swap_variables(index_type u, index_type v) {
//         base_type::swap_variables(u, v);  // also handles asserts
//         std::swap(this->varinfo_[u], this->varinfo_[v]);
//     }

 private:
    struct varinfo_type {
        Vartype vartype;
        bias_type lb;
        bias_type ub;

        varinfo_type(Vartype vartype, bias_type lb, bias_type ub) : vartype(vartype), lb(lb), ub(ub) {}

        varinfo_type(Vartype vartype): vartype(vartype) {
            this->lb = vartype_info<bias_type>::default_min(vartype);
            this->ub = vartype_info<bias_type>::default_max(vartype);
        }
    };

    std::vector<varinfo_type> varinfo_;
};

// template <class B, class N>
// std::ostream& operator<<(std::ostream& os, const BinaryQuadraticModel<B, N>& bqm) {
//     os << "BinaryQuadraticModel\n";

//     if (bqm.vartype() == Vartype::SPIN) {
//         os << "  vartype: spin\n";
//     } else if (bqm.vartype() == Vartype::BINARY) {
//         os << "  vartype: binary\n";
//     } else {
//         os << "  vartype: unkown\n";
//     }

//     os << "  offset: " << bqm.offset() << "\n";

//     os << "  linear (" << bqm.num_variables() << " variables):\n";
//     for (size_t v = 0; v < bqm.num_variables(); ++v) {
//         auto bias = bqm.linear(v);
//         if (bias) {
//             os << "    " << v << " " << bias << "\n";
//         }
//     }

//     os << "  quadratic (" << bqm.num_interactions() << " interactions):\n";
//     for (size_t u = 0; u < bqm.num_variables(); ++u) {
//         auto span = bqm.neighborhood(u);
//         for (auto nit = span.first; nit != span.second && (*nit).first < u; ++nit) {
//             os << "    " << u << " " << (*nit).first << " " << (*nit).second << "\n";
//         }
//     }

//     return os;
// }

}  // namespace dimod
