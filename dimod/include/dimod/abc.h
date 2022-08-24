// Copyright 2022 D-Wave Systems Inc.
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
#include <memory>
#include <utility>
#include <vector>

#include "dimod/vartypes.h"

namespace dimod {
namespace abc {

/**
 * Used internally by QuadraticModelBase to sparsely encode the neighborhood of
 * a variable.
 */
template <class Bias, class Index>
class Neighborhood {
 public:
    /// The first template parameter (Bias).
    using bias_type = Bias;

    /// The second template parameter (Index).
    using index_type = Index;

    /// Unsigned integral type that can represent non-negative values.
    using size_type = std::size_t;

    // /// Exactly `pair<index_type, bias_type>`.
    // using value_type = typename std::pair<index_type, bias_type>;

    /// The linear terms.
    struct value_type {
        index_type v;
        bias_type bias;

        value_type(index_type v, bias_type bias): v(v), bias(bias) {}

        bool operator<(const value_type& other) const { return this->v < other.v; }
    };

    /// A random access iterator to value_type.
    using iterator = typename std::vector<value_type>::iterator;

    /// A random access iterator to value_type`
    using const_iterator = typename std::vector<value_type>::const_iterator;

    /**
     * Return a reference to the bias associated with `v`.
     *
     * This function automatically checks whether `v` is a variable in the
     * neighborhood and throws a `std::out_of_range` exception if it is not.
     */
    const bias_type& at(index_type v) const {
        auto it = this->lower_bound(v);
        if (it != this->cend() && it->v == v) {
            // it exists
            return it->bias;
        } else {
            // it doesn't exist
            throw std::out_of_range("given variable has no interaction");
        }
    }

    value_type& back() { return this->neighborhood_.back(); }

    const value_type& back() const { return this->neighborhood_.back(); }

    /// Returns an iterator to the beginning.
    iterator begin() { return this->neighborhood_.begin(); }

    /// Returns an iterator to the end.
    iterator end() { return this->neighborhood_.end(); }

    /// Returns a const_iterator to the beginning.
    const_iterator cbegin() const { return this->neighborhood_.cbegin(); }

    /// Returns a const_iterator to the end.
    const_iterator cend() const { return this->neighborhood_.cend(); }

    /**
     * Insert a neighbor, bias pair at the end of the neighborhood.
     *
     * Note that this does not keep the neighborhood self-consistent and should
     * only be used when you know that the neighbor is greater than the current
     * last element.
     */
    void emplace_back(index_type v, bias_type bias) { this->neighborhood_.emplace_back(v, bias); }

    /// Returns whether the neighborhood is empty
    bool empty() const { return !this->size(); }

    /**
     * Erase an element from the neighborhood.
     *
     * Returns the number of element removed, either 0 or 1.
     */
    size_type erase(index_type v) {
        auto it = this->lower_bound(v);
        if (it != this->end() && it->v == v) {
            // is there to erase
            this->neighborhood_.erase(it);
            return 1;
        } else {
            return 0;
        }
    }

    /// Erase elements from the neighborhood.
    void erase(iterator first, iterator last) { this->neighborhood_.erase(first, last); }

    /// Return an iterator to the first element that does not come before `v`.
    iterator lower_bound(index_type v) {
        return std::lower_bound(this->begin(), this->end(), v, this->cmp);
    }

    /// Return an iterator to the first element that does not come before `v`.
    const_iterator lower_bound(index_type v) const {
        return std::lower_bound(this->cbegin(), this->cend(), v, this->cmp);
    }

    /**
     * Total bytes consumed by the biases and indices.
     *
     * If `capacity` is true, use the capacity of the underlying vectors rather
     * than the size.
     */
    size_type nbytes(bool capacity = false) const noexcept {
        // so there is no guaruntee that the compiler will not implement
        // pair as pointers or whatever, but this seems like a reasonable
        // assumption.
        if (capacity) {
            return this->neighborhood_.capacity() * sizeof(value_type);
        } else {
            return this->neighborhood_.size() * sizeof(value_type);
        }
    }

    /**
     * Return the bias at neighbor `v` or the default value.
     *
     * Return the bias of `v` if `v` is in the neighborhood, otherwise return
     * the `value` provided without inserting `v`.
     */
    bias_type get(index_type v, bias_type value = 0) const {
        auto it = this->lower_bound(v);

        if (it != this->cend() && it->v == v) {
            // it exists
            return it->bias;
        } else {
            // it doesn't exist
            return value;
        }
    }

    /// Request that the neighborhood capacity be at least enough to contain `n`
    /// elements.
    void reserve(index_type n) { this->neighborhood_.reserve(n); }

    /// Return the size of the neighborhood.
    size_type size() const { return this->neighborhood_.size(); }

    /// Sort the neighborhood and sum the biases of duplicate variables.
    void sort_and_sum() {
        if (!std::is_sorted(this->begin(), this->end())) {
            std::sort(this->begin(), this->end());
        }

        // now remove any duplicates, summing the biases of duplicates
        size_type i = 0;
        size_type j = 1;

        // walk quickly through the neighborhood until we find a duplicate
        while (j < this->neighborhood_.size() &&
               this->neighborhood_[i].v != this->neighborhood_[j].v) {
            ++i;
            ++j;
        }

        // if we found one, move into de-duplication
        while (j < this->neighborhood_.size()) {
            if (this->neighborhood_[i].v == this->neighborhood_[j].v) {
                this->neighborhood_[i].bias += this->neighborhood_[j].bias;
                ++j;
            } else {
                ++i;
                this->neighborhood_[i] = this->neighborhood_[j];
                ++j;
            }
        }

        // finally resize to contain only the unique values
        this->neighborhood_.resize(i + 1);
    }

    /**
     * Access the bias of `v`.
     *
     * If `v` is in the neighborhood, the function returns a reference to
     * its bias. If `v` is not in the neighborhood, it is inserted and a
     * reference is returned to its bias.
     */
    bias_type& operator[](index_type v) {
        auto it = this->lower_bound(v);
        if (it == this->end() || it->v != v) {
            // it doesn't exist so insert
            it = this->neighborhood_.emplace(it, v, 0);
        }
        return it->bias;
    }

 protected:
    std::vector<value_type> neighborhood_;

    static inline bool cmp(value_type ub, index_type v) { return ub.v < v; }
};


template <class Bias, class Index>
class QuadraticModelBase {
 public:
    /// The first template parameter (Bias).
    using bias_type = Bias;

    /// The second template parameter (Index).
    using index_type = Index;

    /// Unsigned integral type that can represent non-negative values.
    using size_type = std::size_t;

    class const_quadratic_iterator {
     public:
        struct value_type {
            index_type u;
            index_type v;
            bias_type bias;

            explicit value_type(index_type u) : u(u), v(-1), bias(NAN) {}
        };
        using pointer = const value_type*;
        using reference = const value_type&;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::forward_iterator_tag;

        const_quadratic_iterator() : adj_ptr_(nullptr), term_(0) {}

        const_quadratic_iterator(const std::vector<Neighborhood<bias_type, index_type>>* adj_ptr,
                                 index_type u)
                : adj_ptr_(adj_ptr), term_(u), vi_(0) {
            if (adj_ptr_ != nullptr) {
                // advance through the neighborhoods until we find one on the
                // lower triangle
                for (index_type& u = this->term_.u;
                     static_cast<size_type>(u) < this->adj_ptr_->size(); ++u) {
                    auto& neighborhood = (*adj_ptr_)[u];

                    if (neighborhood.size() && neighborhood.cbegin()->v <= u) {
                        // we found one
                        this->term_.v = neighborhood.cbegin()->v;
                        this->term_.bias = neighborhood.cbegin()->bias;
                        return;
                    }
                }
            }
        }

        const reference operator*() const { return this->term_; }

        const pointer operator->() const { return &(this->term_); }

        const_quadratic_iterator& operator++() {
            index_type& vi = this->vi_;

            ++vi;  // advance to the next

            for (index_type& u = this->term_.u; static_cast<size_type>(u) < this->adj_ptr_->size();
                 ++u) {
                auto& neighborhood = (*adj_ptr_)[u];

                auto it = neighborhood.cbegin() + vi;

                if (it != neighborhood.cend() && it->v <= u) {
                    // we found one
                    this->term_.v = it->v;
                    this->term_.bias = it->bias;
                    return *this;
                }

                vi = 0;
            }

            return *this;
        }

        const_quadratic_iterator operator++(int) {
            const_quadratic_iterator tmp(*this);
            ++(*this);
            return tmp;
        }

        friend bool operator==(const const_quadratic_iterator& a,
                               const const_quadratic_iterator& b) {
            return (a.adj_ptr_ == nullptr && b.adj_ptr_ == nullptr) ||
                   (a.adj_ptr_ == b.adj_ptr_ && a.term_.u == b.term_.u && a.vi_ == b.vi_);
        }

        friend bool operator!=(const const_quadratic_iterator& a,
                               const const_quadratic_iterator& b) {
            return !(a == b);
        }

     private:
        // note that unlike QuandraticModelBase, we use a regular pointer
        // because the QuadraticModelBase owns the adjacency structure
        const std::vector<Neighborhood<bias_type, index_type>>* adj_ptr_;

        value_type term_;  // the current term

        index_type vi_;  // term_.v's location in the neighborhood of term_.u
    };

    friend class const_quadratic_iterator;

    QuadraticModelBase() : linear_biases_(), adj_ptr_(), offset_(0) {}

    /// Copy constructor
    QuadraticModelBase(const QuadraticModelBase& other)
            : linear_biases_(other.linear_biases_), adj_ptr_(), offset_(other.offset_) {
        // need to handle the adj if present
        if (!other.is_linear()) {
            this->adj_ptr_ = std::unique_ptr<std::vector<Neighborhood<bias_type, index_type>>>(
                    new std::vector<Neighborhood<bias_type, index_type>>(*other.adj_ptr_));
        }
    }

    using const_neighborhood_iterator = typename Neighborhood<bias_type, index_type>::const_iterator;

    /// Move constructor
    QuadraticModelBase(QuadraticModelBase&& other) noexcept { *this = std::move(other); }

    /// Copy assignment operator
    QuadraticModelBase& operator=(const QuadraticModelBase& other) {
        if (this != &other) {
            this->linear_biases_ = other.linear_biases_;  // copy
            if (other.has_adj()) {
                this->adj_ptr_ = std::unique_ptr<std::vector<Neighborhood<bias_type, index_type>>>(
                        new std::vector<Neighborhood<bias_type, index_type>>(*other.adj_ptr_));
            } else {
                this->adj_ptr_.reset(nullptr);
            }

            this->offset_ = other.offset_;
        }
        return *this;
    }

    /// Move assignment operator
    QuadraticModelBase& operator=(QuadraticModelBase&& other) noexcept {
        if (this != &other) {
            this->linear_biases_ = std::move(other.linear_biases_);
            this->adj_ptr_ = std::move(other.adj_ptr_);
            this->offset_ = other.offset_;
        }
        return *this;
    }

    virtual ~QuadraticModelBase() {}

    void add_linear(index_type v, bias_type bias) {
        assert(v >= 0 && static_cast<size_type>(v) < this->num_variables());
        this->linear_biases_[v] += bias;
    }

    void add_quadratic(index_type u, index_type v, bias_type bias) {
        assert(0 <= u && static_cast<size_t>(u) < this->num_variables());
        assert(0 <= v && static_cast<size_t>(v) < this->num_variables());

        this->enforce_adj();

        if (u == v) {
            switch (this->vartype(u)) {
                case Vartype::BINARY: {
                    // 1*1 == 1 and 0*0 == 0 so this is linear
                    this->linear_biases_[u] += bias;
                    break;
                }
                case Vartype::SPIN: {
                    // -1*-1 == +1*+1 == 1 so this is a constant offset
                    this->offset_ += bias;
                    break;
                }
                default: {
                    // self-loop
                    (*adj_ptr_)[u][v] += bias;
                    break;
                }
            }
        } else {
            (*adj_ptr_)[u][v] += bias;
            (*adj_ptr_)[v][u] += bias;
        }
    }

    void add_quadratic(std::initializer_list<index_type> row, std::initializer_list<index_type> col,
                       std::initializer_list<bias_type> biases) {
        auto rit = row.begin();
        auto cit = col.begin();
        auto bit = biases.begin();
        for (; rit < row.end() && cit < col.end() && bit < biases.end(); ++rit, ++cit, ++bit) {
            this->add_quadratic(*rit, *cit, *bit);
        }
    }

    /**
     * Add quadratic bias for the given variables at the end of eachother's neighborhoods.
     *
     * # Parameters
     * - `u` - a variable.
     * - `v` - a variable.
     * - `bias` - the quadratic bias associated with `u` and `v`.
     *
     * # Exceptions
     * When `u` is less than the largest neighbor in `v`'s neighborhood,
     * `v` is less than the largest neighbor in `u`'s neighborhood, or either
     * `u` or `v` is greater than ``num_variables()`` then the behavior of
     * this method is undefined.
     */
    void add_quadratic_back(index_type u, index_type v, bias_type bias) {
        assert(0 <= u && static_cast<size_t>(u) <= this->num_variables());
        assert(0 <= v && static_cast<size_t>(v) <= this->num_variables());

        this->enforce_adj();

        // check the condition for adding at the back
        assert((*adj_ptr_)[v].empty() || (*adj_ptr_)[v].back().v <= u);
        assert((*adj_ptr_)[u].empty() || (*adj_ptr_)[u].back().v <= v);

        if (u == v) {
            switch (this->vartype(u)) {
                case Vartype::BINARY: {
                    // 1*1 == 1 and 0*0 == 0 so this is linear
                    this->add_linear(u, bias);
                    break;
                }
                case Vartype::SPIN: {
                    // -1*-1 == +1*+1 == 1 so this is a constant offset
                    this->offset_ += bias;
                    break;
                }
                default: {
                    // self-loop
                    (*adj_ptr_)[u].emplace_back(v, bias);
                    break;
                }
            }
        } else {
            (*adj_ptr_)[u].emplace_back(v, bias);
            (*adj_ptr_)[v].emplace_back(u, bias);
        }
    }

    /*
     * Add quadratic biases from a dense matrix.
     *
     * `dense` must be an array of length `num_variables^2`.
     *
     * Values on the diagonal are treated differently depending on the variable
     * type.
     *
     * # Exceptions
     * The behavior of this method is undefined when the model has fewer than
     * `num_variables` variables.
     */
    template<class T>
    void add_quadratic_from_dense(const T dense[], index_type num_variables) {
        assert(0 <= num_variables);
        assert(static_cast<size_type>(num_variables) <= this->num_variables());

        this->enforce_adj();

        if (this->is_linear()) {
            for (index_type u = 0; u < num_variables; ++u) {
                // diagonal
                this->add_quadratic_back(u, u, dense[u * (num_variables + 1)]);

                // off-diagonal
                for (index_type v = u + 1; v < num_variables; ++v) {
                    bias_type qbias = dense[u * num_variables + v] + dense[v * num_variables + u];

                    if (qbias) {
                        this->add_quadratic_back(u, v, qbias);
                    }
                }
            }
        } else {
            // we cannot rely on the ordering
            for (index_type u = 0; u < num_variables; ++u) {
                // diagonal
                this->add_quadratic(u, u, dense[u * (num_variables + 1)]);

                // off-diagonal
                for (index_type v = u + 1; v < num_variables; ++v) {
                    bias_type qbias = dense[u * num_variables + v] + dense[v * num_variables + u];

                    if (qbias) {
                        this->add_quadratic(u, v, qbias);
                    }
                }
            }
        }
    }

    const_neighborhood_iterator cbegin_neighborhood(index_type v) const {
        if (has_adj()) {
            return (*adj_ptr_)[v].begin();
        } else {
            // I am a bit suspicious that this works, but it seems to
            return const_neighborhood_iterator(nullptr);
        }
    }

    const_neighborhood_iterator cend_neighborhood(index_type v) const {
        if (has_adj()) {
            return (*adj_ptr_)[v].end();
        } else {
            // I am a bit suspicious that this works, but it seems to
            return const_neighborhood_iterator(nullptr);
        }
    }

    const_quadratic_iterator cbegin_quadratic() const {
        return const_quadratic_iterator(this->adj_ptr_.get(), 0);
    }

    const_quadratic_iterator cend_quadratic() const {
        return const_quadratic_iterator(this->adj_ptr_.get(), this->num_variables());
    }

    /**
     * Return the energy of the given sample.
     *
     * The `sample_start` must be random access iterator pointing to the
     * beginning of the sample.
     *
     * The behavior of this function is undefined when the sample is not
     * `num_variables()` long.
     */
    template <class Iter>  // todo: allow different return types
    bias_type energy(Iter sample_start) const {
        bias_type en = this->offset();

        if (has_adj()) {
            for (index_type u = 0; static_cast<size_type>(u) < this->num_variables(); ++u) {
                auto u_val = *(sample_start + u);


                en += u_val * this->linear(u);

                for (auto& term : (*adj_ptr_)[u]) {
                    if (term.v > u) break;
                    en += term.bias * u_val * *(sample_start + term.v);
                }
            }
        } else {
            throw std::logic_error("todo 353");
        }

        return en;
    }

    /**
     * Return the energy of the given sample.
     *
     * The `sample` must be a vector containing the sample.
     *
     * The behavior of this function is undefined when the sample is not
     * `num_variables()` long.
     */
    template <class T>
    bias_type energy(const std::vector<T>& sample) const {
        // todo: check length?
        return energy(sample.cbegin());
    }

    template <class B, class I>
    bool is_equal(const QuadraticModelBase<B, I>& other) const {
        // easy checks first
        if (this->num_variables() != other.num_variables() ||
            this->num_interactions() != other.num_interactions() ||
            this->offset_ != other.offset_ ||
            !std::equal(this->linear_biases_.begin(), this->linear_biases_.end(),
                        other.linear_biases_.begin())) {
            return false;
        }

        // check the vartype
        for (size_type v = 0; v < this->num_variables(); ++v) {
            if (this->vartype(v) != other.vartype(v)) {
                return false;
            }
        }

        // check quadratic. We already checked the number of interactions so
        // if one is present we can assume it's empty
        if (this->has_adj() && other.has_adj()) {
            // check that the neighborhoods are equal
            throw std::logic_error("todo 142");
        }

        return true;
    }

    bool is_linear() const {
        if (this->has_adj()) {
            for (const auto& n : *adj_ptr_) {
                if (n.size()) return false;
            }
        }
        return true;
    }

    bias_type linear(index_type v) const { return this->linear_biases_[v]; }

    // bias_type linear_at(index_type v)

    virtual bias_type lower_bound(index_type) const = 0;

    size_type num_interactions() const {
        size_type count = 0;
        if (this->has_adj()) {
            index_type v = 0;
            for (auto& n : *(this->adj_ptr_)) {
                count += n.size();

                // account for self-loops
                auto lb = n.lower_bound(v);
                if (lb != n.cend() && lb->v == v) {
                    count += 1;
                }

                ++v;
            }
        }
        return count / 2;
    }

    size_type num_interactions(index_type v) const {
        if (this->has_adj()) {
            return (*adj_ptr_)[v].size();
        } else {
            return 0;
        }
    }

    size_type num_variables() const { return this->linear_biases_.size(); }

    bias_type offset() const { return this->offset_; }

    bias_type quadratic(index_type u, index_type v) const {
        if (!adj_ptr_) {
            return 0;
        }
        return (*adj_ptr_)[u].get(v);
    }

    bias_type quadratic_at(index_type u, index_type v) const {
        if (!adj_ptr_) {
            throw std::out_of_range("given variables have no interaction");
        }
        return (*adj_ptr_)[u].at(v);
    }

    bool remove_interaction(index_type u, index_type v) {
        if (!this->has_adj()) {
            return false;
        }
        if ((*adj_ptr_)[u].erase(v)) {
            if (u != v) {
                (*adj_ptr_)[v].erase(u);
            }
            return true;
        } else {
            return false;
        }
    }

    /// Resize model to contain n variables.
    void resize(index_type n) {
        assert(n >= 0);
        this->linear_biases_.resize(n);

        if (has_adj()) {
            if (static_cast<size_type>(n) < this->num_variables()) {
                // Clean out any of the to-be-deleted variables from the
                // neighborhoods.
                // This approach is better in the dense case. In the sparse case
                // we could determine which neighborhoods need to be trimmed rather
                // than just doing them all.
                throw std::logic_error("not done");
                // for (index_type v = 0; v < n; ++v) {
                //     this->adj_[v].erase(this->adj_[v].lower_bound(n), this->adj_[v].end());
                // }
            }
            this->adj_ptr_->resize(n);
        }
    }

    void scale(bias_type scalar) {
        this->offset_ *= scalar;

        // linear biases
        for (bias_type& bias : this->linear_biases_) {
            bias *= scalar;
        }

        if (!this->has_adj()) {
            // we're done
            return;
        }

        // quadratic biases
        for (auto& n : (*adj_ptr_)) {
            for (auto& term : n) {
                term.bias *= scalar;
            }
        }
    }

    void set_linear(index_type v, bias_type bias) {
        assert(v >= 0 && static_cast<size_type>(v) < this->num_variables());
        this->linear_biases_[v] = bias;
    }

    void set_linear(index_type v, std::initializer_list<bias_type> biases) {
        assert(v >= 0 && static_cast<size_type>(v) < this->num_variables());
        assert(biases.size() + v <= this->num_variables());
        for (const bias_type& b : biases) {
            this->linear_biases_[v] = b;
            ++v;
        }
    }

    void set_offset(bias_type offset) { this->offset_ = offset; }

    void set_quadratic(index_type u, index_type v, bias_type bias) {
        assert(0 <= u && static_cast<size_t>(u) < this->num_variables());
        assert(0 <= v && static_cast<size_t>(v) < this->num_variables());

        this->enforce_adj();

        if (u == v) {
            switch (this->vartype(u)) {
                // unlike add_quadratic, setting is not really defined.
                case Vartype::BINARY: {
                    throw std::domain_error(
                            "Cannot set the quadratic bias of a binary variable with itself");
                }
                case Vartype::SPIN: {
                    throw std::domain_error(
                            "Cannot set the quadratic bias of a spin variable with itself");
                }
                default: {
                    (*adj_ptr_)[u][v] = bias;
                    break;
                }
            }
        } else {
            (*adj_ptr_)[u][v] = bias;
            (*adj_ptr_)[v][u] = bias;
        }
    }

    virtual bias_type upper_bound(index_type) const = 0;

    virtual Vartype vartype(index_type) const = 0;

 protected:
    explicit QuadraticModelBase(std::vector<bias_type>&& biases)
            : linear_biases_(biases), adj_ptr_(), offset_(0) {}

    explicit QuadraticModelBase(index_type n) : linear_biases_(n), adj_ptr_(), offset_(0) {}

    index_type add_variable() { return this->add_variables(1); }

    index_type add_variables(index_type n) {
        assert(n >= 0);
        index_type size = this->num_variables();

        this->linear_biases_.resize(size + n);
        if (this->has_adj()) {
            this->adj_ptr_->resize(size + n);
        }

        return size;
    }

 private:
    std::vector<bias_type> linear_biases_;

    std::unique_ptr<std::vector<Neighborhood<bias_type, index_type>>> adj_ptr_;

    bias_type offset_;

    inline void enforce_adj() {
        if (!this->adj_ptr_) {
            this->adj_ptr_ = std::unique_ptr<std::vector<Neighborhood<bias_type, index_type>>>(
                    new std::vector<Neighborhood<bias_type, index_type>>(this->num_variables()));
        }
    }

    inline bool has_adj() const { return static_cast<bool>(this->adj_ptr_); }
};

}  // namespace abc
}  // namespace dimod
