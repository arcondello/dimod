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
#include <memory>
#include <utility>
#include <vector>

#include "dimod/neighborhood.h"
#include "dimod/vartypes.h"

namespace dimod {
namespace abc {

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

        const_quadratic_iterator() : adj_ptr_(nullptr) {}

        const_quadratic_iterator(const std::vector<Neighborhood<bias_type, index_type>>* adj_ptr,
                                 index_type u)
                : adj_ptr_(adj_ptr), term_(u), vi_(0) {
            if (adj_ptr_ != nullptr) {
                // advance through the neighborhoods until we find one on the
                // lower triangle
                for (index_type& u = this->term_.u;
                     static_cast<size_type>(u) < this->adj_ptr_->size(); ++u) {
                    auto& neighborhood = (*adj_ptr_)[u];

                    if (neighborhood.size() && neighborhood.cbegin()->first <= u) {
                        // we found one
                        this->term_.v = neighborhood.cbegin()->first;
                        this->term_.bias = neighborhood.cbegin()->second;
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

                if (it != neighborhood.cend() && it->first <= u) {
                    // we found one
                    this->term_.v = it->first;
                    this->term_.bias = it->second;
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
    QuadraticModelBase(const QuadraticModelBase& qm)
            : linear_biases_(qm.linear_biases_), adj_ptr_(), offset_(qm.offset_) {
        // need to handle the adj if present
        if (!qm.is_linear()) {
            throw std::logic_error("todo 56");
        }
    }

    /// Move constructor
    QuadraticModelBase(QuadraticModelBase&& other) noexcept { *this = std::move(other); }

    /// Copy assignment operator
    QuadraticModelBase& operator=(const QuadraticModelBase& other) {
        if (this != &other) {
            this->linear_biases_ = other.linear_biases_;  // copy
            if (other.has_adj()) {
                this->enforce_adj();
                (*adj_ptr_) = (*other.adj_ptr_);  // copy
            } else {
                this->adj_ptr_.reset(nullptr);
            }

            this->offset_ = other.offset_;
        }
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

    const_quadratic_iterator cbegin_quadratic() const {
        return const_quadratic_iterator(this->adj_ptr_.get(), 0);
    }

    const_quadratic_iterator cend_quadratic() const {
        return const_quadratic_iterator(this->adj_ptr_.get(), this->num_variables());
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
                if (lb != n.cend() && lb->first == v) {
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
                term.second *= scalar;
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

    friend void swap(QuadraticModelBase& first, QuadraticModelBase& second) {
        using std::swap;
        swap(first.linear_biases_, second.linear_biases_);
        swap(first.adj_ptr_, second.adj_ptr_);
        swap(first.offset_, second.offset_);
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
