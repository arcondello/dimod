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

#include <memory>
#include <vector>

#include "dimod/quadratic_model.h"  // Neighborhood

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

 private:
    std::vector<bias_type> linear_biases_;

    std::unique_ptr<std::vector<Neighborhood<bias_type, index_type>>> adj_ptr_;

    bias_type offset_;

 public:
    QuadraticModelBase() : linear_biases_(), adj_ptr_(), offset_(0) {}

    // QuadraticModelBase()

    // virtual ~QuadraticModelBase() {}

 protected:
    explicit QuadraticModelBase(std::vector<bias_type>&& biases)
            : linear_biases_(biases), adj_ptr_(), offset_(0) {}

    index_type add_variable() {
        this->linear_biases_.push_back(0);
        if (this->has_adj()) {
            this->adj_ptr_->resize(this->adj_ptr_->size() + 1);
        }
        return this->linear_biases_.size() - 1;
    }

 public:
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

    bias_type quadratic(index_type u, index_type v) const {
        if (!adj_ptr_) {
            return 0;
        }
        return (*adj_ptr_)[u].get(v);
    }

    bias_type quadratic_at(index_type u, index_type v) const {
        if (!adj_ptr_) {
            return 0;
        }
        return (*adj_ptr_)[u].at(v);
    }

    // void set_linear(index_type v, bias_type bias) { this->linear_biases_[v] = bias; }

    virtual bias_type upper_bound(index_type) const = 0;

    virtual Vartype vartype(index_type) const = 0;

 private:
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
