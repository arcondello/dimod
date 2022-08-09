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

 public:
    // void add_linear(index_type v, bias_type bias) { this->linear_biases_[v] += bias; }

    bias_type linear(index_type v) const { return this->linear_biases_[v]; }

    // bias_type linear_at(index_type v)

    virtual bias_type lower_bound(index_type) const = 0;

    // void set_linear(index_type v, bias_type bias) { this->linear_biases_[v] = bias; }

    virtual bias_type upper_bound(index_type) const = 0;

    virtual Vartype vartype(index_type) const = 0;
};

}  // namespace abc
}  // namespace dimod
