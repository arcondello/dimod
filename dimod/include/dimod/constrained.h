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
#include <unordered_map>
#include <utility>
#include <vector>

#include "dimod/quadratic_model.h"

namespace dimod {

enum Sense { Le, Ge, Eq };

template <class Bias, class Index>
struct LinearTerm {
    Index v;
    Bias bias;

    LinearTerm(Index v, Bias bias) : v(v), bias(bias) {}
};

template <class Bias, class Index>
class ConstraintBase {
 public:
    /// The first template parameter (Bias).
    using bias_type = Bias;

    /// The second template parameter (Index).
    using index_type = Index;

    /// Unsigned integral type that can represent non-negative values.
    using size_type = std::size_t;

 private:
    Sense sense_;
    bias_type rhs_;

 public:
    ConstraintBase(Sense sense, bias_type rhs) : sense_(sense), rhs_(rhs) {}

    virtual ~ConstraintBase() {}
};

template <class Bias, class Index>
class LinearConstraint : public ConstraintBase<Bias, Index> {
    /// The type of the base class.
    using base_type = ConstraintBase<Bias, Index>;

    /// The first template parameter (Bias).
    using bias_type = Bias;

    /// The second template parameter (Index).
    using index_type = Index;

    /// Unsigned integral type that can represent non-negative values.
    using size_type = std::size_t;

    using linear_type = LinearTerm<bias_type, index_type>;

 private:
    std::vector<linear_type> terms_;

 public:
    LinearConstraint(std::vector<linear_type>&& terms, Sense sense, bias_type rhs)
            : base_type(sense, rhs), terms_(std::move(terms)) {}
};

// template <class Bias, class Index>
// class QuadraticConstraint: public ConstraintBase<Bias, Index> {
//     /// The type of the base class.
//     using base_type = ConstraintBase<Bias, Index>;

//     /// The first template parameter (Bias).
//     using bias_type = Bias;

//     /// The second template parameter (Index).
//     using index_type = Index;

//     /// Unsigned integral type that can represent non-negative values.
//     using size_type = std::size_t;

//     using linear_type = LinearTerm<bias_type, index_type>;

//  private:
//     std::vector<linear_type> terms_;
//     std::vector<Neighborhood<bias_type, index_type>> adj_;

//  public:
//     QuadraticConstraint
// };

template <class Bias, class Index = int>
class ConstrainedQuadraticModel {
 public:
    /// The first template parameter (Bias).
    using bias_type = Bias;

    /// The second template parameter (Index).
    using index_type = Index;

    /// Unsigned integral type that can represent non-negative values.
    using size_type = std::size_t;

    using linear_type = LinearTerm<bias_type, index_type>;

 private:
    // The objective holds the information about the variables, including
    // the vartype, bounds, number etc.
    QuadraticModel<bias_type, index_type> objective_;

    std::vector<std::unique_ptr<ConstraintBase<bias_type, index_type>>> constraints_;

 public:
    template <class B, class I>
    index_type add_linear_constraint(const std::vector<I>& variables, const std::vector<B>& biases,
                                     Sense sense, bias_type rhs) {
        std::vector<linear_type> terms;

        auto vit = variables.begin();
        auto bit = biases.begin();
        for (; vit != variables.end() && bit != biases.end(); ++vit, ++bit) {
            assert(*vit >= 0 && static_cast<size_type>(*vit) < this->num_variables());
            terms.emplace_back(*vit, *bit);
        }

        return this->add_linear_constraint(std::move(terms), sense, rhs);
    }

    index_type add_linear_constraint(std::vector<linear_type>&& terms, Sense sense, bias_type rhs) {
        this->constraints_.emplace_back(std::unique_ptr<ConstraintBase<bias_type, index_type>>(
                new LinearConstraint<bias_type, index_type>(std::move(terms), sense, rhs)));
        return this->constraints_.size() - 1;
    }

    index_type add_variable(Vartype vartype) { return this->objective_.add_variable(vartype); }

    index_type add_variable(Vartype vartype, bias_type lb, bias_type ub) {
        return this->objective_.add_variable(vartype, lb, ub);
    }

    void add_variables(index_type n, Vartype vartype) {
        for (index_type i = 0; i < n; ++i) {
            this->objective_.add_variable(vartype);
        }
    }

    void add_variables(index_type n, Vartype vartype, bias_type lb, bias_type ub) {
        for (index_type i = 0; i < n; ++i) {
            this->objective_.add_variable(vartype, lb, ub);
        }
    }

    /// Return the lower bound for variable `v`.
    bias_type& lower_bound(index_type v) { return this->objective_.lower_bound(v); }

    /// Return the lower bound for variable `v`.
    const bias_type& lower_bound(index_type v) const { return this->objective_.lower_bound(v); }

    size_type num_constraints() const { return this->constraints_.size(); }

    /// Return the number of variables in the model.
    size_type num_variables() const { return this->objective_.num_variables(); }

    const QuadraticModel<bias_type, index_type>& objective() const { return this->objective_; }

    /// Return the upper bound for variable `v`.
    bias_type& upper_bound(index_type v) { return this->objective_.upper_bound(v); }

    /// Return the upper bound for variable `v`.
    const bias_type& upper_bound(index_type v) const { return this->objective_.upper_bound(v); }

    template <class B, class I>
    void set_objective(const QuadraticModelBase<B, I>& objective) {
        if (!this->num_variables()) {
            // the objective is empty, so we can just add, easy peasy
            for (size_type i = 0; i < objective.num_variables(); ++i) {
                this->objective_.add_variable(objective.vartype(i), objective.lower_bound(i),
                                              objective.upper_bound(i));
                this->objective_.linear(i) = objective.linear(i);
            }

            for (auto qit = objective.cbegin_quadratic(); qit != objective.cend_quadratic();
                 ++qit) {
                this->objective_.add_quadratic_back(qit->u, qit->v, qit->bias);
            }

            this->objective_.offset() = objective.offset();

            return;
        }

        throw std::logic_error("not implemented");
    }

    // TODO: objectives from vectors?

    /// Return the vartype of `v`.
    const Vartype& vartype(index_type v) const { return this->objective_.vartype(v); }
};

}  // namespace dimod
