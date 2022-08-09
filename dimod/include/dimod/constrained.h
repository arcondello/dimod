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

#include <unordered_map>
#include <utility>
#include <vector>

#include "dimod/abc.h"
#include "dimod/quadratic_model.h"

namespace dimod {

enum Sense { Le, Ge, Eq };

template <class Bias, class Index>
class ConstrainedQuadraticModel;

template <class Bias, class Index>
class Constraint : public abc::QuadraticModelBase<Bias, Index> {
 public:
    /// The first template parameter (Bias).
    using bias_type = Bias;

    /// The second template parameter (Index).
    using index_type = Index;

    /// Unsigned integral type that can represent non-negative values.
    using size_type = std::size_t;

    using base_type = abc::QuadraticModelBase<bias_type, index_type>;

    using parent_type = ConstrainedQuadraticModel<bias_type, index_type>;

 private:
    const parent_type* parent_;

    std::vector<index_type> variables_;
    Sense sense_;
    bias_type rhs_;

    std::unordered_map<index_type, index_type> indices_;

 public:
    Constraint(const parent_type* parent, std::vector<index_type>&& variables,
               std::vector<bias_type>&& biases, Sense sense, bias_type rhs)
            : base_type(std::move(biases)),
              parent_(parent),
              variables_(std::move(variables)),
              sense_(sense),
              rhs_(rhs) {
        index_type i = 0;
        for (const index_type& v : this->variables_) {
            this->indices_[v] = i++;
        }
    }

    bias_type linear(index_type v) const {
        assert(v >= 0 && static_cast<size_type>(v) < this->parent_->num_variables());

        auto it = this->indices_.find(v);
        if (it == this->indices_.end()) {
            return 0;
        }
        return base_type::linear(it->second);
    }

    bias_type lower_bound(index_type v) const { return this->parent_->lower_bound(v); }

    bias_type upper_bound(index_type v) const { return this->parent_->upper_bound(v); }

    Vartype vartype(index_type v) const { return this->parent_->vartype(v); }
};

template <class Bias, class Index = int>
class ConstrainedQuadraticModel {
 public:
    /// The first template parameter (Bias).
    using bias_type = Bias;

    /// The second template parameter (Index).
    using index_type = Index;

    /// Unsigned integral type that can represent non-negative values.
    using size_type = std::size_t;

    // using linear_type = LinearTerm<bias_type, index_type>;

    using constraint_type = Constraint<bias_type, index_type>;

 private:
    // The objective holds the information about the variables, including
    // the vartype, bounds, number etc.
    QuadraticModel<bias_type, index_type> objective_;

    std::vector<constraint_type> constraints_;

 public:
    index_type add_linear_constraint(std::vector<index_type>&& variables,
                                     std::vector<bias_type>&& biases, Sense sense, bias_type rhs) {
        this->constraints_.emplace_back(this, std::move(variables), std::move(biases), sense, rhs);
        return this->constraints_.size() - 1;
    }

    template <class B, class I>
    index_type add_linear_constraint(const std::vector<I>& variables, const std::vector<B>& biases,
                                     Sense sense, bias_type rhs) {
        // make some vectors we can move
        std::vector<index_type> v(variables.begin(), variables.end());
        std::vector<bias_type> b(biases.begin(), biases.end());
        return this->add_linear_constraint(std::move(v), std::move(b), sense, rhs);
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

    constraint_type& constraint(index_type c) { return this->constraints_[c]; }

    const constraint_type& constraint(index_type c) const { return this->constraints_[c]; }

    /// Return the lower bound for variable `v`.
    bias_type lower_bound(index_type v) const { return this->objective_.lower_bound(v); }

    size_type num_constraints() const { return this->constraints_.size(); }

    /// Return the number of variables in the model.
    size_type num_variables() const { return this->objective_.num_variables(); }

    const QuadraticModel<bias_type, index_type>& objective() const { return this->objective_; }

    /// Return the upper bound for variable `v`.
    bias_type upper_bound(index_type v) const { return this->objective_.upper_bound(v); }

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

    /// Return the vartype of `v`.
    Vartype vartype(index_type v) const { return this->objective_.vartype(v); }
};

}  // namespace dimod
