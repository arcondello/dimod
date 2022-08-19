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

#include "dimod/constrained.h"

namespace dimod {
namespace presolve {

enum Status { FEASIBLE, INFEASIBLE, UNBOUNDED };

template <class Bias, class Index>
class PreSolveBase;

template <class Bias, class Index>
class PostSolver {
 public:
    /// The first template parameter (Bias).
    using bias_type = Bias;

    /// The second template parameter (Index).
    using index_type = Index;

    /// Unsigned integral type that can represent non-negative values.
    using size_type = std::size_t;
};

template <class Bias, class Index = int>
class PreSolver {
 public:
    /// The first template parameter (Bias).
    using bias_type = Bias;

    /// The second template parameter (Index).
    using index_type = Index;

    /// Unsigned integral type that can represent non-negative values.
    using size_type = std::size_t;

    using model_type = ConstrainedQuadraticModel<bias_type, index_type>;

    using postsolver_type = PostSolver<bias_type, index_type>;

 private:
    model_type model_;

    std::vector<std::unique_ptr<PreSolveBase<bias_type, index_type>>> presolvers_;

    postsolver_type postsolver_;

 public:
    template <class B, class I>
    explicit PreSolver(const ConstrainedQuadraticModel<B, I>& cqm) : model_(cqm) {}

    explicit PreSolver(model_type&& cqm) : model_(std::move(cqm)) {}

    template <class T>
    void add_presolver() {
        this->presolvers_.emplace_back(
                std::unique_ptr<PreSolveBase<bias_type, index_type>>(new T()));
    }

    void apply() {
        for (auto&& presolver : this->presolvers_) {
            presolver->apply(this);
        }
    }

    /// Remove any variables where the upper and lower bounds are the same
    size_type remove_fixed_variables() {
        return 0;  // todo
    }

    // Remove any constraints with 0 or 1 variables.
    size_type remove_trivial_constraints() {
        size_type original_size = this->model_.num_constraints();

        size_type i = 0;
        while (i < this->model_.num_constraints()) {
            auto& c = this->model_.constraints[i];

            if (c.num_variables() == 0) {
                // remove after checking feasibity
                throw std::logic_error("hello 78");
            } else if (c.num_variables() == 1) {
                index_type v = c.variable(0);

                // a * x + b ◯ c
                bias_type a = c.linear(v);
                if (!a) {
                    throw std::logic_error("hello 98");
                }
                bias_type rhs = (c.rhs() - c.offset()) / a;

                switch (c.sense()) {
                    case Sense::EQ:
                        this->model_.set_lower_bound(v, std::max(rhs, this->model_.lower_bound(v)));
                        this->model_.set_upper_bound(v, std::min(rhs, this->model_.upper_bound(v)));
                        break;
                    case Sense::LE:
                        this->model_.set_upper_bound(v, std::min(rhs, this->model_.upper_bound(v)));
                        break;
                    case Sense::GE:
                        this->model_.set_lower_bound(v, std::max(rhs, this->model_.lower_bound(v)));
                        break;
                }

                this->model_.remove_constraint(i);
            } else {
                // advance
                ++i;
            }
        }

        return original_size - this->model_.num_constraints();
    }

    const model_type& model() const { return this->model_; }

    const postsolver_type& postsolver() const { return this->postsolver_; }
};

template <class Bias, class Index>
class PreSolveBase {
 public:
    /// The first template parameter (Bias).
    using bias_type = Bias;

    /// The second template parameter (Index).
    using index_type = Index;

    /// Unsigned integral type that can represent non-negative values.
    using size_type = std::size_t;

    virtual void apply(PreSolver<bias_type, index_type>* presolver) = 0;
};

namespace techniques {

template <class Bias, class Index = int>
class TrivialPresolver : public PreSolveBase<Bias, Index> {
 public:
    /// The first template parameter (Bias).
    using bias_type = Bias;

    /// The second template parameter (Index).
    using index_type = Index;

    /// Unsigned integral type that can represent non-negative values.
    using size_type = std::size_t;

    void apply(PreSolver<bias_type, index_type>* presolver) {
        // any constraints of the form a*x <= b can just be removed
        while (presolver->remove_trivial_constraints() || presolver->remove_fixed_variables()) {}
    }
};
}  // namespace techniques
}  // namespace presolve
}  // namespace dimod
