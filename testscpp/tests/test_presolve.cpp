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

#include "catch2/catch.hpp"
#include "dimod/constrained_quadratic_model.h"
#include "dimod/presolve.h"

namespace dimod {

SCENARIO("constrained quadratic models can be presolved") {
    GIVEN("a cqm with some trivial issues") {
        auto cqm = ConstrainedQuadraticModel<double>();

        // v0 has default bounds, but a constraint restricting it
        auto v0 = cqm.add_variable(Vartype::INTEGER);
        auto c0 = cqm.add_linear_constraint({v0}, {1}, Sense::LE, 5);  // v0 <= 5

        // v1 has bounds that implicitly fix it
        auto v1 = cqm.add_variable(Vartype::INTEGER, 5, 5);

        // v2 has default bounds, but an equality constraint that fixes it
        auto v2 = cqm.add_variable(Vartype::INTEGER);
        auto c2 = cqm.add_linear_constraint({v2}, {1}, Sense::EQ, 7);

        // v3 has default bounds, but two inequality constraints that fix it
        auto v3 = cqm.add_variable(Vartype::INTEGER);
        auto c3a = cqm.add_linear_constraint({v3}, {1}, Sense::LE, 5.5);
        auto c3b = cqm.add_linear_constraint({v3}, {1}, Sense::GE, 4.5);

        WHEN("passed into a presolver") {
            auto presolver = presolve::PreSolver<double>(std::move(cqm));

            THEN("it has moved correctly") {
                const auto& newcqm = presolver.model();

                REQUIRE(newcqm.num_variables() == 4);
                REQUIRE(newcqm.num_constraints() == 4);

                CHECK(newcqm.vartype(v0) == Vartype::INTEGER);
                CHECK(newcqm.lower_bound(v0) == 0);
                CHECK(newcqm.upper_bound(v0) ==
                      vartype_limits<double, Vartype::INTEGER>::default_max());
                CHECK(newcqm.constraint_ref(c0).linear(v0) == 1);
                CHECK(newcqm.constraint_ref(c0).sense() == Sense::LE);
                CHECK(newcqm.constraint_ref(c0).rhs() == 5);

                CHECK(newcqm.vartype(v1) == Vartype::INTEGER);
                CHECK(newcqm.lower_bound(v1) == 5);
                CHECK(newcqm.upper_bound(v1) == 5);

                CHECK(newcqm.vartype(v2) == Vartype::INTEGER);
                CHECK(newcqm.lower_bound(v2) == 0);
                CHECK(newcqm.upper_bound(v2) ==
                      vartype_limits<double, Vartype::INTEGER>::default_max());
                CHECK(newcqm.constraint_ref(c2).linear(v2) == 1);
                CHECK(newcqm.constraint_ref(c2).sense() == Sense::EQ);
                CHECK(newcqm.constraint_ref(c2).rhs() == 7);

                CHECK(newcqm.vartype(v3) == Vartype::INTEGER);
                CHECK(newcqm.lower_bound(v3) == 0);
                CHECK(newcqm.upper_bound(v3) ==
                      vartype_limits<double, Vartype::INTEGER>::default_max());
                CHECK(newcqm.constraint_ref(c3a).linear(v3) == 1);
                CHECK(newcqm.constraint_ref(c3a).sense() == Sense::LE);
                CHECK(newcqm.constraint_ref(c3a).rhs() == 5.5);
                CHECK(newcqm.constraint_ref(c3b).linear(v3) == 1);
                CHECK(newcqm.constraint_ref(c3b).sense() == Sense::GE);
                CHECK(newcqm.constraint_ref(c3b).rhs() == 4.5);
            }

            AND_WHEN("the default presolving is applied") {
                presolver.load_default_presolvers();
                presolver.apply();

                THEN("most of the constraints/variables are removed") {
                    CHECK(presolver.model().num_constraints() == 0);
                    CHECK(presolver.model().num_variables() == 1);
                }

                AND_WHEN("we then undo the transformation") {
                    auto original = presolver.postsolver().apply(std::vector<int>{3});
                    CHECK(original == std::vector<int>{3, 5, 7, 5});
                }
            }
        }
    }

    GIVEN("a cqm with some spin variables") {
        auto cqm = ConstrainedQuadraticModel<double>();
        cqm.add_variables(Vartype::SPIN, 5);
        for (size_t v = 0; v < 5; ++v) {
            for (size_t u = v + 1; u < 5; ++u) {
                cqm.objective.set_quadratic(u, v, 1);
            }
        }

        WHEN("the default presolving is applied") {
            auto presolver = presolve::PreSolver<double>(std::move(cqm));
            presolver.load_default_presolvers();
            presolver.apply();

            THEN("most of the constraints/variables are removed") {
                CHECK(presolver.model().num_constraints() == 0);
                CHECK(presolver.model().num_variables() == 5);

                for (size_t v = 0; v < 5; ++v) {
                    CHECK(presolver.model().vartype(v) == Vartype::BINARY);
                }

                AND_WHEN("we then undo the transformation") {
                    auto original = presolver.postsolver().apply(std::vector<int>{0, 1, 0, 1, 0});
                    CHECK(original == std::vector<int>{-1, +1, -1, +1, -1});
                }
            }
        }
    }

    GIVEN("a CQM with self-loops") {
        auto cqm = ConstrainedQuadraticModel<double>();
        cqm.add_variables(Vartype::INTEGER, 5);

        cqm.objective.set_quadratic(0, 0, 1.5);
        cqm.objective.set_quadratic(3, 3, 3.5);

        cqm.add_constraint();
        cqm.constraint_ref(0).add_quadratic(4, 4, 5);
        cqm.constraint_ref(0).add_quadratic(3, 3, 6);

        WHEN("presolving is applied") {
            auto presolver = presolve::PreSolver<double>(std::move(cqm));
            presolver.load_default_presolvers();
            presolver.apply();

            THEN("the self-loops are removed") {
                auto& model = presolver.model();

                REQUIRE(model.num_variables() == 8);
                REQUIRE(model.num_constraints() == 1);
                CHECK(model.objective.num_interactions() == 2);
                CHECK(model.objective.quadratic(0, 5) == 1.5);
                CHECK(model.objective.quadratic(3, 6) == 3.5);

                CHECK(model.constraint_ref(0).num_interactions() == 2);
                CHECK(model.constraint_ref(0).quadratic(4, 7) == 5);
                CHECK(model.constraint_ref(0).quadratic(3, 6) == 6);
            }

                AND_WHEN("we then undo the transformation") {
                    auto original = presolver.postsolver().apply(std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8});
                    CHECK(original == std::vector<int>{1, 2, 3, 4, 5});
                }
        }
    }

    GIVEN("a CQM with 5 binary variables") {
        auto cqm = ConstrainedQuadraticModel<double>();
        cqm.add_variables(Vartype::BINARY, 5);

        WHEN("we add a discrete constraint that has an offset but not a rhs") {
            auto& c1 = cqm.constraint_ref(cqm.add_constraint());
            c1.set_linear(0, 1);
            c1.set_linear(1, 1);
            c1.set_offset(-1);
            c1.mark_discrete();

            AND_WHEN("we presolve is applied") {
                auto presolver = presolve::PreSolver<double>(std::move(cqm));
                presolver.load_default_presolvers();
                presolver.apply();

                THEN("the constraint is still marked as discrete") {
                    REQUIRE(presolver.model().num_constraints() == 1);
                    CHECK(presolver.model().constraint_ref(0).marked_discrete());
                }
            }
        }
        WHEN("we add a discrete constraint that has an offset and a rhs") {
            auto& c1 = cqm.constraint_ref(cqm.add_constraint());
            c1.set_linear(0, 1);
            c1.set_linear(1, 1);
            c1.set_rhs(1);
            c1.set_offset(-1);
            c1.mark_discrete();

            AND_WHEN("we presolve is applied") {
                auto presolver = presolve::PreSolver<double>(std::move(cqm));
                presolver.load_default_presolvers();
                presolver.apply();

                THEN("the constraint is still marked as discrete") {
                    REQUIRE(presolver.model().num_constraints() == 1);
                    CHECK(!presolver.model().constraint_ref(0).marked_discrete());
                }
            }
        }

        WHEN("we add two overlapping discrete constraints and one non-overlapping") {
            auto& c1 = cqm.constraint_ref(cqm.add_constraint());
            c1.set_linear(0, 1);
            c1.set_linear(1, 1);
            c1.set_rhs(1);
            c1.mark_discrete();
            auto& c2 = cqm.constraint_ref(cqm.add_constraint());
            c2.set_linear(2, 1);
            c2.set_linear(1, 1);
            c2.set_rhs(1);
            c2.mark_discrete();
            auto& c3 = cqm.constraint_ref(cqm.add_constraint());
            c3.set_linear(3, 1);
            c3.set_linear(4, 1);
            c3.set_rhs(1);
            c3.mark_discrete();

            AND_WHEN("we presolve is applied") {
                auto presolver = presolve::PreSolver<double>(std::move(cqm));
                presolver.load_default_presolvers();
                presolver.apply();

                THEN("two will still be marked as discrete") {
                    REQUIRE(presolver.model().num_constraints() == 3);
                    CHECK(presolver.model().constraint_ref(0).marked_discrete() !=
                          presolver.model().constraint_ref(1).marked_discrete());
                    CHECK(presolver.model().constraint_ref(2).marked_discrete());
                }
            }
        }
    }

    GIVEN("a CQM with one integer variable and two bound constraitns") {
        auto cqm = ConstrainedQuadraticModel<double>();
        cqm.add_variable(Vartype::INTEGER);  // unbounded

        auto& constraint = cqm.constraint_ref(cqm.add_linear_constraint({0}, {1}, Sense::LE, 5));
        constraint.set_weight(5);  // make soft

        CHECK(constraint.is_soft());

        cqm.add_linear_constraint({0}, {1}, Sense::LE, 10);  // hard constraint

        WHEN("we presolve is applied") {
            auto presolver = presolve::PreSolver<double>(std::move(cqm));
            presolver.load_default_presolvers();
            presolver.apply();

            THEN("the hard constraint is removed") {
                REQUIRE(presolver.model().num_constraints() == 1);
                CHECK(presolver.model().upper_bound(0) == 10);
            }
        }
    }
}

}  // namespace dimod
