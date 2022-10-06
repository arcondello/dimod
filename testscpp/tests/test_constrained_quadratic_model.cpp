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
#include "dimod/quadratic_model.h"

namespace dimod {

SCENARIO("ConstrainedQuadraticModel  tests") {
    GIVEN("an empty CQM") {
        auto cqm = ConstrainedQuadraticModel<double>();

        THEN("some basic properties can be discovered") {
            CHECK(cqm.num_variables() == 0);
            CHECK(cqm.num_constraints() == 0);
        }

        THEN("the cqm can be copied") {
            auto cqm2 = cqm;

            THEN("some basic properties can be discovered") {
                CHECK(cqm2.num_variables() == 0);
                CHECK(cqm2.num_constraints() == 0);
            }
        }

        WHEN("we add variables") {
            cqm.add_variables(Vartype::INTEGER, 5);
            cqm.add_variables(Vartype::REAL, 3, -10, 10);
            cqm.add_variable(Vartype::SPIN);
            cqm.add_variable(Vartype::BINARY, 0, 1);

            THEN("we can read them off appropriately") {
                REQUIRE(cqm.num_variables() == 10);

                for (int i = 0; i < 5; ++i) CHECK(cqm.vartype(i) == Vartype::INTEGER);
                for (int i = 5; i < 8; ++i) CHECK(cqm.vartype(i) == Vartype::REAL);
                for (int i = 5; i < 8; ++i) CHECK(cqm.lower_bound(i) == -10);
                for (int i = 5; i < 8; ++i) CHECK(cqm.upper_bound(i) == +10);
                CHECK(cqm.vartype(8) == Vartype::SPIN);
                CHECK(cqm.vartype(9) == Vartype::BINARY);
            }

            AND_WHEN("we try to set the linear biases in the objective") {
                cqm.objective.set_linear(0, 1.5);

                THEN("it is reflected in the model") {
                    REQUIRE(cqm.objective.num_variables() == 1);
                    CHECK(cqm.objective.linear(0) == 1.5);
                }
            }
        }

        AND_GIVEN("a quadratic model") {
            auto qm = QuadraticModel<double>();
            auto u = qm.add_variable(Vartype::INTEGER, -5, 5);
            auto v = qm.add_variable(Vartype::BINARY);
            qm.set_linear(u, 1);
            qm.set_linear(v, 2);
            qm.set_quadratic(u, v, 1.5);
            qm.set_offset(10);

            WHEN("we set the objective via the set_objective() method") {
                cqm.set_objective(qm);

                THEN("the objective updates appropriately") {
                    REQUIRE(cqm.objective.num_variables() == 2);
                    REQUIRE(cqm.objective.num_interactions() == 1);
                    CHECK(cqm.objective.linear(0) == 1);
                    CHECK(cqm.objective.linear(1) == 2);
                    CHECK(cqm.objective.quadratic(0, 1) == 1.5);
                    CHECK(cqm.objective.offset() == 10);
                    CHECK(cqm.lower_bound(0) == -5);
                    CHECK(cqm.upper_bound(0) == 5);
                    CHECK(cqm.vartype(0) == Vartype::INTEGER);
                    CHECK(cqm.vartype(1) == Vartype::BINARY);
                }
            }

            WHEN("we set the objective via the set_objective() with a relabel") {
                cqm.add_variable(Vartype::BINARY);
                cqm.add_variable(Vartype::INTEGER, -5, 5);
                cqm.set_objective(qm, std::vector<int>{1, 0});

                THEN("the objective updates appropriately") {
                    REQUIRE(cqm.objective.num_variables() == 2);
                    REQUIRE(cqm.objective.num_interactions() == 1);
                    CHECK(cqm.objective.linear(1) == 1);
                    CHECK(cqm.objective.linear(0) == 2);
                    CHECK(cqm.objective.quadratic(0, 1) == 1.5);
                    CHECK(cqm.objective.offset() == 10);
                    CHECK(cqm.lower_bound(1) == -5);
                    CHECK(cqm.upper_bound(1) == 5);
                    CHECK(cqm.vartype(1) == Vartype::INTEGER);
                    CHECK(cqm.vartype(0) == Vartype::BINARY);
                }

                AND_WHEN("that objective is overwritten") {
                    auto qm2 = QuadraticModel<double>();
                    qm2.add_variable(Vartype::BINARY);
                    qm2.set_linear(0, 10);
                    cqm.set_objective(qm2, std::vector<int>{0});

                    THEN("everything updates as expected") {
                        REQUIRE(cqm.objective.num_variables() == 1);
                        REQUIRE(cqm.objective.num_interactions() == 0);
                        CHECK(cqm.objective.linear(1) == 0);
                        CHECK(cqm.objective.linear(0) == 10);
                        CHECK(cqm.objective.quadratic(0, 1) == 0);
                        CHECK(cqm.objective.offset() == 0);
                        CHECK(cqm.lower_bound(1) == -5);
                        CHECK(cqm.upper_bound(1) == 5);
                        CHECK(cqm.vartype(1) == Vartype::INTEGER);
                        CHECK(cqm.vartype(0) == Vartype::BINARY);
                    }
                }
            }
        }

        WHEN("we add an empty constraint") {
            auto c0 = cqm.add_constraint();

            THEN("that constraint is well formed") {
                CHECK(cqm.constraint_ref(c0).num_variables() == 0);
                CHECK(cqm.constraint_ref(c0).num_interactions() == 0);
                CHECK(!cqm.constraint_ref(c0).is_soft());
                CHECK(cqm.constraint_ref(c0).rhs() == 0);
            }

            AND_WHEN("we manipulate that constraint") {
                auto u = cqm.add_variable(Vartype::BINARY);
                auto v = cqm.add_variable(Vartype::INTEGER, -5, 5);
                cqm.constraint_ref(c0).add_linear(u, 1.5);
                cqm.constraint_ref(c0).add_linear(v, 12.5);
            }
        }

        WHEN("the objective is set via a QM") {
            auto qm = QuadraticModel<float>();
            qm.add_variable(Vartype::SPIN);
            qm.add_variable(Vartype::REAL, -5, +17);
            qm.set_linear(0, {1.5, -5});
            qm.add_quadratic(0, 1, -2);
            qm.set_offset(5);

            cqm.set_objective(qm);

            REQUIRE(cqm.num_variables() == 2);
            CHECK(cqm.objective.num_variables() == 2);
            REQUIRE(cqm.objective.num_interactions() == 1);
            CHECK(cqm.objective.linear(0) == 1.5);
            CHECK(cqm.objective.linear(1) == -5);
            CHECK(cqm.objective.quadratic_at(0, 1) == -2);
            CHECK(cqm.objective.offset() == 5);
            CHECK(cqm.objective.lower_bound(0) == qm.lower_bound(0));
            CHECK(cqm.objective.upper_bound(0) == qm.upper_bound(0));
            CHECK(cqm.objective.lower_bound(1) == qm.lower_bound(1));
            CHECK(cqm.objective.upper_bound(1) == qm.upper_bound(1));
        }

        WHEN("10 variables and two empty constraints are added") {
            cqm.add_constraints(2);
            cqm.add_variable(Vartype::INTEGER);
            cqm.add_variable(Vartype::INTEGER);

            THEN("quadratic biases can be added") {
                REQUIRE(cqm.num_constraints() == 2);

                cqm.constraint_ref(0).add_quadratic(0, 1, 1.5);

                CHECK(cqm.constraint_ref(0).num_interactions() == 1);
                CHECK(cqm.constraint_ref(0).num_interactions(0) == 1);
                CHECK(cqm.constraint_ref(0).num_interactions(1) == 1);
                CHECK(cqm.constraint_ref(0).quadratic(0, 1) == 1.5);
                CHECK(cqm.constraint_ref(0).quadratic(1, 0) == 1.5);
            }
        }
    }

    GIVEN("a CQM with an objective and a constraint") {
        auto cqm = ConstrainedQuadraticModel<double>();

        // auto qm = QuadraticModel<double>();
        cqm.add_variables(Vartype::INTEGER, 5, 0, 5);
        cqm.add_variables(Vartype::BINARY, 5);
        // qm.set_linear(0, {0, 1, -2, 3, -4, 5, -6, 7, -8, 9});
        // qm.add_quadratic(0, 1, 1);
        // qm.add_quadratic(1, 2, 2);
        // qm.set_offset(5);
        // cqm.set_objective(qm);

        // cqm.add_constraints(1, Sense::EQ);
        // cqm.constraint_ref(0).set_linear(2, 2);
        // cqm.constraint_ref(0).set_linear(5, -5);
        // cqm.constraint_ref(0).add_quadratic(2, 4, 8);
        // cqm.constraint_ref(0).set_offset(4);

        // todo: test Move, copy, constructors and operators
    }

    GIVEN("a CQM with 30 variables") {
        auto cqm = ConstrainedQuadraticModel<double>();

        cqm.add_variables(Vartype::BINARY, 10);
        cqm.add_variables(Vartype::INTEGER, 10, -5, 5);
        cqm.add_variables(Vartype::REAL, 10, -100, 105);

        WHEN("we add a constraint over only a few variables") {
            auto c0 = cqm.add_constraint();

            auto& constraint = cqm.constraint_ref(c0);

            constraint.add_linear(5, 15);
            constraint.add_linear(7, 105);
            constraint.add_linear(3, 5);
            constraint.add_quadratic(5, 7, 56);
            constraint.add_quadratic(3, 7, 134);

            THEN("we can iterate over the quadratic interactions") {
                auto it = constraint.cbegin_quadratic();

                CHECK(it->u == 7);
                CHECK(it->v == 5);
                CHECK(it->bias == 56);
                it++;
                CHECK(it->u == 3);
                CHECK(it->v == 7);
                CHECK(it->bias == 134);
                it++;
                CHECK(it == constraint.cend_quadratic());
            }

            THEN("we can iterate over the neighborhoods") {
                auto it = constraint.cbegin_neighborhood(7);

                CHECK(it->v == 5);
                CHECK(it->bias == 56);
                it++;
                CHECK(it->v == 3);
                CHECK(it->bias == 134);
                it++;
                CHECK(it == constraint.cend_neighborhood(7));
            }
        }
    }

    GIVEN("a CQM with a contraint") {
        auto cqm = ConstrainedQuadraticModel<double>();
        cqm.add_variables(Vartype::BINARY, 3);

        auto c0 = cqm.add_constraint();
        cqm.constraint_ref(c0).set_linear(2, 1);
        cqm.constraint_ref(c0).set_sense(Sense::GE);
        cqm.constraint_ref(c0).set_rhs(1);

        std::vector<int> sample{0, 0, 1};
        CHECK(cqm.constraint_ref(c0).energy(sample.begin()) == 1);
    }

    GIVEN("A CQM with several constraints") {
        auto cqm = ConstrainedQuadraticModel<double>();
        cqm.add_variables(Vartype::BINARY, 7);

        auto c0 = cqm.add_constraint();
        auto c1 = cqm.add_constraint();
        auto c2 = cqm.add_constraint();

        cqm.constraint_ref(c0).add_linear(0, 1.5);
        cqm.constraint_ref(c0).add_linear(1, 2.5);
        cqm.constraint_ref(c0).add_linear(2, 3.5);

        cqm.constraint_ref(c1).add_linear(2, 4.5);
        cqm.constraint_ref(c1).add_linear(3, 5.5);
        cqm.constraint_ref(c1).add_linear(4, 6.5);

        cqm.constraint_ref(c2).add_linear(5, 8.5);
        cqm.constraint_ref(c2).add_linear(4, 7.5);

        THEN("we can read off the values as expected") {
            const auto& const2 = cqm.constraint_ref(c2);
            REQUIRE(const2.num_variables() == 2);
            CHECK(const2.linear(4) == 7.5);
            CHECK(const2.linear(5) == 8.5);
            CHECK(const2.variables() == std::vector<int>{5, 4});
        }

        WHEN("constraint c1 is removed") {
            cqm.remove_constraint(c1);

            THEN("all variables are preserved, but one constraint is removed") {
                CHECK(cqm.num_variables() == 7);
                CHECK(cqm.num_constraints() == 2);
                CHECK(cqm.constraint_ref(1).linear(4) == 7.5);
            }
        }

        WHEN("a variable is removed") {
            cqm.remove_variable(3);

            THEN("everything is updated appropriately") {
                REQUIRE(cqm.num_variables() == 6);

                const auto& const0 = cqm.constraint_ref(c0);
                REQUIRE(const0.num_variables() == 3);
                CHECK(const0.linear(0) == 1.5);
                CHECK(const0.linear(1) == 2.5);
                CHECK(const0.linear(2) == 3.5);
                CHECK(const0.variables() == std::vector<int>{0, 1, 2});

                const auto& const1 = cqm.constraint_ref(c1);
                REQUIRE(const1.num_variables() == 2);
                CHECK(const1.linear(2) == 4.5);
                CHECK(const1.linear(3) == 6.5);                       // reindexed
                CHECK(const1.variables() == std::vector<int>{2, 3});  // partly reindexed

                const auto& const2 = cqm.constraint_ref(c2);
                REQUIRE(const2.num_variables() == 2);
                CHECK(const2.linear(3) == 7.5);                       // reindexed
                CHECK(const2.linear(4) == 8.5);                       // reindexed
                CHECK(const2.variables() == std::vector<int>{4, 3});  // reindexed
            }
        }
    }

    GIVEN("A CQM with two constraints") {
        auto cqm = ConstrainedQuadraticModel<double>();
        auto x = cqm.add_variable(Vartype::BINARY);
        auto y = cqm.add_variable(Vartype::BINARY);
        auto i = cqm.add_variable(Vartype::INTEGER);
        auto j = cqm.add_variable(Vartype::INTEGER);
        // auto z = cqm.add_variable(Vartype::BINARY);

        cqm.objective.set_linear(x, 1);
        cqm.objective.set_linear(y, 2);
        cqm.objective.set_linear(i, 3);
        cqm.objective.set_linear(j, 4);

        auto& const0 = cqm.constraint_ref(cqm.add_constraint());
        const0.set_linear(i, 3);
        const0.set_quadratic(x, j, 2);
        const0.set_quadratic(i, j, 5);

        CHECK(cqm.objective.num_variables() == 4);
        CHECK(const0.num_variables() == 3);

        WHEN("we substitute a variable with a 0 multiplier") {
            cqm.substitute_variable(x, 0, 0);

            THEN("the biases are updated in the objective and constraint") {
                REQUIRE(cqm.objective.num_variables() == 4);
                CHECK(cqm.objective.linear(x) == 0);
                CHECK(cqm.objective.linear(y) == 2);
                CHECK(cqm.objective.linear(i) == 3);
                CHECK(cqm.objective.linear(j) == 4);

                auto& const0 = cqm.constraint_ref(0);
                REQUIRE(const0.num_variables() == 3);
                REQUIRE(const0.num_interactions() == 2);
                CHECK(const0.linear(i) == 3);
                CHECK(const0.quadratic_at(x, j) == 0);
                CHECK(const0.quadratic_at(i, j) == 5);

            }
        }

        // WHEN("we fix a variable") {
        //     cqm.fix_variable(x, 0);

        //     THEN("everything is updated correctly") {
        //         REQUIRE(cqm.num_variables() == 3);

        //         REQUIRE(const0.num_variables() == 2);
        //         REQUIRE(const0.linear(i-1) == 3);
        //     }
        // }
    }

    GIVEN("A constraint with one-hot constraints") {
        auto cqm = ConstrainedQuadraticModel<double>();
        cqm.add_variables(Vartype::BINARY, 10);
        auto c0 = cqm.add_linear_constraint({0, 1, 2, 3, 4}, {1, 1, 1, 1, 1}, Sense::EQ, 1);
        auto c1 = cqm.add_linear_constraint({5, 6, 7, 8, 9}, {2, 2, 2, 2, 2}, Sense::EQ, 2);

        THEN("the constraints can be tests for one-hotness") {
            CHECK(cqm.constraint_ref(c0).is_onehot());
            CHECK(cqm.constraint_ref(c1).is_onehot());
            CHECK(cqm.constraint_ref(c0).is_disjoint(cqm.constraint_ref(c1)));
        }

        WHEN("we change a linear bias") {
            cqm.constraint_ref(c0).set_linear(0, 1.5);

            THEN("it's no longer one-hot") { CHECK(!cqm.constraint_ref(c0).is_onehot()); }
        }

        WHEN("we add an overlapping variable") {
            cqm.constraint_ref(c0).set_linear(5, 1);
            THEN("they are no longer disjoint") {
                CHECK(cqm.constraint_ref(c0).is_onehot());
                CHECK(cqm.constraint_ref(c1).is_onehot());
                CHECK(!cqm.constraint_ref(c0).is_disjoint(cqm.constraint_ref(c1)));
            }
        }
    }
}
}  // namespace dimod
