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

#include "../Catch2/single_include/catch2/catch.hpp"
#include "dimod/constrained.h"


namespace dimod {

SCENARIO("constrained quadratic models") {
    GIVEN("an empty CQM") {
        auto cqm = ConstrainedQuadraticModel<double>();

        THEN("some basic properties can be discovered") {
            CHECK(cqm.num_variables() == 0);
            // CHECK(cqm.num_constraints() == 0);
        }

        WHEN("the objective is set via a QM") {
            auto qm = QuadraticModel<float>();
            qm.add_variable(Vartype::SPIN);
            qm.add_variable(Vartype::REAL, -5, +17);
            qm.linear(0) = 1.5;
            qm.linear(1) = -5;
            qm.add_quadratic(0, 1, -2);
            qm.offset() = 5;

            cqm.set_objective(qm);

            REQUIRE(cqm.num_variables() == 2);
            CHECK(cqm.objective().num_variables() == 2);
            REQUIRE(cqm.objective().num_interactions() == 1);
            CHECK(cqm.objective().linear(0) == 1.5);
            CHECK(cqm.objective().linear(1) == -5);
            CHECK(cqm.objective().quadratic_at(0, 1) == -2);
            CHECK(cqm.objective().offset() == 5);
            CHECK(cqm.objective().lower_bound(0) == qm.lower_bound(0));
            CHECK(cqm.objective().upper_bound(0) == qm.upper_bound(0));
            CHECK(cqm.objective().lower_bound(1) == qm.lower_bound(1));
            CHECK(cqm.objective().upper_bound(1) == qm.upper_bound(1));
        }

        WHEN("one constraint is copied") {
            cqm.add_variables(10, Vartype::BINARY);
            CHECK(cqm.num_variables() == 10);

            std::vector<int> variables {2, 4, 7};
            std::vector<float> biases {20, 40, 70};

            cqm.add_linear_constraint(variables, biases, Sense::Le, 5);

            REQUIRE(cqm.num_constraints() == 1);
            // CHECK(cqm.constraint(0).linear(0) == 0);
        }
    }
}


}  // namespace dimod
