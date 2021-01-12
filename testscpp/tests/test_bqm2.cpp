// Copyright 2020 D-Wave Systems Inc.
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

#include <cinttypes>
#include <stdexcept>
#include <vector>

#include "../Catch2/single_include/catch2/catch.hpp"
#include "dimod/adjvectorbqm2.h"

namespace dimod {

// make it easier to swap in later
template <class B>
using AdjVectorBQM = AdjVectorBQM2<B>;

SCENARIO("bqms can be constructed from COO iterators") {
    GIVEN("some COO arrays") {
        int irow[4] = {0, 2, 0, 1};
        int icol[4] = {0, 2, 1, 2};
        float bias[4] = {.5, -2, 2, -3};

        WHEN("a bqm is constructed with 0 length") {
            auto bqm = AdjVectorBQM<float>(&irow[0], &icol[0], &bias[0], 0);

            THEN("the bqm is empty") { REQUIRE(bqm.num_variables() == 0); }
        }

        WHEN("a bqm is constructed") {
            auto bqm = AdjVectorBQM<float>(&irow[0], &icol[0], &bias[0], 4);

            THEN("is takes its values form the arrays") {
                REQUIRE(bqm.num_variables() == 3);
                REQUIRE(bqm.linear(0) == .5);
                REQUIRE(bqm.linear(1) == 0.);
                REQUIRE(bqm.linear(2) == -2);

                REQUIRE(bqm.num_interactions() == 2);
                REQUIRE(bqm.quadratic(0, 1) == 2);
                REQUIRE(bqm.quadratic(2, 1) == -3);
                REQUIRE_THROWS_AS(bqm.quadratic_at(0, 2), std::out_of_range);
            }
        }
    }

    GIVEN("some COO arrays with duplicates") {
        int irow[6] = {0, 2, 0, 1, 0, 0};
        int icol[6] = {0, 2, 1, 2, 1, 0};
        float bias[6] = {.5, -2, 2, -3, 4, 1};

        WHEN("a bqm is constructed") {
            auto bqm = AdjVectorBQM<float>(&irow[0], &icol[0], &bias[0], 6);

            THEN("it combines duplicate values") {
                REQUIRE(bqm.num_variables() == 3);
                REQUIRE(bqm.linear(0) == 1.5);
                REQUIRE(bqm.linear(1) == 0.);
                REQUIRE(bqm.linear(2) == -2);

                REQUIRE(bqm.num_interactions() == 2);
                REQUIRE(bqm.quadratic(0, 1) == 6);
                REQUIRE(bqm.quadratic(2, 1) == -3);
                REQUIRE_THROWS_AS(bqm.quadratic_at(0, 2), std::out_of_range);
            }
        }
    }

    GIVEN("some COO arrays with multiple duplicates") {
        int irow[4] = {0, 1, 0, 1};
        int icol[4] = {1, 2, 1, 0};
        float bias[4] = {-1, 1, -2, -3};

        WHEN("a bqm is constructed") {
            auto bqm = AdjVectorBQM<float>(&irow[0], &icol[0], &bias[0], 4);

            THEN("it combines duplicate values") {
                REQUIRE(bqm.num_variables() == 3);
                REQUIRE(bqm.linear(0) == 0);
                REQUIRE(bqm.linear(1) == 0);
                REQUIRE(bqm.linear(2) == 0);

                REQUIRE(bqm.num_interactions() == 2);
                REQUIRE(bqm.quadratic(0, 1) == -6);
                REQUIRE(bqm.quadratic(1, 0) == -6);
                REQUIRE(bqm.quadratic(2, 1) == 1);
                REQUIRE(bqm.quadratic(1, 2) == 1);
                REQUIRE_THROWS_AS(bqm.quadratic_at(0, 2), std::out_of_range);
                REQUIRE_THROWS_AS(bqm.quadratic_at(2, 0), std::out_of_range);
            }
        }
    }
}

SCENARIO("bqms can be resized", "[bqm]") {
    GIVEN("A bqm with some entries") {
        float Q[16] = {1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1};
        auto bqm = AdjVectorBQM<float>(Q, 4);

        REQUIRE(bqm.num_variables() == 4);
        REQUIRE(bqm.num_interactions() == 4);

        WHEN("the size is increased") {
            bqm.resize(10);

            THEN("more disconnected variables are added") {
                REQUIRE(bqm.num_variables() == 10);
                REQUIRE(bqm.num_interactions() == 4);

                REQUIRE(bqm.linear(5) == 0);
                REQUIRE(bqm.linear(6) == 0);
                REQUIRE(bqm.linear(7) == 0);
            }
        }

        WHEN("the size is decreased by 1") {
            bqm.resize(3);

            THEN("a variable and its interactions are removed") {
                REQUIRE(bqm.num_variables() == 3);
                REQUIRE(bqm.num_interactions() == 2);
            }
        }

        WHEN("the size is descreased by 2") {
            bqm.resize(2);

            THEN("two variables and their interactions are removed") {
                REQUIRE(bqm.num_variables() == 2);
                REQUIRE(bqm.num_interactions() == 0);
            }
        }
    }
}

SCENARIO("quadratic biases can be manipulated", "[bqm]") {
    GIVEN("A bqm with some variables but no interactions") {
        auto bqm = AdjVectorBQM<float>();
        bqm.resize(10);

        REQUIRE(bqm.num_variables() == 10);
        REQUIRE(bqm.num_interactions() == 0);

        WHEN("set_quadratic() is called") {
            bqm.set_quadratic(2, 3, .5);

            THEN("the quadratic bias is set") {
                REQUIRE(bqm.quadratic(2, 3) == .5);
                REQUIRE(bqm.quadratic(3, 2) == .5);
                REQUIRE(bqm.num_interactions() == 1);
            }

            AND_WHEN("add_interaction() is applied to the same interaction") {
                bqm.add_interaction(2, 3, 1);

                THEN("the biases are summed") {
                    REQUIRE(bqm.quadratic(2, 3) == 1.5);
                    REQUIRE(bqm.quadratic(3, 2) == 1.5);
                    REQUIRE(bqm.num_interactions() == 1);
                }
            }

            AND_WHEN("remove_interaction() is applied") {
                REQUIRE(bqm.remove_interaction(2, 3));

                THEN("there are no interactions") {
                    REQUIRE(bqm.num_interactions() == 0);
                }
            }
        }

        WHEN("add_interaction() is called on a non-existant interaciton") {
            bqm.add_interaction(2, 3, .5);

            THEN("the quadratic bias is set") {
                REQUIRE(bqm.quadratic(2, 3) == .5);
                REQUIRE(bqm.quadratic(3, 2) == .5);
                REQUIRE(bqm.num_interactions() == 1);
            }

            AND_WHEN("set_quadratic() is applied to the same interaction") {
                bqm.set_quadratic(2, 3, 1);

                THEN("the biases are overwritten") {
                    REQUIRE(bqm.quadratic(2, 3) == 1);
                    REQUIRE(bqm.quadratic(3, 2) == 1);
                    REQUIRE(bqm.num_interactions() == 1);
                }
            }
        }
    }

    GIVEN("A dense bqm with some variables") {
        std::vector<int> Q(100, 1);
        auto bqm = AdjVectorBQM<float>(&Q[0], 10);

        REQUIRE(bqm.num_variables() == 10);
        REQUIRE(bqm.num_interactions() == 45);

        WHEN("set_quadratic() is called") {
            bqm.set_quadratic(2, 3, .5);

            THEN("the quadratic bias is set") {
                REQUIRE(bqm.quadratic(2, 3) == .5);
                REQUIRE(bqm.quadratic(3, 2) == .5);
                REQUIRE(bqm.num_interactions() == 45);
            }

            AND_WHEN("add_interaction() is applied to the same interaction") {
                bqm.add_interaction(2, 3, 1);

                THEN("the biases are summed") {
                    REQUIRE(bqm.quadratic(2, 3) == 1.5);
                    REQUIRE(bqm.quadratic(3, 2) == 1.5);
                    REQUIRE(bqm.num_interactions() == 45);
                }
            }

            AND_WHEN("remove_interaction() is applied") {
                REQUIRE(bqm.remove_interaction(2, 3));

                THEN("there is one less interactions") {
                    REQUIRE(bqm.num_interactions() == 44);
                }

                AND_WHEN(
                        "remove_interaction() is applied again the edge is "
                        "gone") {
                    REQUIRE(!bqm.remove_interaction(2, 3));
                }
            }
        }

        WHEN("add_interaction() is called") {
            bqm.add_interaction(2, 3, .5);

            THEN("the quadratic bias is added") {
                REQUIRE(bqm.quadratic(2, 3) == 2.5);
                REQUIRE(bqm.quadratic(3, 2) == 2.5);
                REQUIRE(bqm.num_interactions() == 45);
            }

            AND_WHEN("set_quadratic() is applied to the same interaction") {
                bqm.set_quadratic(2, 3, 1);

                THEN("the biases are overwritten") {
                    REQUIRE(bqm.quadratic(2, 3) == 1);
                    REQUIRE(bqm.quadratic(3, 2) == 1);
                    REQUIRE(bqm.num_interactions() == 45);
                }
            }
        }

        WHEN("remove_interaction is called for each edge") {
            for (auto u = 0; u < bqm.num_variables() - 1; ++u) {
                for (auto v = u + 1; v < bqm.num_variables(); ++v) {
                    REQUIRE(bqm.remove_interaction(u, v));
                }
            }

            THEN("the bqm is now empty") {
                REQUIRE(bqm.num_interactions() == 0);
            }
        }
    }
}

TEMPLATE_TEST_CASE("Scenario: bqms can be constructed from dense arrays",
                   "[bqm]", std::int8_t, std::int16_t, std::int32_t,
                   std::int64_t, float, double) {
    GIVEN("A matrix as an array") {
        TestType Q[9] = {1, 0, 3, 2, 1, 0, 1, 0, 0};
        int num_variables = 3;

        WHEN("a bqm with float biases is constructed") {
            auto bqm = AdjVectorBQM<float>(Q, num_variables);

            THEN("it gets its linear from the diagonal") {
                REQUIRE(bqm.num_variables() == 3);
                REQUIRE(bqm.linear(0) == 1);
                REQUIRE(bqm.linear(1) == 1);
                REQUIRE(bqm.linear(2) == 0);
            }

            THEN("it gets its quadratic from the off-diagonal") {
                REQUIRE(bqm.num_interactions() == 2);

                // test both forward and backward
                REQUIRE(bqm.quadratic(0, 1) == 2);
                REQUIRE(bqm.quadratic(1, 0) == 2);
                REQUIRE(bqm.quadratic(0, 2) == 4);
                REQUIRE(bqm.quadratic(2, 0) == 4);
                REQUIRE(bqm.quadratic(1, 2) == 0);
                REQUIRE(bqm.quadratic(2, 1) == 0);

                // ignores 0s
                REQUIRE_THROWS_AS(bqm.quadratic_at(1, 2), std::out_of_range);
                REQUIRE_THROWS_AS(bqm.quadratic_at(2, 1), std::out_of_range);
            }
        }

        WHEN("a bqm with float biases is constructed ignoring the diagonal") {
            auto bqm = AdjVectorBQM<float>(Q, num_variables, true);

            THEN("it has 0 bias on the diagonal") {
                REQUIRE(bqm.num_variables() == 3);
                REQUIRE(bqm.linear(0) == 0);
                REQUIRE(bqm.linear(1) == 0);
                REQUIRE(bqm.linear(2) == 0);
            }

            AND_THEN("it gets its quadratic from the off-diagonal") {
                REQUIRE(bqm.num_interactions() == 2);

                // test both forward and backward
                REQUIRE(bqm.quadratic(0, 1) == 2);
                REQUIRE(bqm.quadratic(1, 0) == 2);
                REQUIRE(bqm.quadratic(0, 2) == 4);
                REQUIRE(bqm.quadratic(2, 0) == 4);
                REQUIRE(bqm.quadratic(1, 2) == 0);
                REQUIRE(bqm.quadratic(2, 1) == 0);

                // ignores 0s
                REQUIRE_THROWS_AS(bqm.quadratic_at(1, 2), std::out_of_range);
                REQUIRE_THROWS_AS(bqm.quadratic_at(2, 1), std::out_of_range);
            }
        }

        WHEN("a bqm with int biases is constructed") {
            auto bqm = AdjVectorBQM<int>(Q, num_variables);

            THEN("it gets its linear from the diagonal") {
                REQUIRE(bqm.num_variables() == 3);
                REQUIRE(bqm.linear(0) == 1);
                REQUIRE(bqm.linear(1) == 1);
                REQUIRE(bqm.linear(2) == 0);
            }

            THEN("it gets its quadratic from the off-diagonal") {
                REQUIRE(bqm.num_interactions() == 2);

                // test both forward and backward
                REQUIRE(bqm.quadratic(0, 1) == 2);
                REQUIRE(bqm.quadratic(1, 0) == 2);
                REQUIRE(bqm.quadratic(0, 2) == 4);
                REQUIRE(bqm.quadratic(2, 0) == 4);
                REQUIRE(bqm.quadratic(1, 2) == 0);
                REQUIRE(bqm.quadratic(2, 1) == 0);

                // ignores 0s
                REQUIRE_THROWS_AS(bqm.quadratic_at(1, 2), std::out_of_range);
                REQUIRE_THROWS_AS(bqm.quadratic_at(2, 1), std::out_of_range);
            }
        }
    }
}

TEMPLATE_TEST_CASE("Scenario: energies can be calculated", "[bqm]", std::int8_t,
                   float) {
    GIVEN("a binary quadratic model") {
        auto bqm = AdjVectorBQM<float>();
        bqm.resize(5);
        bqm.linear(0) = 1.5;
        bqm.linear(1) = -1;
        bqm.linear(2) = 2;
        bqm.linear(3) = 8;
        bqm.set_quadratic(0, 3, -1);
        bqm.set_quadratic(3, 2, 5);
        bqm.offset = -5;

        AND_GIVEN("a sample as an array") {
            TestType sample[5] = {-1, 1, -1, 1, 1};

            THEN("the energy can be calculated") {
                REQUIRE(bqm.energy(&sample[0]) == -5.5);
            }
        }

        AND_GIVEN("a sample as a vector") {
            std::vector<TestType> sample{-1, 1, -1, 1, 1};

            THEN("the energy can be calculated") {
                REQUIRE(bqm.energy(&sample[0]) == -5.5);
            }
        }
    }
}

}  // namespace dimod
