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
#include "dimod/quadratic_model.h"

namespace dimod {

TEST_CASE("BinaryQuadraticModel constructors") {
    GIVEN("a linear BQM") {
        auto bqm = BinaryQuadraticModel<double>(5, Vartype::SPIN);
        bqm.set_linear(0, {0, -1, +2, -3, +4});
        bqm.set_offset(5);

        AND_GIVEN("another BQM with the same values") {
            auto expected = BinaryQuadraticModel<double>(5, Vartype::SPIN);
            expected.set_linear(0, {0, -1, +2, -3, +4});
            expected.set_offset(5);

            WHEN("another BQM is constructed from it using the assignment operator") {
                auto cpy = BinaryQuadraticModel<double>(bqm);

                THEN("all of the values are copied appropriately") {
                    CHECK(cpy.is_equal(expected));
                }

                AND_WHEN("the original BQM is modified") {
                    bqm.set_linear(1, 1);
                    bqm.add_quadratic(1, 2, 2);
                    bqm.set_offset(-5);

                    THEN("the copy is not also modified") {
                        CHECK(cpy.is_equal(expected));
                        CHECK(!bqm.is_equal(expected));  // sanity check
                    }
                }
            }

            WHEN("another BQM is constructed from it using the copy constructor") {
                BinaryQuadraticModel<double> cpy(bqm);

                THEN("all of the values are copied appropriately") {
                    CHECK(cpy.is_equal(expected));
                }

                AND_WHEN("the original BQM is modified") {
                    bqm.set_linear(1, 1);
                    bqm.add_quadratic(1, 2, 2);
                    bqm.set_offset(-5);

                    THEN("the copy is not also modified") {
                        CHECK(cpy.is_equal(expected));
                        CHECK(!bqm.is_equal(expected));  // sanity check
                    }
                }
            }
        }
    }
}

TEST_CASE("BinaryQuadraticModel quadratic iteration") {
    GIVEN("a linear BQM") {
        auto bqm = BinaryQuadraticModel<double>(5, Vartype::SPIN);
        bqm.set_linear(0, {0, -1, +2, -3, +4});
        bqm.set_offset(5);

        THEN("quadric iteration should return nothing") {
            CHECK(bqm.cbegin_quadratic() == bqm.cend_quadratic());
        }
    }

    GIVEN("a linear BQM that was once quadratic") {
        // as a technical detail, as of August 2022, this leaves the structure
        // for an adjacency in-place. If that changes in the future and this
        // test isn't removed then it does no harm
        auto bqm = BinaryQuadraticModel<double>(5, Vartype::SPIN);
        bqm.set_linear(0, {0, -1, +2, -3, +4});
        bqm.set_offset(5);
        bqm.add_quadratic(0, 1, 1);
        bqm.remove_interaction(0, 1);

        THEN("quadric iteration should return nothing") {
            CHECK(bqm.cbegin_quadratic() == bqm.cend_quadratic());
        }
    }

    GIVEN("a BQM with one quadratic interaction") {
        auto bqm = BinaryQuadraticModel<double>(5, Vartype::SPIN);
        bqm.set_linear(0, {0, -1, +2, -3, +4});
        bqm.set_offset(5);
        bqm.add_quadratic(1, 3, 5);

        THEN("our quadratic iterator has one value in it") {
            std::vector<int> row;
            std::vector<int> col;
            std::vector<double> biases;
            for (auto it = bqm.cbegin_quadratic(); it != bqm.cend_quadratic(); ++it) {
                row.push_back(it->u);
                col.push_back(it->v);
                biases.push_back(it->bias);
            }
            CHECK(row == std::vector<int>{3});
            CHECK(col == std::vector<int>{1});
            CHECK(biases == std::vector<double>{5});
        }
    }

    GIVEN("a BQM with several quadratic interactions") {
        auto bqm = BinaryQuadraticModel<double>(5, Vartype::SPIN);
        bqm.set_linear(0, {0, -1, +2, -3, +4});
        bqm.set_offset(5);
        bqm.add_quadratic({0, 0, 0, 3, 3}, {2, 3, 4, 4, 1}, {1, 2, 3, 4, 5});

        THEN("our quadratic iterator will correctly return the lower triangle") {
            std::vector<int> row;
            std::vector<int> col;
            std::vector<double> biases;
            // also sneak in a post-fix operator test
            for (auto it = bqm.cbegin_quadratic(); it != bqm.cend_quadratic(); it++) {
                row.push_back(it->u);
                col.push_back(it->v);
                biases.push_back(it->bias);
            }
            CHECK(row == std::vector<int>{2, 3, 3, 4, 4});
            CHECK(col == std::vector<int>{0, 0, 1, 0, 3});
            CHECK(biases == std::vector<double>{1, 2, 5, 3, 4});
        }
    }
}

TEST_CASE("BinaryQuadraticModel scale") {
    GIVEN("a bqm with linear, quadratic interactions and an offset") {
        auto bqm = BinaryQuadraticModel<double>(3, Vartype::BINARY);
        bqm.set_offset(2);
        bqm.set_linear(0, {1, 2, 3});
        bqm.set_quadratic(0, 1, 1);
        bqm.set_quadratic(1, 2, 2);

        WHEN("we scale it") {
            bqm.scale(5.5);

            THEN("all biases and the offset are scaled") {
                REQUIRE(bqm.offset() == 11.0);
                REQUIRE(bqm.linear(0) == 5.5);
                REQUIRE(bqm.linear(1) == 11.0);
                REQUIRE(bqm.linear(2) == 16.5);
                REQUIRE(bqm.quadratic(0, 1) == 5.5);
                REQUIRE(bqm.quadratic(1, 2) == 11.0);
            }
        }
    }
}

TEST_CASE("BinaryQuadraticModel vartype and bounds") {
    GIVEN("an empty BINARY BQM") {
        auto bqm = BinaryQuadraticModel<double>(Vartype::BINARY);

        THEN("the vartype and bounds are binary") {
            CHECK(bqm.vartype() == Vartype::BINARY);
            CHECK(bqm.lower_bound() == 0);
            CHECK(bqm.upper_bound() == 1);
        }
    }

    GIVEN("an empty SPIN BQM") {
        auto bqm = BinaryQuadraticModel<double>(Vartype::SPIN);

        THEN("the vartype and bounds are spin") {
            CHECK(bqm.vartype() == Vartype::SPIN);
            CHECK(bqm.lower_bound() == -1);
            CHECK(bqm.upper_bound() == +1);
        }
    }

    GIVEN("an empty default BQM") {
        auto bqm = BinaryQuadraticModel<double>();

        THEN("the vartype and bounds are binary") {
            CHECK(bqm.vartype() == Vartype::BINARY);
            CHECK(bqm.lower_bound() == 0);
            CHECK(bqm.upper_bound() == 1);
        }
    }
}

}  // namespace dimod
