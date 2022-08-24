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

        THEN("neighborhood iteration should return nothing") {
            for (std::size_t v = 0; v < bqm.num_variables(); ++v) {
                CHECK(bqm.cbegin_neighborhood(v) == bqm.cend_neighborhood(v));
            }
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

TEST_CASE("BinaryQuadraticModel swap") {
    GIVEN("two BQMs of different vartypes") {
        auto bqm1 = BinaryQuadraticModel<double>(3, Vartype::BINARY);
        bqm1.set_offset(2);
        bqm1.set_linear(0, {1, 2, 3});
        bqm1.set_quadratic(0, 1, 1);
        bqm1.set_quadratic(1, 2, 2);

        auto bqm2 = BinaryQuadraticModel<double>(3, Vartype::SPIN);
        bqm2.set_offset(-2);
        bqm2.set_linear(0, {-1, -2, -3});
        bqm2.set_quadratic(0, 1, -1);
        bqm2.set_quadratic(1, 2, -2);

        WHEN("we swap them") {
            std::swap(bqm1, bqm2);

            THEN("everything is moved appropriately") {
                CHECK(bqm1.linear(0) == -1);
                CHECK(bqm1.linear(1) == -2);
                CHECK(bqm1.linear(2) == -3);
                CHECK(bqm2.linear(0) == 1);
                CHECK(bqm2.linear(1) == 2);
                CHECK(bqm2.linear(2) == 3);

                CHECK(bqm1.quadratic(0, 1) == -1);
                CHECK(bqm1.quadratic(1, 2) == -2);
                CHECK(bqm2.quadratic(0, 1) == 1);
                CHECK(bqm2.quadratic(1, 2) == 2);

                CHECK(bqm1.offset() == -2);
                CHECK(bqm2.offset() == 2);

                CHECK(bqm1.vartype() == Vartype::SPIN);
                CHECK(bqm2.vartype() == Vartype::BINARY);
            }
        }
    }
}

TEMPLATE_TEST_CASE_SIG("Scenario: BinaryQuadraticModel tests", "[qmbase][bqm]",
                       ((typename Bias, Vartype vartype), Bias, vartype), (double, Vartype::BINARY),
                       (double, Vartype::SPIN), (float, Vartype::BINARY), (float, Vartype::SPIN)) {
    GIVEN("an empty BQM") {
        auto bqm = BinaryQuadraticModel<Bias>(vartype);

        WHEN("the bqm is resized") {
            bqm.resize(10);

            THEN("it will have the correct number of variables with 0 bias") {
                REQUIRE(bqm.num_variables() == 10);
                REQUIRE(bqm.num_interactions() == 0);
                for (auto v = 0u; v < bqm.num_variables(); ++v) {
                    REQUIRE(bqm.linear(v) == 0);
                }
            }
        }

        // AND_GIVEN("some COO-formatted arrays") {
        //     int irow[4] = {0, 2, 0, 1};
        //     int icol[4] = {0, 2, 1, 2};
        //     float bias[4] = {.5, -2, 2, -3};
        //     std::size_t length = 4;

        //     WHEN("we add the biases with add_quadratic") {
        //         bqm.add_quadratic(&irow[0], &icol[0], &bias[0], length);

        //         THEN("it takes its values from the arrays") {
        //             REQUIRE(bqm.num_variables() == 3);

        //             if (bqm.vartype() == Vartype::SPIN) {
        //                 REQUIRE(bqm.linear(0) == 0);
        //                 REQUIRE(bqm.linear(1) == 0);
        //                 REQUIRE(bqm.linear(2) == 0);
        //                 REQUIRE(bqm.offset() == -1.5);
        //             } else {
        //                 REQUIRE(bqm.vartype() == Vartype::BINARY);
        //                 REQUIRE(bqm.linear(0) == .5);
        //                 REQUIRE(bqm.linear(1) == 0);
        //                 REQUIRE(bqm.linear(2) == -2);
        //                 REQUIRE(bqm.offset() == 0);
        //             }

        //             REQUIRE(bqm.num_interactions() == 2);
        //             REQUIRE(bqm.quadratic(0, 1) == 2);
        //             REQUIRE(bqm.quadratic(2, 1) == -3);
        //             REQUIRE_THROWS_AS(bqm.quadratic_at(0, 2), std::out_of_range);
        //         }
        //     }
        // }

        // AND_GIVEN("some COO-formatted arrays with duplicates") {
        //     int irow[6] = {0, 2, 0, 1, 0, 0};
        //     int icol[6] = {0, 2, 1, 2, 1, 0};
        //     float bias[6] = {.5, -2, 2, -3, 4, 1};
        //     std::size_t length = 6;

        //     WHEN("we add the biases with add_quadratic") {
        //         bqm.add_quadratic(&irow[0], &icol[0], &bias[0], length);

        //         THEN("it combines duplicate values") {
        //             REQUIRE(bqm.num_variables() == 3);

        //             if (bqm.vartype() == Vartype::SPIN) {
        //                 REQUIRE(bqm.linear(0) == 0);
        //                 REQUIRE(bqm.linear(1) == 0);
        //                 REQUIRE(bqm.linear(2) == 0);
        //                 REQUIRE(bqm.offset() == -.5);
        //             } else {
        //                 REQUIRE(bqm.vartype() == Vartype::BINARY);
        //                 REQUIRE(bqm.linear(0) == 1.5);
        //                 REQUIRE(bqm.linear(1) == 0.);
        //                 REQUIRE(bqm.linear(2) == -2);
        //                 REQUIRE(bqm.offset() == 0);
        //             }

        //             REQUIRE(bqm.num_interactions() == 2);
        //             REQUIRE(bqm.quadratic(0, 1) == 6);
        //             REQUIRE(bqm.quadratic(2, 1) == -3);
        //             REQUIRE_THROWS_AS(bqm.quadratic_at(0, 2), std::out_of_range);
        //         }
        //     }
        // }

        // AND_GIVEN("some COO-formatted arrays with multiple duplicates") {
        //     int irow[4] = {0, 1, 0, 1};
        //     int icol[4] = {1, 2, 1, 0};
        //     float bias[4] = {-1, 1, -2, -3};
        //     std::size_t length = 4;

        //     WHEN("we add the biases with add_quadratic") {
        //         bqm.add_quadratic(&irow[0], &icol[0], &bias[0], length);

        //         THEN("it combines duplicate values") {
        //             REQUIRE(bqm.num_variables() == 3);
        //             REQUIRE(bqm.linear(0) == 0);
        //             REQUIRE(bqm.linear(1) == 0);
        //             REQUIRE(bqm.linear(2) == 0);

        //             REQUIRE(bqm.num_interactions() == 2);
        //             REQUIRE(bqm.quadratic(0, 1) == -6);
        //             REQUIRE(bqm.quadratic(1, 0) == -6);
        //             REQUIRE(bqm.quadratic(2, 1) == 1);
        //             REQUIRE(bqm.quadratic(1, 2) == 1);
        //             REQUIRE_THROWS_AS(bqm.quadratic_at(0, 2), std::out_of_range);
        //             REQUIRE_THROWS_AS(bqm.quadratic_at(2, 0), std::out_of_range);
        //         }
        //     }
        // }
    }

    GIVEN("a BQM constructed from a dense array") {
        float Q[9] = {1, 0, 3, 2, 1, 0, 1, 0, 0};
        int num_variables = 3;

        auto bqm = BinaryQuadraticModel<Bias>(Q, num_variables, vartype);

        THEN("it handles the diagonal according to its vartype") {
            REQUIRE(bqm.num_variables() == 3);

            if (bqm.vartype() == Vartype::SPIN) {
                REQUIRE(bqm.linear(0) == 0);
                REQUIRE(bqm.linear(1) == 0);
                REQUIRE(bqm.linear(2) == 0);
                REQUIRE(bqm.offset() == 2);
            } else {
                REQUIRE(bqm.vartype() == Vartype::BINARY);
                REQUIRE(bqm.linear(0) == 1);
                REQUIRE(bqm.linear(1) == 1);
                REQUIRE(bqm.linear(2) == 0);
                REQUIRE(bqm.offset() == 0);
            }
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

        THEN("we can iterate over the neighborhood of 0") {
            std::size_t count = 0;
            for (auto it = bqm.cbegin_neighborhood(0); it != bqm.cend_neighborhood(0);
                 ++it, ++count) {
                CHECK(bqm.quadratic_at(0, it->v) == it->bias);
            }
            CHECK(count == 2);
        }

        WHEN("we iterate over the quadratic biases") {
            auto first = bqm.cbegin_quadratic();
            auto last = bqm.cend_quadratic();
            THEN("we read out the lower triangle") {
                CHECK(first->u == 1);
                CHECK(first->v == 0);
                CHECK(first->bias == 2);
                CHECK((*first).u == 1);
                CHECK((*first).v == 0);
                CHECK((*first).bias == 2);

                ++first;

                CHECK(first->u == 2);
                CHECK(first->v == 0);
                CHECK(first->bias == 4);

                first++;

                CHECK(first == last);
            }
        }
    }

    // GIVEN("a BQM with five variables, two interactions and an offset") {
    //     auto bqm = BinaryQuadraticModel<Bias>(5, vartype);
    //     bqm.set_linear({1, -3.25, 0, 3, -4.5});
    //     bqm.set_quadratic(0, 3, -1);
    //     bqm.set_quadratic(3, 1, 5.6);
    //     bqm.set_quadratic(0, 1, 1.6);
    //     bqm.set_offset(-3.8);

    //     AND_GIVEN("the set of all possible five variable samples") {
    //         // there are smarter ways to do this but it's simple
    //         std::vector<std::vector<int>> spn_samples;
    //         std::vector<std::vector<int>> bin_samples;
    //         for (auto i = 0; i < 1 << bqm.num_variables(); ++i) {
    //             std::vector<int> bin_sample;
    //             std::vector<int> spn_sample;
    //             for (size_t v = 0; v < bqm.num_variables(); ++v) {
    //                 bin_sample.push_back((i >> v) & 1);
    //                 spn_sample.push_back(2 * ((i >> v) & 1) - 1);
    //             }

    //             bin_samples.push_back(bin_sample);
    //             spn_samples.push_back(spn_sample);
    //         }

    //         std::vector<double> energies;
    //         if (vartype == Vartype::SPIN) {
    //             for (auto& sample : spn_samples) {
    //                 energies.push_back(bqm.energy(sample));
    //             }
    //         } else {
    //             for (auto& sample : bin_samples) {
    //                 energies.push_back(bqm.energy(sample));
    //             }
    //         }

    //         WHEN("we change the vartype to spin") {
    //             bqm.change_vartype(Vartype::SPIN);

    //             THEN("the energies will match") {
    //                 for (size_t si = 0; si < energies.size(); ++si) {
    //                     REQUIRE(energies[si] == Approx(bqm.energy(spn_samples[si])));
    //                 }
    //             }
    //         }

    //         WHEN("we change the vartype to binary") {
    //             bqm.change_vartype(Vartype::BINARY);
    //             THEN("the energies will match") {
    //                 for (size_t si = 0; si < energies.size(); ++si) {
    //                     REQUIRE(energies[si] == Approx(bqm.energy(bin_samples[si])));
    //                 }
    //             }
    //         }
    //     }
    // }
}

// TEMPLATE_TEST_CASE_SIG("Scenario: BQMs can be combined", "[bqm]",
//                        ((typename B0, typename B1, Vartype vartype), B0, B1, vartype),
//                        (float, float, Vartype::BINARY), (double, float, Vartype::BINARY),
//                        (float, double, Vartype::BINARY), (double, double, Vartype::BINARY),
//                        (float, float, Vartype::SPIN), (double, float, Vartype::SPIN),
//                        (float, double, Vartype::SPIN), (double, double, Vartype::SPIN)) {
//     GIVEN("a BQM with 3 variables") {
//         auto bqm0 = BinaryQuadraticModel<B0>(3, vartype);
//         bqm0.linear(2) = -1;
//         bqm0.set_quadratic(0, 1, 1.5);
//         bqm0.set_quadratic(0, 2, -2);
//         bqm0.set_quadratic(1, 2, 7);
//         bqm0.offset() = -4;

//         AND_GIVEN("a BQM with 5 variables and the same vartype") {
//             auto bqm1 = BinaryQuadraticModel<B1>(5, vartype);
//             bqm1.linear(0) = 1;
//             bqm1.linear(1) = -3.25;
//             bqm1.linear(2) = 2;
//             bqm1.linear(3) = 3;
//             bqm1.linear(4) = -4.5;
//             bqm1.set_quadratic(0, 1, 5.6);
//             bqm1.set_quadratic(0, 3, -1);
//             bqm1.set_quadratic(1, 2, 1.6);
//             bqm1.set_quadratic(3, 4, -25);
//             bqm1.offset() = -3.8;

//             WHEN("the first is updated with the second") {
//                 bqm0.add_bqm(bqm1);

//                 THEN("the biases are added") {
//                     REQUIRE(bqm0.num_variables() == 5);
//                     REQUIRE(bqm0.num_interactions() == 5);

//                     CHECK(bqm0.offset() == Approx(-7.8));

//                     CHECK(bqm0.linear(0) == Approx(1));
//                     CHECK(bqm0.linear(1) == Approx(-3.25));
//                     CHECK(bqm0.linear(2) == Approx(1));
//                     CHECK(bqm0.linear(3) == Approx(3));
//                     CHECK(bqm0.linear(4) == Approx(-4.5));

//                     CHECK(bqm0.quadratic(0, 1) == Approx(7.1));
//                     CHECK(bqm0.quadratic(0, 2) == Approx(-2));
//                     CHECK(bqm0.quadratic(0, 3) == Approx(-1));
//                     CHECK(bqm0.quadratic(1, 2) == Approx(8.6));
//                     CHECK(bqm0.quadratic(3, 4) == Approx(-25));
//                 }
//             }

//             WHEN("the second is updated with the first") {
//                 bqm1.add_bqm(bqm0);

//                 THEN("the biases are added") {
//                     REQUIRE(bqm1.num_variables() == 5);
//                     REQUIRE(bqm1.num_interactions() == 5);

//                     CHECK(bqm1.offset() == Approx(-7.8));

//                     CHECK(bqm1.linear(0) == Approx(1));
//                     CHECK(bqm1.linear(1) == Approx(-3.25));
//                     CHECK(bqm1.linear(2) == Approx(1));
//                     CHECK(bqm1.linear(3) == Approx(3));
//                     CHECK(bqm1.linear(4) == Approx(-4.5));

//                     CHECK(bqm1.quadratic(0, 1) == Approx(7.1));
//                     CHECK(bqm1.quadratic(0, 2) == Approx(-2));
//                     CHECK(bqm1.quadratic(0, 3) == Approx(-1));
//                     CHECK(bqm1.quadratic(1, 2) == Approx(8.6));
//                     CHECK(bqm1.quadratic(3, 4) == Approx(-25));
//                 }
//             }
//         }

//         AND_GIVEN("a BQM with 5 variables and a different vartype") {
//             Vartype vt;
//             if (vartype == Vartype::SPIN) {
//                 vt = Vartype::BINARY;
//             } else {
//                 vt = Vartype::SPIN;
//             }

//             auto bqm1 = BinaryQuadraticModel<B1>(5, vt);
//             bqm1.linear(0) = 1;
//             bqm1.linear(1) = -3.25;
//             bqm1.linear(2) = 2;
//             bqm1.linear(3) = 3;
//             bqm1.linear(4) = -4.5;
//             bqm1.set_quadratic(0, 1, 5.6);
//             bqm1.set_quadratic(0, 3, -1);
//             bqm1.set_quadratic(1, 2, 1.6);
//             bqm1.set_quadratic(3, 4, -25);
//             bqm1.offset() = -3.8;

//             WHEN("the first is updated with the second") {
//                 auto bqm0_cp = BinaryQuadraticModel<B0>(bqm0);
//                 auto bqm1_cp = BinaryQuadraticModel<B1>(bqm1);

//                 bqm0.add_bqm(bqm1);

//                 THEN("it was as if the vartype was changed first") {
//                     bqm1_cp.change_vartype(vartype);
//                     bqm0_cp.add_bqm(bqm1_cp);

//                     REQUIRE(bqm0.num_variables() == bqm0_cp.num_variables());
//                     REQUIRE(bqm0.num_interactions() == bqm0_cp.num_interactions());
//                     REQUIRE(bqm0.offset() == Approx(bqm0_cp.offset()));
//                     for (auto u = 0u; u < bqm0.num_variables(); ++u) {
//                         REQUIRE(bqm0.linear(u) == Approx(bqm0_cp.linear(u)));

//                         auto span = bqm0.neighborhood(u);
//                         for (auto it = span.first; it != span.second; ++it) {
//                             REQUIRE((*it).second == Approx(bqm0_cp.quadratic_at(u, (*it).first)));
//                         }
//                     }
//                 }
//             }
//         }
//     }
// }

// SCENARIO("One bqm can be added to another") {
//     GIVEN("Two BQMs of different vartypes") {
//         auto bin = BinaryQuadraticModel<double>(2, Vartype::BINARY);
//         bin.linear(0) = .3;
//         bin.set_quadratic(0, 1, -1);

//         auto spn = BinaryQuadraticModel<double>(2, Vartype::SPIN);
//         spn.linear(1) = -1;
//         spn.set_quadratic(0, 1, 1);
//         spn.offset() = 1.2;

//         WHEN("the spin one is added to the binary") {
//             bin.add_bqm(spn);

//             THEN("the combined model is correct") {
//                 REQUIRE(bin.num_variables() == 2);
//                 CHECK(bin.num_interactions() == 1);

//                 CHECK(bin.linear(0) == -1.7);
//                 CHECK(bin.linear(1) == -4);

//                 CHECK(bin.quadratic(0, 1) == 3);
//             }
//         }

//         WHEN("the spin one is added to the binary one with an offset") {
//             std::vector<int> mapping = {1, 2};
//             bin.add_bqm(spn, mapping);

//             THEN("the combined model is correct") {
//                 REQUIRE(bin.num_variables() == 3);
//                 CHECK(bin.num_interactions() == 2);

//                 CHECK(bin.linear(0) == .3);
//                 CHECK(bin.linear(1) == -2);
//                 CHECK(bin.linear(2) == -4);

//                 CHECK(bin.quadratic(0, 1) == -1);
//                 CHECK(bin.quadratic(1, 2) == 4);
//                 CHECK_THROWS_AS(bin.quadratic_at(0, 2), std::out_of_range);
//             }
//         }
//     }

//     GIVEN("Two BQMs of the same vartype") {
//         auto bqm0 = BinaryQuadraticModel<double>(2, Vartype::SPIN);
//         bqm0.linear(0) = -1;
//         bqm0.set_quadratic(0, 1, 1.5);
//         bqm0.offset() = 1.5;

//         auto bqm1 = BinaryQuadraticModel<double>(3, Vartype::SPIN);
//         bqm1.linear(0) = -2;
//         bqm1.linear(2) = 3;
//         bqm1.set_quadratic(0, 1, 5);
//         bqm1.set_quadratic(2, 1, 1);
//         bqm1.offset() = 1.5;

//         WHEN("the spin one is added to the binary one with a mapping") {
//             std::vector<int> mapping = {0, 1, 2};
//             bqm0.add_bqm(bqm1, mapping);

//             THEN("the biases are summed") {
//                 REQUIRE(bqm0.num_variables() == 3);
//                 CHECK(bqm0.num_interactions() == 2);

//                 CHECK(bqm0.offset() == 3);

//                 CHECK(bqm0.linear(0) == -3);
//                 CHECK(bqm0.linear(1) == 0);
//                 CHECK(bqm0.linear(2) == 3);

//                 CHECK(bqm0.quadratic(0, 1) == 6.5);
//                 CHECK(bqm0.quadratic(1, 2) == 1);
//                 CHECK_THROWS_AS(bqm0.quadratic_at(0, 2), std::out_of_range);
//             }
//         }
//     }
// }


}  // namespace dimod
