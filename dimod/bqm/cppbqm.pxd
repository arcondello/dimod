# distutils: language = c++
# cython: language_level=3
#
# Copyright 2019 D-Wave Systems Inc.
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.
#
# =============================================================================

from libcpp cimport bool
from libcpp.map cimport map as cppmap
from libcpp.pair cimport pair
from libcpp.vector cimport vector

from dimod.bqm.common cimport VarIndex, Bias


# cdef extern from "src/adjarray.cc":
#     pass


#     size_t num_variables[V, B](const pair[vector[pair[size_t, B]],
#                                           vector[pair[V, B]]]&)
#     size_t num_interactions[V, B](const pair[vector[pair[size_t, B]],
#                                              vector[pair[V, B]]]&)

#     B get_linear[V, B](const pair[vector[pair[size_t, B]], vector[pair[V, B]]]&,
#                        V)
#     pair[B, bool] get_quadratic[V, B](const pair[vector[pair[size_t, B]], vector[pair[V, B]]]&,
#                                       V, V)

#     size_t degree[V, B](const pair[vector[pair[size_t, B]], vector[pair[V, B]]]&, V)

#     pair[vector[pair[V, B]].iterator, vector[pair[V, B]].iterator] neighborhood[V, B](pair[vector[pair[size_t, B]], vector[pair[V, B]]]&, V)
#     pair[vector[pair[V, B]].iterator, vector[pair[V, B]].iterator] neighborhood[V, B](pair[vector[pair[size_t, B]], vector[pair[V, B]]]&, V, bool)

#     void set_linear[V, B](pair[vector[pair[size_t, B]], vector[pair[V, B]]]&,
#                           V, B)
#     bool set_quadratic[V, B](pair[vector[pair[size_t, B]], vector[pair[V, B]]]&,
#                              V, V, B)

#     void copy_bqm[V, B, BQM](BQM&, pair[vector[pair[size_t, B]], vector[pair[V, B]]]&)


cdef extern from "src/shapeable.cc":
    pass

cdef extern from "src/shapeable.h" namespace "dimod" nogil:

    # some of these should have const bqm inputs but cython seems to have
    # trouble with that

    size_t num_variables[N, B](vector[pair[N, B]]&)
    size_t num_interactions[N, B](vector[pair[N, B]]&)

    B get_linear[N, V, B](vector[pair[N, B]]&, V)
    pair[B, bool] get_quadratic[N, V, B](vector[pair[N, B]]&, V, V)

    size_t degree[N, V, B](vector[pair[N, B]]&, V)

    pair[vector[pair[V, B]].iterator,
         vector[pair[V, B]].iterator] neighborhood[V, B](
        vector[pair[vector[pair[V, B]], B]]&, V)
    pair[cppmap[V, B].iterator,
         cppmap[V, B].iterator] neighborhood[V, B](
        vector[pair[cppmap[V, B], B]]&, V)

    void set_linear[N, V, B](vector[pair[N, B]]&, V, B)

    void set_quadratic[V, B](vector[pair[vector[pair[V, B]], B]]&, V, V, B)
    void set_quadratic[V, B](vector[pair[cppmap[V, B], B]]&, V, V, B)

    size_t add_variable[N, B](vector[pair[N, B]]&)


    V add_variable[V, B](vector[pair[vector[pair[V, B]], B]]&)

    void copy_bqm[V, B, BQM](BQM&, vector[pair[vector[pair[V, B]], B]]&)
    void copy_bqm[V, B, BQM](BQM&, vector[pair[cppmap[V, B], B]]&)

    size_t pop_variable[N, B](vector[pair[N, B]]&)

    bool remove_interaction[V, B](vector[pair[vector[pair[V, B]], B]]&, V, V)
    bool remove_interaction[V, B](vector[pair[cppmap[V, B], B]]&, V, V)


cdef extern from "dimod/adjarraybqm.h" namespace "dimod" nogil:

    cdef cppclass AdjArrayBQM[V, B]:
        ctypedef V variable_type
        ctypedef size_t neighborhood_type
        ctypedef B bias_type
        ctypedef size_t size_type

        vector[pair[neighborhood_type, bias_type]] invars
        vector[pair[variable_type, bias_type]] outvars

        cppclass outvars_iterator:
            pair[variable_type, bias_type]& operator*()
            outvars_iterator operator++()
            outvars_iterator operator--()
            outvars_iterator operator+(size_type)
            outvars_iterator operator-(size_type)
            size_t operator-(outvars_iterator)
            bint operator==(outvars_iterator)
            bint operator!=(outvars_iterator)
            bint operator<(outvars_iterator)
            bint operator>(outvars_iterator)
            bint operator<=(outvars_iterator)
            bint operator>=(outvars_iterator)


        AdjArrayBQM() except +
        AdjArrayBQM(AdjArrayBQM&) except +  # cython cannot handle template here

        size_type num_interactions() except +
        size_type num_variables() except +

        bias_type get_linear(variable_type) except +
        pair[bias_type, bool] get_quadratic(variable_type, variable_type) except +

        size_type degree(variable_type) except +

        pair[outvars_iterator, outvars_iterator] neighborhood(variable_type) except +

        void set_linear(variable_type, bias_type) except +
        bool set_quadratic(variable_type, variable_type, bias_type) except +

cdef extern from "dimod/shapeablebqm.h" namespace "dimod" nogil:

    cdef cppclass cppAdjMapBQM[V, B]:
        ctypedef V variable_type
        ctypedef B bias_type
        ctypedef size_t size_type

        vector[pair[cppmap[variable_type, bias_type], bias_type]] adj

        cppclass outvars_iterator:
            pair[variable_type, bias_type]& operator*()
            outvars_iterator operator++()
            outvars_iterator operator--()
            outvars_iterator operator+(size_type)
            outvars_iterator operator-(size_type)
            size_t operator-(outvars_iterator)
            bint operator==(outvars_iterator)
            bint operator!=(outvars_iterator)
            bint operator<(outvars_iterator)
            bint operator>(outvars_iterator)
            bint operator<=(outvars_iterator)
            bint operator>=(outvars_iterator)

        cppAdjMapBQM() except +

        # cython cannot handle template here, so we fix them. It's really
        # annoying
        cppAdjMapBQM(AdjArrayBQM[VarIndex, Bias]&) except +

        size_type num_interactions() except +
        size_type num_variables() except +

    cdef cppclass cppAdjVectorBQM[V, B]:
        ctypedef V variable_type
        ctypedef B bias_type
        ctypedef size_t size_type

        vector[pair[vector[pair[variable_type, bias_type]], bias_type]] adj

        cppclass outvars_iterator:
            pair[variable_type, bias_type]& operator*()
            outvars_iterator operator++()
            outvars_iterator operator--()
            outvars_iterator operator+(size_type)
            outvars_iterator operator-(size_type)
            size_t operator-(outvars_iterator)
            bint operator==(outvars_iterator)
            bint operator!=(outvars_iterator)
            bint operator<(outvars_iterator)
            bint operator>(outvars_iterator)
            bint operator<=(outvars_iterator)
            bint operator>=(outvars_iterator)

        cppAdjVectorBQM() except +
        # cppAdjVectorBQM(AdjArrayBQM&) except +  # cython cannot handle template here

        size_type num_interactions() except +
        size_type num_variables() except +
