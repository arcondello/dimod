# distutils: language = c++
# cython: language_level=3

# Copyright 2022 D-Wave Systems Inc.
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

cimport numpy as np

from dimod.binary cimport cyBQM_float32, cyBQM_float64
from dimod.cyvariables cimport cyVariables
from dimod.libcpp cimport *
from dimod.quadratic cimport cyQM_float32, cyQM_float64


cdef extern from "dimod/constrained.h" namespace "dimod" nogil:
    cdef cppclass cppConstrainedQuadraticModel "dimod::ConstrainedQuadraticModel" [Bias, Index]:
        ctypedef Bias bias_type
        ctypedef size_t size_type
        ctypedef Index index_type

        index_type add_variable(cppVartype, bias_type, bias_type)
        bias_type& lower_bound(index_type)
        # const bias_type& lower_bound(index_type) const
        const cppQuadraticModel[Bias, Index]& objective() const
        void set_objective[B, I](const cppQuadraticModel[B, I]&)
        void set_objective[B, I](const cppBinaryQuadraticModel[B, I]&)
        const cppVartype& vartype(index_type) const
        bias_type& upper_bound(index_type)
        # const bias_type& upper_bound(index_type) const


ctypedef np.float64_t bias_type
ctypedef np.int32_t index_type


ctypedef fused cyQM:
    cyBQM_float32
    cyBQM_float64
    cyQM_float32
    cyQM_float64


cdef class cyConstrainedQuadraticModel:
    cdef cppConstrainedQuadraticModel[bias_type, index_type] data

    cdef readonly constraint_labels
    cdef readonly object dtype
    cdef readonly object index_dtype
    cdef readonly cyVariables variables

    cdef public int REAL_INTERACTIONS

    # def _set_objective(self, cyQM objective)
