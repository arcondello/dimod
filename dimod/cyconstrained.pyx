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

import warnings

from cython.operator cimport preincrement as inc, dereference as deref

from libc.math cimport ceil, floor

import numpy as np

import dimod

from dimod.cyutilities cimport cppvartype
from dimod.libcpp cimport cppvartype_info
from dimod.variables import Variables
from dimod.vartypes import Vartype, as_vartype, VartypeLike


BIAS_DTYPE = np.float64
INDEX_DTYPE = np.int32


cdef class cyConstrainedQuadraticModel:
    def __init__(self):
        self.dtype = np.dtype(BIAS_DTYPE)
        self.index_dtype = np.dtype(INDEX_DTYPE)
        self.variables = Variables()

        self.REAL_INTERACTIONS = dimod.REAL_INTERACTIONS


    def add_variable(self, vartype, v=None, *, lower_bound=None, upper_bound=None):
        cdef cppVartype vt = cppvartype(as_vartype(vartype, extended=True))

        cdef bias_type lb
        cdef bias_type ub

        cdef Py_ssize_t vi
        if v is not None and self.variables.count(v):
            # it already exists, so check that vartype matches
            vi = self.variables.index(v)
            if self.data.vartype(vi) != vt:
                raise TypeError(f"variable {v!r} already exists with a different vartype")
            if vt != cppVartype.BINARY and vt != cppVartype.SPIN:
                if lower_bound is not None:
                    lb = lower_bound
                    if lb != self.data.lower_bound(vi):
                        raise ValueError(
                            f"the specified lower bound, {lower_bound}, for "
                            f"variable {v!r} is different than the existing lower "
                            f"bound, {int(self.data.lower_bound(vi))}")
                if upper_bound is not None:
                    ub = upper_bound
                    if ub != self.data.upper_bound(vi):
                        raise ValueError(
                            f"the specified upper bound, {upper_bound}, for "
                            f"variable {v!r} is different than the existing upper "
                            f"bound, {int(self.data.upper_bound(vi))}")

            return v
        
        # the variable does not exist yet
        if vt == cppVartype.BINARY or vt == cppVartype.SPIN:
            # just ignore the provided values
            lb = cppvartype_info[bias_type].default_min(vt)
            ub = cppvartype_info[bias_type].default_max(vt)
        elif vt == cppVartype.INTEGER or vt == cppVartype.REAL:
            if lower_bound is None:
                lb = cppvartype_info[bias_type].default_min(vt)
            else:
                lb = lower_bound
                if lb < cppvartype_info[bias_type].min(vt):
                    raise ValueError(f"lower_bound cannot be less than {cppvartype_info[bias_type].min(vt)}")

            if upper_bound is None:
                ub = cppvartype_info[bias_type].default_max(vt)
            else:
                ub = upper_bound
                if ub > cppvartype_info[bias_type].max(vt):
                    raise ValueError(f"upper_bound cannot be greater than {cppvartype_info[bias_type].max(vt)}")
            
            if lb > ub:
                raise ValueError("lower_bound must be less than or equal to upper_bound")

            if vt == cppVartype.INTEGER and ceil(lb) > floor(ub):
                raise ValueError("there must be at least one valid integer between lower_bound and upper_bound")
        else:
            raise RuntimeError(f"unknown vartype: {vartype!s}")

        self.data.add_variable(vt, lb, ub)

        self.variables._append(v)

        return self.variables.at(-1)


    def set_lower_bound(self, v, bias_type lb):
        cdef Py_ssize_t vi = self.variables.index(v)
        cdef cppVartype cppvartype = self.data.vartype(vi)

        if cppvartype == cppVartype.BINARY or cppvartype == cppVartype.SPIN:
            raise ValueError(
                "cannot set the lower bound for BINARY or SPIN variables, "
                f"{v!r} is a {self.vartype(v).name} variable")

        if lb < cppvartype_info[bias_type].min(cppvartype):
            raise ValueError(f"lower_bound cannot be less than {cppvartype_info[bias_type].min(cppvartype)}")
            
        if lb > self.data.upper_bound(vi):
            raise ValueError(
                f"the specified lower bound, {lb}, cannot be set greater than the "
                f"current upper bound, {self.data.upper_bound(vi)}"
                )

        if cppvartype == cppVartype.INTEGER:
            if ceil(lb) > floor(self.data.upper_bound(vi)):
                raise ValueError(
                    "there must be at least one integer value between "
                    f"the specified lower bound, {lb} and the "
                    f"current upper bound, {self.data.upper_bound(vi)}"
                    )

        cdef bias_type *b = &(self.data.lower_bound(vi))
        b[0] = lb


    def set_objective(self, cyQM objective):

        if not self.data.objective().num_variables():
            self.data.set_objective(deref(objective.data()))
            self.variables._extend(objective.variables)
            return

        # there is something there already and we need to handle it

        raise NotImplementedError


    def set_upper_bound(self, v, bias_type ub):
        cdef Py_ssize_t vi = self.variables.index(v)
        cdef cppVartype cppvartype = self.data.vartype(vi)

        if cppvartype == cppVartype.BINARY or cppvartype == cppVartype.SPIN:
            raise ValueError(
                "cannot set the upper bound for BINARY or SPIN variables, "
                f"{v!r} is a {self.vartype(v).name} variable")

        if ub > cppvartype_info[bias_type].max(cppvartype):
            raise ValueError(f"upper_bound cannot be more than {cppvartype_info[bias_type].max(cppvartype)}")
            
        if ub < self.data.lower_bound(vi):
            raise ValueError(
                f"the specified upper bound, {ub}, cannot be set less than the "
                f"current lower bound, {self.data.lower_bound(vi)}"
                )

        if cppvartype == cppVartype.INTEGER:
            if ceil(self.data.lower_bound(vi)) > floor(ub):
                raise ValueError(
                    "there must be at least one integer value between "
                    f"the specified upper bound, {ub} and the "
                    f"current lower bound, {self.data.lower_bound(vi)}"
                    )

        cdef bias_type *b = &(self.data.upper_bound(vi))
        b[0] = ub


    def vartype(self, v):
        cdef Py_ssize_t vi = self.variables.index(v)
        cdef cppVartype cppvartype = self.data.vartype(vi)

        if cppvartype == cppVartype.BINARY:
            return Vartype.BINARY
        elif cppvartype == cppVartype.SPIN:
            return Vartype.SPIN
        elif cppvartype == cppVartype.INTEGER:
            return Vartype.INTEGER
        elif cppvartype == cppVartype.REAL:
            return Vartype.REAL
        else:
            raise RuntimeError("unexpected vartype")
