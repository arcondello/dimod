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

import numbers
import typing
import warnings

import dimod

from dimod.binary import BinaryQuadraticModel
from dimod.cyconstrained import cyConstrainedQuadraticModel
from dimod.quadratic import QuadraticModel
from dimod.sym import Comparison, Sense
from dimod.vartypes import as_vartype, Vartype


class ConstrainedQuadraticModel(cyConstrainedQuadraticModel):

    def _iterable_to_qm(self, iterable: typing.Iterable) -> QuadraticModel:
        qm = QuadraticModel()

        def _add_variable(v):
            # handles vartype, and bounds
            vartype = self.vartype(v)

            if vartype is not Vartype.SPIN and vartype is not Vartype.BINARY:
                # need to worry about bounds
                qm.add_variable(vartype, v,
                                lower_bound=self.lower_bound(v),
                                upper_bound=self.upper_bound(v))
            else:
                qm.add_variable(vartype, v)

        for *variables, bias in iterable:
            if len(variables) == 0:
                qm.offset += bias
            elif len(variables) == 1:
                v, = variables
                _add_variable(v)
                qm.add_linear(v, bias)
            elif len(variables) == 2:
                u, v = variables
                _add_variable(u)
                _add_variable(v)
                qm.add_quadratic(u, v, bias)
            else:
                raise ValueError("terms must be constant, linear or quadratic")
        return qm

    def add_constraint(self, data, *args, **kwargs) -> typing.Hashable:
        # in python 3.8+ we can use singledispatchmethod
        if isinstance(data, (BinaryQuadraticModel, QuadraticModel)):
            return self.add_constraint_from_model(data, *args, **kwargs)
        elif isinstance(data, Comparison):
            return self.add_constraint_from_comparison(data, *args, **kwargs)
        elif isinstance(data, typing.Iterable):
            return self.add_constraint_from_iterable(data, *args, **kwargs)
        else:
            raise TypeError("unexpected data format")

    def add_constraint_from_comparison(self,
                                       comp: Comparison,
                                       *,
                                       label: typing.Optional[typing.Hashable] = None,
                                       copy: bool = True) -> typing.Hashable:
        if not isinstance(comp.rhs, numbers.Number):
            raise TypeError("comparison should have a numeric rhs")

        if isinstance(comp.lhs, (BinaryQuadraticModel, QuadraticModel)):
            return self.add_constraint_from_model(comp.lhs, comp.sense, rhs=comp.rhs,
                                                  label=label, copy=copy)
        else:
            raise ValueError("comparison should have a binary quadratic model "
                             "or quadratic model lhs.")

    def add_constraint_from_iterable(self, iterable: typing.Iterable,
                                     sense: typing.Union[Sense, str],
                                     rhs: dimod.typing.Bias = 0,
                                     *,
                                     label: typing.Optional[typing.Hashable] = None,
                                     ) -> typing.Hashable:
        # we could cythonize this one, but let's keep it simple for now
        qm = self._iterable_to_qm(iterable)

        # use quadratic model in the future
        return self.add_constraint_from_model(
            qm, sense, rhs=rhs, label=label, copy=False)

    def add_constraint_from_model(self,
                                  qm: typing.Union[dimod.BinaryQuadraticModel,
                                                   dimod.QuadraticModel],
                                  sense: typing.Union[dimod.sym.Sense, str],
                                  rhs: dimod.typing.Bias = 0,
                                  *,
                                  label: typing.Optional[typing.Hashable] = None,
                                  copy: bool = True) -> typing.Hashable:

        if qm.dtype == object:
            raise NotImplementedError

        return super().add_constraint_from_model(qm.data, sense, rhs, label=label, copy=copy)

    def add_variable(self,
                     vartype: dimod.typing.VartypeLike,
                     v: typing.Optional[dimod.typing.Variable] = None,
                     *,
                     lower_bound: typing.Optional[float] = None,
                     upper_bound: typing.Optional[float] = None,
                     ) -> dimod.typing.Variable:
        try:
            vartype = as_vartype(vartype, extended=True)
        except TypeError:
            # in dimod<0.11 the argument order was v, vartype so let's allow that case
            warnings.warn(
                "Parameter order CQM.add_variable(v, vartype) "
                "is deprecated since dimod 0.11.0 and will be removed in 0.13.0. "
                "Use CQM.add_variable(vartype, v) instead.",
                DeprecationWarning, stacklevel=2)
            v, vartype = vartype, v

        return super().add_variable(vartype, v, lower_bound=lower_bound, upper_bound=upper_bound)

    def set_objective(self, objective: typing.Union[dimod.BinaryQuadraticModel,
                                                    dimod.QuadraticModel,
                                                    typing.Iterable],
                      ) -> None:
        if isinstance(objective, typing.Iterable):
            raise NotImplementedError

        elif objective.dtype == object:
            # TODO: test and handle
            raise NotImplementedError

        super().set_objective(objective.data)


CQM = ConstrainedQuadraticModel
