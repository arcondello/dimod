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

import typing
import warnings

import dimod

from dimod.cyconstrained import cyConstrainedQuadraticModel
from dimod.vartypes import as_vartype


class ConstrainedQuadraticModel(cyConstrainedQuadraticModel):

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
